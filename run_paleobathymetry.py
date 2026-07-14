#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# simple_paleobathymetry
#
# A streamlined, config-driven workflow to reconstruct paleobathymetry of ocean
# crust (isostatically compensated for sediment load) through geological time.
#
# The pipeline runs four steps:
#   1. convert seafloor age to basement depth        (choice of subsidence model)
#   2. generate distance-to-passive-margin grids     (ocean_basin_proximity)
#   3. predict compacted sediment thickness          (Dutkiewicz et al., 2017)
#   4. compute paleobathymetry with sediments        (Sykes, 1996 isostasy)
#
# Step 1 offers four age -> depth (thermal subsidence) models, selected with
# `age_depth.model` in the config:
#   gdh1             - Stein & Stein (1992)      GDH1 plate-cooling analytic model
#   rhcw18           - Richards et al. (2018)    tabulated plate-cooling model
#   parsons_sclater  - Parsons & Sclater (1977)  analytic model
#   crosby09         - Crosby et al. (2009)      analytic model
#
# Everything is controlled by a single YAML config file:
#
#     python run_paleobathymetry.py config.yml
#
# Copyright (C) 2026 The University of Sydney / EarthByte
# Distributed under the GNU General Public License, version 2 (see LICENSE).
# -----------------------------------------------------------------------------

import argparse
import os
import sys
import math
import warnings

import numpy as np
import yaml


# =============================================================================
# Small helpers
# =============================================================================

def log(msg):
    print("[simple_paleobathymetry] {}".format(msg), flush=True)


def load_config(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return cfg


def target_grid(spacing_degrees):
    """Gridline-registered global lon/lat vectors at the given spacing.

    Matches the registration used by the distance / sediment-thickness grids
    (GMT gridline registration on -R-180/180/-90/90), so all grids align.
    """
    n_lon = int(round(360.0 / spacing_degrees)) + 1
    n_lat = int(round(180.0 / spacing_degrees)) + 1
    lon = np.linspace(-180.0, 180.0, n_lon)
    lat = np.linspace(-90.0, 90.0, n_lat)
    return lon, lat


def _load_zgrid(path):
    """Load a netCDF grid as an xarray DataArray with dims ('lat', 'lon').

    Handles the common coordinate naming conventions (x/y, lon/lat,
    longitude/latitude) used by GMT, GPlately and age-grid products.
    """
    import xarray as xr

    ds = xr.open_dataset(path)
    zname = "z" if "z" in ds.data_vars else list(ds.data_vars)[0]
    da = ds[zname]

    rename = {}
    for d in da.dims:
        dl = str(d).lower()
        if dl in ("x", "lon", "longitude"):
            rename[d] = "lon"
        elif dl in ("y", "lat", "latitude"):
            rename[d] = "lat"
    da = da.rename(rename)

    if "lat" not in da.dims or "lon" not in da.dims:
        raise ValueError(
            "Could not identify lon/lat coordinates in {} (dims={})".format(path, list(da.dims))
        )
    return da.transpose("lat", "lon")


def _write_zgrid(path, z, lon, lat):
    """Write a 2-D array to netCDF as variable 'z' with dims ('lat', 'lon')."""
    import xarray as xr

    ds = xr.Dataset(
        {"z": (("lat", "lon"), np.asarray(z, dtype="float64"))},
        coords={"lat": lat, "lon": lon},
    )
    ds.to_netcdf(path, encoding={"z": {"zlib": True}})
    ds.close()


# =============================================================================
# Plate model (GPlately Plate Model Manager, or local files)
# =============================================================================

def load_plate_model(cfg):
    """Return a dict describing the plate model, from PMM or from local files.

    Keys returned:
        rotation_files, topology_files, coastline_files, cob_files,
        anchor_plate_id, big_time, age_grid(time)->path
    """
    pm = cfg["plate_model"]
    anchor = int(pm.get("anchor_plate_id", 0))
    data_dir = pm.get("data_dir", "plate_model")
    age_cfg = cfg["age_grids"]

    if pm.get("use_pmm", True):
        from plate_model_manager import PlateModelManager

        name = pm["name"]
        pmm = PlateModelManager()
        available = pmm.get_available_model_names()
        if name not in available:
            raise ValueError(
                "Plate model '{}' is not available in the Plate Model Manager.\n"
                "Available models: {}".format(name, ", ".join(sorted(available)))
            )
        model = pmm.get_model(name, data_dir=data_dir)
        log("Using PMM plate model '{}' (valid {}-{} Ma).".format(
            name, model.get_small_time(), model.get_big_time()))

        info = dict(
            rotation_files=model.get_rotation_model(),
            topology_files=model.get_topologies(),
            coastline_files=model.get_coastlines(return_none_if_not_exist=True) or [],
            cob_files=model.get_COBs(return_none_if_not_exist=True) or [],
            anchor_plate_id=anchor,
            big_time=int(model.get_big_time()),
        )

        if str(age_cfg.get("source", "pmm")).lower() == "pmm":
            def age_grid(t, _model=model):
                return _model.get_age_grid(t)
        else:
            template = age_cfg["local_template"]
            def age_grid(t, _tmpl=template):
                return _tmpl.format(time=t)
        info["age_grid"] = age_grid
        return info

    # --- local plate model files -------------------------------------------
    local = pm.get("local", {}) or {}
    info = dict(
        rotation_files=list(local.get("rotation_files", []) or []),
        topology_files=list(local.get("topology_files", []) or []),
        coastline_files=list(local.get("coastline_files", []) or []),
        cob_files=list(local.get("cob_features", []) or []),
        anchor_plate_id=anchor,
        big_time=None,
    )
    if not info["rotation_files"] or not info["topology_files"]:
        raise ValueError(
            "use_pmm is false: you must supply plate_model.local.rotation_files "
            "and plate_model.local.topology_files."
        )
    template = age_cfg.get("local_template")
    if str(age_cfg.get("source", "local")).lower() == "pmm":
        raise ValueError("age_grids.source is 'pmm' but plate_model.use_pmm is false.")
    if not template:
        raise ValueError("age_grids.local_template must be set when using a local plate model.")
    def age_grid(t, _tmpl=template):
        return _tmpl.format(time=t)
    info["age_grid"] = age_grid
    return info


# =============================================================================
# Step 1: seafloor age -> basement depth  (thermal subsidence model)
# =============================================================================
#
# As oceanic lithosphere moves away from a mid-ocean ridge it cools, contracts
# and sinks. "Depth to basement" is the depth of the top of the igneous crust
# (i.e. the seafloor BEFORE any sediment is added on top). Every model below
# turns a seafloor age (Ma) into that depth. Depths are returned in metres and
# are NEGATIVE downwards (e.g. -2600 m = 2600 m below sea level). Cells with no
# age (NaN, i.e. continents / no ocean floor) stay NaN.
#
# Four models are provided so results can be compared; pick one with
# `age_depth.model` in the config.

# Cache of loaded age-depth lookup tables, keyed by file path (so the Richards
# table is only read from disk once even across many time steps).
_AGE_DEPTH_TABLE_CACHE = {}


def _load_age_depth_table(path):
    """Load a two-column (age_Ma, depth_m) age-depth lookup table.

    Whitespace-separated; supports Fortran-style E-notation. Returns
    (ages, depths) as ascending 1-D arrays of positive depths (m).
    """
    import pandas as pd

    path = os.path.abspath(path)
    if path in _AGE_DEPTH_TABLE_CACHE:
        return _AGE_DEPTH_TABLE_CACHE[path]
    if not os.path.isfile(path):
        raise FileNotFoundError(
            "age-depth table not found: {}\n"
            "The 'rhcw18' model needs its lookup table (shipped in data/).".format(path)
        )
    tbl = pd.read_csv(path, sep=r"\s+", header=None, names=["age", "depth"],
                      engine="python")
    tbl = tbl.dropna().sort_values("age")
    ages = tbl["age"].to_numpy(dtype="float64")
    depths = tbl["depth"].to_numpy(dtype="float64")
    _AGE_DEPTH_TABLE_CACHE[path] = (ages, depths)
    return ages, depths


def convert_age_to_depth(age, model="gdh1", richards_table_path=None):
    """Convert seafloor age (Ma) to depth-to-basement (m, negative down).

    Parameters
    ----------
    age : array-like
        Seafloor age in Ma (NaN where there is no ocean floor).
    model : str
        One of:
          'gdh1'            Stein & Stein (1992) GDH1 plate-cooling model.
          'rhcw18'          Richards, Hoggard, Cowton & White (2018) tabulated
                            plate-cooling model, interpolated from a lookup table.
          'parsons_sclater' Parsons & Sclater (1977) analytic model.
          'crosby09'        Crosby, McKenzie & Sclater (2009) analytic model.
        Common alternative spellings are also accepted.
    richards_table_path : str, optional
        Path to the RHCW18 age-depth lookup table (required for 'rhcw18').

    Returns
    -------
    numpy.ndarray
        Depth to basement in metres, negative downwards, NaN where age is NaN.
    """
    age = np.asarray(age, dtype="float64")
    depth = np.full(age.shape, np.nan, dtype="float64")
    valid = ~np.isnan(age)
    key = str(model).strip().lower()

    # ---- GDH1: Stein & Stein (1992), Nature 359, 123-129 --------------------
    if key in ("gdh1", "ghd1", "stein_stein", "steinstein"):
        neg = valid & (age < 0)
        young = valid & (age >= 0) & (age <= 20)
        old = valid & (age > 20)
        depth[neg] = -2600.0                              # zero-age ridge depth
        depth[young] = -(2600.0 + 365.0 * np.sqrt(age[young]))
        depth[old] = -(5651.0 - 2473.0 * np.exp(-0.0278 * age[old]))

    # ---- RHCW18: Richards et al. (2018), JGR Solid Earth 123, 9136-9161 -----
    # Tabulated plate-cooling model (preferred parameters from the RHCW18_Plate
    # _Model repository). We linearly interpolate the age-depth table; ages
    # beyond the table are clamped to its deepest value, and depths are made
    # negative to match the sign convention here.
    elif key in ("rhcw18", "rchw18", "richards", "richards18", "r18"):
        if richards_table_path is None:
            raise ValueError(
                "age_depth.model is 'rhcw18' but no age_depth.richards_table "
                "path was supplied."
            )
        ages_tbl, depths_tbl = _load_age_depth_table(richards_table_path)
        depth[valid] = -1.0 * np.interp(age[valid], ages_tbl, depths_tbl, left=0.0)

    # ---- Parsons & Sclater (1977), JGR 82, 803-827 -------------------------
    elif key in ("parsons_sclater", "ps", "ps_tbl", "parsons"):
        neg = valid & (age < 0)
        young = valid & (age >= 0) & (age < 70)
        old = valid & (age >= 70)
        depth[neg] = -2500.0                              # zero-age ridge depth
        depth[young] = -(2500.0 + 350.0 * np.sqrt(age[young]))
        depth[old] = -(6400.0 - 3200.0 * np.exp(-age[old] / 62.8))

    # ---- Crosby et al. (2006) / Crosby & McKenzie (2009) ------------------
    elif key in ("crosby09", "crosby", "crosby_2009", "cmk09"):
        young = valid & (age >= 0) & (age <= 75)
        mid = valid & (age > 75) & (age <= 160)
        old = valid & (age > 160)
        depth[young] = -(2652.0 + 324.0 * np.sqrt(age[young]))
        depth[mid] = -(5028.0 + 5.26 * age[mid]
                       - 250.0 * np.sin((age[mid] - 75.0) / 30.0))
        depth[old] = -5750.0

    else:
        raise ValueError(
            "Unknown age_depth.model '{}'. Choose one of: "
            "gdh1, rhcw18, parsons_sclater, crosby09.".format(model)
        )

    return depth


def step_age_to_depth(cfg, model, times, out_dir):
    ad = cfg.get("age_depth", {}) or {}
    depth_model = ad.get("model", "gdh1")
    richards_path = ad.get("richards_table")
    log("STEP 1: converting seafloor age to basement depth "
        "(model: {}) ...".format(depth_model))
    os.makedirs(out_dir, exist_ok=True)
    spacing = float(cfg["grids"]["output_spacing"])
    lon, lat = target_grid(spacing)

    for t in times:
        age_path = model["age_grid"](t)
        age_da = _load_zgrid(age_path).interp(lon=lon, lat=lat)
        depth = convert_age_to_depth(age_da.data, model=depth_model,
                                     richards_table_path=richards_path)
        out_path = os.path.join(out_dir, "basement_depth_{:.0f}Ma.nc".format(t))
        _write_zgrid(out_path, depth, lon, lat)
    log("STEP 1 done -> {}".format(out_dir))


# =============================================================================
# Step 2: distance-to-passive-margin grids  (ocean_basin_proximity)
# =============================================================================

def _resolve_proximity_and_obstacles(cfg, model):
    """Return (proximity_files, obstacle_files, plate_boundary_obstacles).

    Default: proximity = COB *line segments* from the plate model,
             obstacles = coastlines (shortest-path around continents).
    Option:  continent contouring -> passive-margin features / contour features.
    """
    prox = cfg["proximity"]

    if prox.get("use_continent_contouring", False):
        cc = prox.get("continent_contouring", {}) or {}
        pm_feat = cc.get("passive_margin_features")
        cc_feat = cc.get("continent_contour_features")
        if not pm_feat:
            raise ValueError(
                "proximity.use_continent_contouring is true but "
                "continent_contouring.passive_margin_features is not set."
            )
        proximity_files = [pm_feat]
        obstacle_files = [cc_feat] if (prox.get("use_continent_obstacles", True) and cc_feat) else None
    else:
        proximity_files = list(model["cob_files"])
        if not proximity_files:
            raise ValueError(
                "No COB line-segment features found for the plate model. "
                "Provide plate_model.local.cob_features, or enable continent contouring."
            )
        obstacle_files = list(model["coastline_files"]) if prox.get("use_continent_obstacles", True) else None
        if prox.get("use_continent_obstacles", True) and not obstacle_files:
            log("WARNING: use_continent_obstacles is true but no coastline files were found; "
                "falling back to straight-line great-circle distances.")
            obstacle_files = None

    pbo = prox.get("plate_boundary_obstacles", ["MidOceanRidge", "SubductionZone"])
    return proximity_files, obstacle_files, pbo


def _report_cob_geometry_types(proximity_files):
    """Log the geometry types of the resolved COB proximity features (info only).

    The distance-to-margin computation measures the distance from each ocean
    point to the *continent-ocean boundary*. That target may be supplied either
    as COB line segments (polylines) or as continental COB polygons -- both are
    valid: for a polygon the distance is measured to its boundary, i.e. the COB
    line. This routine just reports which type was loaded so the run log is
    self-documenting; it never warns or fails. Any error is swallowed.
    """
    try:
        import pygplates

        n_polyline = n_polygon = n_other = 0
        for path in proximity_files:
            try:
                feature_collection = pygplates.FeatureCollection(path)
            except Exception:
                continue
            for feature in feature_collection:
                try:
                    geometries = feature.get_geometries()
                except Exception:
                    try:
                        geometries = feature.get_all_geometries()
                    except Exception:
                        geometries = []
                for geometry in geometries:
                    if isinstance(geometry, pygplates.PolygonOnSphere):
                        n_polygon += 1
                    elif isinstance(geometry, pygplates.PolylineOnSphere):
                        n_polyline += 1
                    else:
                        n_other += 1

        if n_polyline or n_polygon or n_other:
            log("  COB proximity geometries: {} polyline(s), {} polygon(s), "
                "{} other -- distance is measured to the continent-ocean "
                "boundary either way.".format(n_polyline, n_polygon, n_other))
    except Exception:
        # Purely informational; never interrupt the workflow.
        pass


def step_distance_grids(cfg, model, times, out_dir):
    log("STEP 2: generating distance-to-passive-margin grids ...")
    os.makedirs(out_dir, exist_ok=True)

    import pygplates
    from ocean_basin_proximity import (
        generate_and_write_proximity_data_parallel,
        generate_input_points_grid,
    )

    internal = float(cfg["grids"]["internal_spacing"])
    output = float(cfg["grids"]["output_spacing"])
    clamp_km = cfg["proximity"].get("clamp_mean_distance_km", None)

    proximity_files, obstacle_files, pbo = _resolve_proximity_and_obstacles(cfg, model)
    log("  proximity target: {}".format(proximity_files))
    log("  continent obstacles: {}".format(obstacle_files))
    if not cfg["proximity"].get("use_continent_contouring", False):
        _report_cob_geometry_types(proximity_files)

    input_points = generate_input_points_grid(internal)[0]

    # Pre-resolve the age grid path for every time.
    age_grid_filenames_and_paleo_times = [(model["age_grid"](t), t) for t in times]

    max_recon = cfg["time"].get("max_reconstruction_time", None)
    if max_recon is None:
        max_recon = model.get("big_time", None)

    clamp_radians = None
    if clamp_km:
        clamp_radians = float(clamp_km) / pygplates.Earth.mean_radius_in_kms

    kwargs = dict(
        input_points=input_points,
        rotation_filenames=model["rotation_files"],
        proximity_filenames=proximity_files,
        proximity_features_are_topological=False,
        proximity_feature_types=None,
        topological_reconstruction_filenames=model["topology_files"],
        age_grid_filenames_and_paleo_times=age_grid_filenames_and_paleo_times,
        time_increment=1,
        output_distance_with_time=False,
        output_mean_distance=True,
        output_standard_deviation_distance=False,
        output_directory=out_dir,
        max_topological_reconstruction_time=max_recon,
        anchor_plate_id=model["anchor_plate_id"],
        proximity_distance_threshold_radians=None,
        clamp_mean_proximity_distance_radians=clamp_radians,
        output_grd_files=(internal, output),
        num_cpus=int(cfg["run"].get("num_cpus", 1)),
    )
    if obstacle_files:
        kwargs["continent_obstacle_filenames"] = obstacle_files
        kwargs["plate_boundary_obstacle_feature_types"] = pbo

    try:
        generate_and_write_proximity_data_parallel(**kwargs)
    except TypeError as exc:
        # Older ocean_basin_proximity.py hard-wires the plate-boundary obstacles
        # (mid-ocean ridges + subduction zones) and does not accept this kwarg.
        if "plate_boundary_obstacle_feature_types" in str(exc):
            kwargs.pop("plate_boundary_obstacle_feature_types", None)
            generate_and_write_proximity_data_parallel(**kwargs)
        else:
            raise
    log("STEP 2 done -> {}".format(out_dir))


# =============================================================================
# Step 3: predicted compacted sediment thickness  (Dutkiewicz et al., 2017)
# =============================================================================

def step_sediment_thickness(cfg, model, times, distance_dir, out_dir):
    log("STEP 3: predicting compacted sediment thickness ...")
    os.makedirs(out_dir, exist_ok=True)

    import joblib
    from ocean_basin_proximity import generate_input_points_grid
    from predict_sediment_thickness import predict_sedimentation, write_grd_file

    st = cfg["sediment_thickness"]
    output = float(cfg["grids"]["output_spacing"])
    input_points = generate_input_points_grid(output)[0]

    distance_template = os.path.join(
        distance_dir, "mean_distance_{:.1f}d".format(output) + "_{:.1f}.nc"
    )

    def _one(t):
        result = predict_sedimentation(
            input_points=input_points,
            age_grid_filename=model["age_grid"](t),
            distance_grid_filename=distance_template.format(t),
            mean_age=st["mean_age"],
            mean_distance=st["mean_distance"],
            variance_age=st["variance_age"],
            variance_distance=st["variance_distance"],
            age_distance_polynomial_coefficients=st["polynomial_coefficients"],
            max_age=st.get("max_age"),
            max_distance=st.get("max_distance"),
        )
        write_grd_file(
            os.path.join(out_dir, "sediment_thickness_{:.0f}Ma.nc".format(t)),
            output_data=result,
            grid_spacing=output,
            num_grid_longitudes=None,
            num_grid_latitudes=None,
        )

    n_jobs = int(cfg["run"].get("num_cpus", 1))
    with joblib.Parallel(n_jobs=n_jobs) as parallel:
        parallel(joblib.delayed(_one)(t) for t in times)
    log("STEP 3 done -> {}".format(out_dir))


# =============================================================================
# Step 4: paleobathymetry with sediments (isostatic)  (Sykes, 1996)
# =============================================================================

def step_paleobathymetry(cfg, times, sedthick_dir, basement_dir, out_dir):
    log("STEP 4: computing isostatically-compensated paleobathymetry ...")
    os.makedirs(out_dir, exist_ok=True)

    for t in times:
        sed = _load_zgrid(os.path.join(sedthick_dir, "sediment_thickness_{:.0f}Ma.nc".format(t)))
        base = _load_zgrid(os.path.join(basement_dir, "basement_depth_{:.0f}Ma.nc".format(t)))
        # Align basement to the sediment-thickness grid (guards against any
        # registration differences between products).
        base = base.interp(lon=sed.lon, lat=sed.lat)

        sed_m = sed.data                       # sediment thickness (m, positive)
        basement_m = base.data                 # basement depth (m, negative down)

        # Isostatic sediment-load correction (Sykes, 1996), thickness in km:
        #   correction(m) = (0.43422 * h_km - 0.010395 * h_km^2) * 1000, >= 0
        h_km = sed_m / 1000.0
        iso = (0.43422 * h_km - 0.010395 * (h_km ** 2)) * 1000.0
        iso = np.where(iso >= 0.0, iso, 0.0)

        # Seafloor depth including the (isostatically compensated) sediment pile.
        paleobath = basement_m + sed_m - iso

        out_path = os.path.join(out_dir, "paleobathymetry_{:.0f}Ma.nc".format(t))
        _write_zgrid(out_path, paleobath, sed.lon.data, sed.lat.data)
    log("STEP 4 done -> {}".format(out_dir))


# =============================================================================
# Driver
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Streamlined paleobathymetry workflow (age -> depth -> "
                    "distance -> sediment thickness -> paleobathymetry).")
    parser.add_argument("config", help="Path to the YAML configuration file.")
    args = parser.parse_args()

    cfg = load_config(args.config)

    # Make the public predicting-sediment-thickness package importable.
    pst_dir = cfg["run"].get("predicting_sediment_thickness_dir",
                             "submodules/predicting-sediment-thickness")
    pst_dir = os.path.abspath(pst_dir)
    if not os.path.isdir(pst_dir):
        log("ERROR: predicting-sediment-thickness not found at '{}'.".format(pst_dir))
        log("Clone it, e.g.:")
        log("  git submodule add https://github.com/EarthByte/predicting-sediment-thickness "
            "submodules/predicting-sediment-thickness")
        sys.exit(1)
    if pst_dir not in sys.path:
        sys.path.insert(0, pst_dir)

    # Feature toggles that this streamlined workflow deliberately does not implement.
    feats = cfg.get("features", {}) or {}
    if feats.get("use_lips") or feats.get("use_seamounts"):
        log("NOTE: LIP / seamount contributions are not included in this "
            "streamlined workflow; proceeding without them.")

    t = cfg["time"]
    times = list(range(int(t["min"]), int(t["max"]) + 1, int(t.get("step", 1))))
    log("Time range: {} to {} Ma (step {}), {} times.".format(
        t["min"], t["max"], t.get("step", 1), len(times)))

    out_root = cfg["run"].get("output_dir", "output")
    basement_dir = os.path.join(out_root, "BasementDepth")
    distance_dir = os.path.join(out_root, "Distances")
    sedthick_dir = os.path.join(out_root, "SedimentThickness")
    paleobath_dir = os.path.join(out_root, "Paleobathymetry")
    os.makedirs(out_root, exist_ok=True)

    steps = cfg["run"].get("steps", {}) or {}
    model = load_plate_model(cfg)

    if steps.get("age_to_depth", True):
        step_age_to_depth(cfg, model, times, basement_dir)
    if steps.get("distance_grids", True):
        step_distance_grids(cfg, model, times, distance_dir)
    if steps.get("sediment_thickness", True):
        step_sediment_thickness(cfg, model, times, distance_dir, sedthick_dir)
    if steps.get("paleobathymetry", True):
        step_paleobathymetry(cfg, times, sedthick_dir, basement_dir, paleobath_dir)

    log("All requested steps complete. Outputs in '{}'.".format(out_root))


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        main()
