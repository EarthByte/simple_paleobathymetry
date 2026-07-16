#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# generate_continent_contours.py
#
# OPTIONAL step of the paleobathymetry workflow: generate *dynamically
# contoured* passive margins through time, instead of using a static COB
# line-segment file. This is the in-workflow equivalent of the EarthByte
# `continent-contouring` workflow (`create_passive_margins.py`), built on the
# same engine, which ships inside gplately: gplately.ptt.continent_contours.
#
# It can be run in two ways:
#
#   1. Standalone, as a pre-step:
#          python generate_continent_contours.py config.yml
#
#   2. Automatically, from run_paleobathymetry.py, when the config has
#          proximity:
#            use_continent_contouring: true
#            continent_contouring:
#              generate:
#                enabled: true
#                continent_polygon_files: [ ... ]
#      In that case run_paleobathymetry.py imports generate_contours() and runs
#      it before the distance step, then feeds the results straight in.
#
# For each requested time it produces, in the contouring output directory:
#     continent_contour_features_<t>.gpml   (COB polylines  -> distance obstacles)
#     passive_margin_features_<t>.gpml      (contours minus active/subduction
#                                            segments -> proximity target)
#     continent_mask_<t>.nc                 (continental-crust mask)
# and aggregated collections used by the workflow:
#     continent_contour_features.gpmlz
#     passive_margin_features.gpmlz
#
# Requires: pygplates + gplately (already workflow dependencies).
#
# NOTE: this reproduces the EarthByte create_passive_margins.py logic.
# -----------------------------------------------------------------------------

import os
import argparse

import yaml


def _load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def _feature_collection(paths):
    import pygplates
    fc = pygplates.FeatureCollection()
    for p in paths:
        fc.add(pygplates.FeatureCollection(p))
    return fc


def _make_buffer_and_gap_ramp(ramp_cfg, earth_r):
    """Deep-time buffer/gap ramp (opt-in), ported from the commented-out
    continent_contouring_buffer_and_gap_distance_radians() function in
    EarthByte's create_passive_margins.py: a larger buffer/gap distance before
    Pangea assembly, linearly interpolated down to a smaller one after, then
    linearly reduced for continents smaller than an area threshold.
    """
    import math

    pre_rad = math.radians(float(ramp_cfg.get("pre_pangea_distance_degrees", 2.5)))
    post_rad = math.radians(float(ramp_cfg.get("post_pangea_distance_degrees", 0.0)))
    pre_time = float(ramp_cfg.get("pre_pangea_time_ma", 300))
    post_time = float(ramp_cfg.get("post_pangea_time_ma", 250))
    area_threshold_sr = (float(ramp_cfg.get("small_continent_area_threshold_square_kms", 500000))
                          / (earth_r * earth_r))

    def buffer_and_gap(time, contoured_continent):
        if time > pre_time:
            distance_rad = pre_rad
        elif time < post_time:
            distance_rad = post_rad
        else:
            interp = float(time - post_time) / (pre_time - post_time)
            distance_rad = interp * pre_rad + (1 - interp) * post_rad

        area_sr = contoured_continent.get_area()
        if area_sr < area_threshold_sr:
            distance_rad *= area_sr / area_threshold_sr
        return distance_rad

    return buffer_and_gap


def _contour_output_dir(cfg):
    """Directory the contour/passive-margin/mask files are written to."""
    prox = cfg["proximity"]
    cc_cfg = (prox.get("continent_contouring", {}) or {})
    gen = (cc_cfg.get("generate", {}) or {})
    out_root = cfg.get("run", {}).get("output_dir", "output")
    return (cc_cfg.get("output_dir")
            or gen.get("output_dir")
            or os.path.join(out_root, "ContinentContours"))


def aggregated_feature_paths(cfg):
    """(passive_margin_features.gpmlz, continent_contour_features.gpmlz) paths."""
    out_dir = _contour_output_dir(cfg)
    return (os.path.join(out_dir, "passive_margin_features.gpmlz"),
            os.path.join(out_dir, "continent_contour_features.gpmlz"))


def _resolve_plate_inputs(cfg):
    """Resolve (rotation_files, topology_files, continent_files, anchor) for
    contouring, from either a PMM plate model or local plate-model files --
    mirrors run_paleobathymetry.load_plate_model().

    Continental polygons: proximity.continent_contouring.generate.continent_polygon_files
    overrides when given; otherwise, for a PMM model, they are auto-sourced from
    model.get_continental_polygons(). A local model always requires them explicitly.
    """
    prox = cfg["proximity"]
    cc_cfg = (prox.get("continent_contouring", {}) or {})
    gen = (cc_cfg.get("generate", {}) or {})

    pm = cfg["plate_model"]
    anchor = int(pm.get("anchor_plate_id", 0))
    continent_files = list(gen.get("continent_polygon_files", []) or [])

    if pm.get("use_pmm", True):
        from plate_model_manager import PlateModelManager

        name = pm["name"]
        pmm = PlateModelManager()
        model = pmm.get_model(name, data_dir=pm.get("data_dir", "plate_model"))
        rotation_files = list(model.get_rotation_model() or [])
        topology_files = list(model.get_topologies() or [])
        if not continent_files:
            continent_files = list(
                model.get_continental_polygons(return_none_if_not_exist=True) or [])
        if not continent_files:
            raise ValueError(
                "continent contouring needs continental polygons, but PMM model "
                "'{}' does not provide any and "
                "proximity.continent_contouring.generate.continent_polygon_files "
                "is empty. Supply continent_polygon_files explicitly.".format(name))
        local_paths = []  # PMM paths are already resolved/cached; skip existence check below.
    else:
        local = pm.get("local", {}) or {}
        rotation_files = list(local.get("rotation_files", []) or [])
        topology_files = list(local.get("topology_files", []) or [])
        if not rotation_files:
            raise ValueError("continent contouring needs plate_model.local.rotation_files.")
        if not continent_files:
            raise ValueError(
                "continent contouring needs "
                "proximity.continent_contouring.generate.continent_polygon_files "
                "(the continental polygons to contour).")
        local_paths = rotation_files + continent_files + topology_files

    missing = [p for p in local_paths if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError("continent-contouring input(s) not found: {}".format(missing))

    return rotation_files, topology_files, continent_files, anchor


def generate_contours(cfg, log=print):
    """Dynamically contour continents through time; write passive-margin /
    contour / mask files. Returns (passive_margin_path, contour_features_path).

    `cfg` is the already-parsed workflow config dict. Uses:
        plate_model.{use_pmm, name, data_dir, anchor_plate_id, local.rotation_files,
        local.topology_files}
        proximity.continent_contouring.generate.*   (parameters below)
        time.{min,max,step}
        run.output_dir
    """
    import pygplates
    import gplately
    from gplately.ptt import continent_contours

    prox = cfg["proximity"]
    cc_cfg = (prox.get("continent_contouring", {}) or {})
    gen = (cc_cfg.get("generate", {}) or {})

    # --- plate-model inputs -------------------------------------------------
    rotation_files, topology_files, continent_files, anchor = _resolve_plate_inputs(cfg)

    # --- time range ---------------------------------------------------------
    t = cfg["time"]
    start, end, step = int(t["min"]), int(t["max"]), int(t.get("step", 1))
    times = list(range(start, end + 1, step))

    # --- contouring parameters (with sensible defaults) --------------------
    earth_r = pygplates.Earth.mean_radius_in_kms
    spacing = float(gen.get("point_spacing_degrees", 0.25))
    area_km2 = gen.get("area_threshold_square_kms", None)
    if area_km2 is not None:
        area_threshold_sr = float(area_km2) / (earth_r * earth_r)
    else:
        area_threshold_sr = float(gen.get("area_threshold_steradians", 0.0))
    buffer_and_gap_mode = str(gen.get("buffer_and_gap_mode", "constant")).lower()
    if buffer_and_gap_mode == "ramp":
        buffer_gap_rad = _make_buffer_and_gap_ramp(gen.get("buffer_and_gap_ramp", {}) or {}, earth_r)
    else:
        buffer_gap_km = gen.get("buffer_and_gap_distance_kms", 0.0)
        buffer_gap_rad = float(buffer_gap_km) / earth_r if buffer_gap_km else 0.0
    excl_km2 = gen.get("exclusion_area_threshold_square_kms", None)
    if excl_km2 is not None:
        exclusion_area_sr = float(excl_km2) / (earth_r * earth_r)
    else:
        exclusion_area_sr = float(gen.get("exclusion_area_threshold_steradians", 800000.0 / (earth_r * earth_r)))
    separation_rad = float(gen.get("separation_distance_threshold_radians",
                                    continent_contours.DEFAULT_CONTINENT_SEPARATION_DISTANCE_THRESHOLD_RADIANS))
    max_sub_km = float(gen.get("max_distance_of_subduction_from_active_margin_kms", 500.0))
    max_sub_rad = max_sub_km / earth_r

    out_dir = _contour_output_dir(cfg)
    os.makedirs(out_dir, exist_ok=True)

    rotation_model = pygplates.RotationModel(rotation_files, default_anchor_plate_id=anchor)
    continent_features = _feature_collection(continent_files)
    topology_features = _feature_collection(topology_files) if topology_files else None

    log("[contours] contouring {} times ({}-{} Ma), spacing {} deg -> {}".format(
        len(times), start, end, spacing, out_dir))

    # --- the contouring engine (in gplately) --------------------------------
    contourer = continent_contours.ContinentContouring(
        rotation_model,
        continent_features,
        continent_contouring_point_spacing_degrees=spacing,
        continent_contouring_area_threshold_steradians=area_threshold_sr,
        continent_contouring_buffer_and_gap_distance_radians=buffer_gap_rad,
        continent_exclusion_area_threshold_steradians=exclusion_area_sr,
        continent_separation_distance_threshold_radians=separation_rad,
    )

    all_contours, all_passive = [], []
    for time in times:
        # subduction-zone lines at this time (to mark active margins)
        subduction_lines = []
        if topology_features is not None:
            resolved, shared = [], []
            pygplates.resolve_topologies(topology_features, rotation_model,
                                         resolved, float(time), shared)
            for sbs in shared:
                if (sbs.get_feature().get_feature_type()
                        == pygplates.FeatureType.gpml_subduction_zone):
                    for seg in sbs.get_shared_sub_segments():
                        subduction_lines.append(seg.get_resolved_geometry())

        continent_mask, contoured_continents = \
            contourer.get_continent_mask_and_contoured_continents(float(time))

        contour_feats, passive_feats = [], []
        for continent in contoured_continents:
            for contour in continent.get_contours():
                cf = pygplates.Feature()
                cf.set_geometry(contour)
                cf.set_valid_time(time + 0.5 * step, time - 0.5 * step)
                contour_feats.append(cf)

                # Passive margin = contour minus segments near a subduction zone.
                for margin in _passive_margin_polylines(contour, subduction_lines,
                                                        max_sub_rad, pygplates):
                    pf = pygplates.Feature()
                    pf.set_geometry(margin)
                    pf.set_valid_time(time + 0.5 * step, time - 0.5 * step)
                    passive_feats.append(pf)

        pygplates.FeatureCollection(contour_feats).write(
            os.path.join(out_dir, "continent_contour_features_{}.gpml".format(time)))
        pygplates.FeatureCollection(passive_feats).write(
            os.path.join(out_dir, "passive_margin_features_{}.gpml".format(time)))
        try:
            gplately.grids.write_netcdf_grid(
                os.path.join(out_dir, "continent_mask_{}.nc".format(time)),
                continent_mask.astype("float"))
        except Exception as exc:
            log("[contours]  (continent mask nc for {} Ma skipped: {})".format(time, exc))

        all_contours.extend(contour_feats)
        all_passive.extend(passive_feats)
        log("[contours]  {:>4} Ma: {} contour + {} passive-margin features".format(
            time, len(contour_feats), len(passive_feats)))

    passive_path, contour_path = aggregated_feature_paths(cfg)
    pygplates.FeatureCollection(all_contours).write(contour_path)
    pygplates.FeatureCollection(all_passive).write(passive_path)
    log("[contours] done. Aggregated:")
    log("  {}".format(passive_path))
    log("  {}".format(contour_path))
    return passive_path, contour_path


def _passive_margin_polylines(contour_polyline, subduction_lines, max_distance_radians, pygplates):
    """Split a continent-contour polyline into passive-margin polylines.

    A great-circle-arc segment of the contour is 'active' if it lies within
    max_distance_radians of any subduction-zone line; otherwise it is a passive
    margin. Consecutive passive arcs are joined into polylines. (Mirrors the
    EarthByte create_passive_margins.py logic.)
    """
    points = list(contour_polyline.get_points())
    if len(points) < 2:
        return []

    def _near_subduction(p0, p1):
        if not subduction_lines:
            return False
        seg = pygplates.PolylineOnSphere([p0, p1])
        for sub in subduction_lines:
            if pygplates.GeometryOnSphere.distance(
                    seg, sub, max_distance_radians) is not None:
                return True
        return False

    margins = []
    run = [points[0]]
    for i in range(1, len(points)):
        if _near_subduction(points[i - 1], points[i]):
            if len(run) >= 2:
                margins.append(pygplates.PolylineOnSphere(run))
            run = [points[i]]
        else:
            run.append(points[i])
    if len(run) >= 2:
        margins.append(pygplates.PolylineOnSphere(run))
    return margins


def main():
    ap = argparse.ArgumentParser(
        description="Generate dynamically-contoured passive margins through time.")
    ap.add_argument("config", nargs="?", default="config.yml")
    args = ap.parse_args()
    cfg = _load_config(args.config)
    generate_contours(cfg)


if __name__ == "__main__":
    main()
