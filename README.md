# simple_paleobathymetry

A streamlined, config-driven workflow that reconstructs **paleobathymetry of
ocean crust** through geological time.

Given a plate model and a set of seafloor-age grids, the workflow produces one
paleobathymetry grid per time step, from the present day back through time. It
is deliberately self-contained: everything you need to understand *what it does*
and *how each number is computed* is in this README, so you should not need to
read any other repository to use it.

**Scope.** Steps 1–4 reconstruct the depth of **oceanic crust** only — crust
whose age is defined in the seafloor-age grids. They do **not** compute
paleobathymetry for submerged **continental** crust (continental shelves,
plateaus, rifted or stretched continental margins); those areas have no seafloor
age and are left undefined (`NaN`) in the Step 4 output. An optional **Step 5**
closes this gap using **[pyBacktrack](https://github.com/EarthByte/pyBacktrack)**,
which backtracks/backstrips drill-site and grid stratigraphy to reconstruct
paleo-water-depths on submerged continental crust and on ocean crust preserved
today, merging in the Step 4 grids to also cover crust that has since subducted.

A typical run is a single command:

```bash
python run_paleobathymetry.py config.yml
```

---

## Contents

- [What the workflow computes](#what-the-workflow-computes)
- [Quick start](#quick-start)
- [The five steps in detail](#the-five-steps-in-detail)
  - [Step 1 — Seafloor age → basement depth](#step-1--seafloor-age--basement-depth)
  - [Step 2 — Distance to passive continental margins](#step-2--distance-to-passive-continental-margins)
  - [Step 3 — Predicted sediment thickness](#step-3--predicted-sediment-thickness)
  - [Step 4 — Paleobathymetry with sediment isostasy](#step-4--paleobathymetry-with-sediment-isostasy)
  - [Step 5 — pyBacktrack paleobathymetry (merged)](#step-5--pybacktrack-paleobathymetry-merged)
- [Configuration reference](#configuration-reference)
- [Inputs and outputs](#inputs-and-outputs)
- [Installation](#installation)
- [Frequently asked questions](#frequently-asked-questions)
- [References](#references)
- [License](#license)

---

## What the workflow computes

The depth of the seafloor in the past is not something we can measure directly —
the ocean floor that existed 100 million years ago has largely been subducted.
Instead we *reconstruct* it from physical principles:

1. **Oceanic lithosphere subsides as it ages.** New crust forms hot at
   mid-ocean ridges. As the plate spreads away it loses heat to the ocean, and
   the mantle lithosphere beneath the crust **thickens and cools**, becoming
   denser. (The crust itself stays essentially unchanged in density, apart from
   minor hydrothermal alteration — it is the growth and cooling of the mantle
   lithosphere that matters.) This progressively denser, thicker lithospheric
   column sinks isostatically, so the seafloor deepens with age. That age–depth
   relationship is well established, so a map of seafloor age can be turned into
   a map of basement (bare-rock) depth.

2. **Sediment piles up on top of that basement**, and the seafloor we care
   about is the *top of the sediment*, not the bare rock. Older, more distal
   ocean floor accumulates thicker sediment. The workflow predicts this
   thickness and adds it on top of the basement, applying an isostatic
   correction because the weight of the sediment pushes the crust down further.

The four steps below carry out exactly this reasoning. The final product,
**paleobathymetry**, is the depth of the sediment-covered **ocean floor**:

```
paleobathymetry  =  basement depth  +  sediment thickness  −  isostatic correction
```

All grids are global netCDF files (one per reconstruction time), on a regular
longitude/latitude grid, with depths in **metres, negative downwards**.

---

## Quick start

```bash
# 1. create the software environment (one time)
conda env create -f environment.yml
conda activate simple_paleobathymetry

# 2. add the sediment-thickness engine as a submodule (one time)
git submodule add https://github.com/EarthByte/predicting-sediment-thickness \
    submodules/predicting-sediment-thickness
git submodule update --init --recursive

# 3. run
python run_paleobathymetry.py config.yml
```

The shipped `config.yml` runs the default configuration end to end: the
**Zahirovic2022** plate model (fetched automatically), the **GDH1** age-depth
model, distances measured to the shipped **COB line segments**, and Dutkiewicz et al. (2017)
sediment thickness — no large igneous provinces or seamounts. Edit `config.yml`
to change any of this; every option is explained inline in that file and in the
[configuration reference](#configuration-reference) below.

---

## The five steps in detail

The steps run in order; each writes grids that the next step reads. You can turn
any step on or off in the config (`run.steps`) to re-run part of the pipeline
without recomputing everything. Step 5 is optional and reads Step 4's output
without modifying it.

### Step 1 — Seafloor age → basement depth

**Input:** a seafloor-age grid for each time step (age of the ocean floor in Ma).
**Output:** `output/BasementDepth/basement_depth_<t>Ma.nc` — depth to igneous
basement (metres, negative down), i.e. the seafloor depth *before* any sediment.

Oceanic lithosphere subsides as it cools with age. This step applies a **thermal
subsidence model** to every ocean cell to turn its age into a depth. Cells with
no age (continents, or where the age grid is undefined) stay `NaN`.

Four published models are available; select one with `age_depth.model`. They
agree closely for young crust and differ mainly for old crust (where different
assumptions about the deep thermal structure matter). Providing all four lets
you test the sensitivity of your results to that choice.

**`gdh1` — Stein & Stein (1992), GDH1** *(default)*. The most widely used
plate-cooling reference model:

```
age ≤ 20 Ma :  depth = 2600 + 365·√age
age > 20 Ma :  depth = 5651 − 2473·exp(−0.0278·age)          (metres below sea level)
```

**`rhcw18` — Richards, Hoggard, Cowton & White (2018).** A more recent
plate-cooling model whose preferred parameters do not reduce to a simple
formula, so it is supplied as a precomputed **age-depth lookup table**
(`data/RHCW18_age_depth.dat`, two columns: age in Ma, depth in m). The workflow
linearly interpolates that table for each age; ages older than the table's range
are held at its deepest value. The table ships with this repository, so the
model works out of the box.

**`parsons_sclater` — Parsons & Sclater (1977).** A classic analytic model:

```
age <  70 Ma :  depth = 2500 + 350·√age
age ≥ 70 Ma :  depth = 6400 − 3200·exp(−age / 62.8)
```

**`crosby09` — Crosby, McKenzie & Sclater (2009).** Adds a long-wavelength
correction for old seafloor:

```
0 ≤ age ≤ 75 Ma :  depth = 2652 + 324·√age
75 < age ≤ 160  :  depth = 5028 + 5.26·age − 250·sin((age − 75)/30)
age > 160 Ma    :  depth = 5750
```

In every case the workflow stores depth as a **negative** number (below sea
level).

### Step 2 — Distance to passive continental margins

**Input:** the plate model (rotations + topologies), the seafloor-age grids, and
a clean set of **passive-margin COB line segments** (the default dataset shipped
in `data/`).
**Output:** `output/Distances/mean_distance_<spacing>d_<t>.nc` — for each ocean
cell, its **mean distance (km) to the nearest passive continental margin,
averaged over the lifetime of that piece of ocean floor**.

This is the least obvious step, so it is worth explaining carefully, because the
sediment-thickness prediction in Step 3 depends on it.

**Why distance to a margin?** Most ocean sediment is supplied from the
continents (rivers, shelves, turbidity currents) and rains down from productive
surface waters near land. A patch of seafloor that has always sat far out in the
open ocean receives little sediment; one that has spent its life close to a
continental margin receives a lot. So "distance to the nearest passive margin"
is a strong predictor of how much sediment a piece of seafloor has accumulated.

**Why the *mean over the lifetime*, not just the present-day distance?**
Sediment accumulates over the *entire history* of a piece of ocean floor, not
just where it happens to be now. A parcel of crust may form near a margin, then
drift into the open ocean as the basin widens. What controls its total sediment
load is where it has been *throughout its life*. The workflow therefore:

1. takes each ocean cell that exists at reconstruction time *t*;
2. uses the topological plate model to **reconstruct that cell backwards in
   time**, step by step (1 Myr increments), from *t* to the moment the crust
   formed at a ridge (its own age);
3. at each step, measures the distance from the reconstructed position to the
   nearest passive margin *as it was at that time*;
4. **averages those distances over the cell's whole lifetime**.

That lifetime-averaged distance is what gets written out and fed to Step 3.

**Shortest path *around* continents.** Sediment cannot travel through
continents, so a straight great-circle line to the nearest margin can be
misleading when land lies in between (for example across a narrow isthmus). By
default the workflow therefore measures the **shortest path that routes around
continental obstacles** (present-day coastlines are used as the obstacles),
approximating how sediment would actually be delivered by water. You can switch
this off (`proximity.use_continent_obstacles: false`) to use plain great-circle
distances instead.

**The 0–3000 km range.** The sediment-thickness relationship in Step 3 was
trained on real ocean data, where essentially all ocean floor lies within a few
thousand kilometres of a margin. The relationship is only meaningful within that
calibrated range, so distances are **clamped to a maximum of 3000 km**
(`proximity.clamp_mean_distance_km`). Beyond ~3000 km the prediction would be an
extrapolation, so the clamp both keeps the input in range and reflects that
"very far from any margin" saturates to a thin, pelagic sediment cover.

**What defines a "passive margin" — and why COB *line segments*, not polygons.**
The distance calculation needs a clean set of **continent-ocean boundary (COB)
line segments (polylines) that trace passive margins only**. This matters:

- A **COB polygon** outlines the *entire* continent, so its boundary also runs
  along **active** margins (subduction zones, transforms). Measuring distance to
  a polygon would therefore report short "passive-margin" distances next to
  active margins that receive little continent-derived sediment — **false
  positives** that corrupt the sediment-thickness prediction. **COB polygons
  must not be used** as the proximity target.
- A curated COB **line-segment** dataset includes only the passive-margin
  boundaries and omits active margins, which is exactly what the algorithm
  needs.

Because the GPlately **Plate Model Manager does not deliver COB line segments**
to end users, this workflow **ships its own default dataset** and does *not* take
COBs from the plate model. The shipped file is the global present-day COB line
segments of **Müller et al. (2016)** (Gee & Kent 2007 timescale) — the GPlates
default set, made compatible with Zahirovic2022 — containing passive-margin
boundaries only:

```
data/Global_EarthByte_GPlates_PresentDay_COBs.gpmlz
```

This default is suitable for **most EarthByte Mesozoic–Cenozoic plate models
published up to 2022**. The exception is models with severely altered COB
outlines (the **Clennett / Alfonso** family); for those, supply the matching COB
line segments via `proximity.cob_line_segments`. To use your own COB line
segments, point that setting at your file (it stays independent of the plate
model). As a safeguard, the workflow inspects the proximity features at run time
and **warns if it finds polygons** rather than line segments.

> **⚠️ Completeness of paleo-COBs — read this before choosing a method.**
> The distance-to-passive-margin grid (and therefore the predicted sediment
> thickness that depends on it) is only as complete as the set of passive
> margins it measures distance to. A *static* COB line-segment file is a **fixed
> list of margins**: any passive margin that is missing from that file — or that
> did not yet exist at present day but existed in the past — is simply never seen
> by the algorithm, so ocean floor near it is assigned a distance to some *other,
> farther* margin and its sediment thickness comes out wrong. Static COB files
> routinely omit paleo-margins (rifts that later closed, margins consumed by
> subduction, or margins the compiler simply did not digitise).
>
> **Continent contouring is the only method that guarantees every paleo-COB is
> found and used**, because it *re-derives* the passive margins from the
> reconstructed continents at each time step rather than reading them from a
> fixed list — so a margin that existed at 120 Ma but not today is still
> contoured and used at 120 Ma. **If you do not use continent contouring, it is
> your responsibility to supply a `cob_line_segments` file that already contains
> every relevant paleo-COB for your time range;** the shipped present-day dataset
> does not, beyond the Mesozoic–Cenozoic. When in doubt, switch continent
> contouring on.

Continent contouring is **off by default**, so no extra preprocessing is needed.

**Optional: dynamically contoured passive margins.** Instead of the static COB
line segments, you can have the workflow trace passive margins *at each time
step* by contouring the reconstructed continental crust. This is the in-workflow
equivalent of EarthByte's separate
[continent-contouring](https://github.com/EarthByte/continent-contouring)
workflow — in fact the contouring engine it uses **is gplately's**
(`gplately.ptt.continent_contours`), which you already have installed. The
workflow drives that engine and then splits each continent contour into
**passive** vs **active** margin segments (a segment near a subduction zone is
active; everything else is passive), exactly as EarthByte's
`create_passive_margins.py` does. To switch it on:

```yaml
proximity:
  use_continent_contouring: true
  continent_contouring:
    generate:
      enabled: true
      continent_polygon_files: [plate_model/continental_polygons.gpml]
      # point_spacing_degrees, area_threshold_square_kms, buffer_and_gap_distance_kms,
      # max_distance_of_subduction_from_active_margin_kms ... (see config.yml)
```

When enabled, generation runs **automatically before the distance step** (via
`generate_continent_contours.py`) and writes, into
`continent_contouring.output_dir` (default `output/ContinentContours/`):

```
continent_contour_features_<t>.gpml   passive_margin_features_<t>.gpml   continent_mask_<t>.nc
continent_contour_features.gpmlz      passive_margin_features.gpmlz   (aggregated, used by Step 2)
```

The aggregated `passive_margin_features.gpmlz` becomes the proximity target and
`continent_contour_features.gpmlz` the continent obstacle — the two
`continent_contouring` paths are filled in for you. Existing outputs are reused
unless you set `generate.force: true`. You can also run it as a standalone
pre-step (`python generate_continent_contours.py config.yml`), or skip generation
entirely and point `passive_margin_features` / `continent_contour_features` at
files you already have (leave `generate.enabled: false`).

> **Deep time (beyond ~250 Ma).** The shipped COB line-segment dataset is a
> *present-day* set reconstructed with the plate model, and it represents
> **Mesozoic–Cenozoic** passive margins. For older reconstructions (pre-Pangea
> breakup) it no longer captures the passive margins that existed then, so
> distances would be measured to the wrong boundaries. For deep-time
> paleobathymetry you should **switch continent contouring on** (as above): it
> identifies passive-margin segments *dynamically at each time step* and produces
> a suitable time-dependent proximity file. The workflow prints a reminder if you
> request times beyond ~250 Ma with contouring still off.

### Step 3 — Predicted sediment thickness

**Input:** the seafloor-age grids (Step 1's input) and the distance grids
(Step 2's output).
**Output:** `output/SedimentThickness/sediment_thickness_<t>Ma.nc` — predicted
**compacted sediment thickness (metres)** for each ocean cell at each time.

Here is exactly how the number is produced — you do not need to consult the
sediment-thickness repository.

**The idea.** Dutkiewicz et al. (2017) showed that the thickness of sediment on
ocean crust can be predicted well from just two variables:

- **seafloor age** — older crust has had more time to accumulate sediment, and
- **distance to the nearest passive margin** — nearer crust receives more
  sediment (this is the Step 2 grid).

They fitted a statistical relationship (a polynomial regression) between these
two predictors and the *logarithm* of observed sediment thickness, using a
global compilation of sediment-thickness measurements. The workflow applies that
same fitted relationship to every ocean cell.

**The exact calculation.** For each cell, with seafloor age `A` (Ma) and
lifetime-mean distance-to-margin `D` (m — note the training used metres):

1. **Clamp to the calibrated range.** If `A > max_age`, set `A = max_age`; if
   `D > max_distance`, set `D = max_distance`. (`max_distance` is 3000 km, the
   same range discussed in Step 2.) This prevents the fit from being used as an
   extrapolation.

2. **Standardise** each predictor by subtracting its training mean and dividing
   by its training standard deviation (the square root of its training
   variance), exactly as the original machine-learning scaler did:

   ```
   a = (A − mean_age)      / √variance_age
   d = (D − mean_distance) / √variance_distance
   ```

3. **Evaluate the degree-3 polynomial** in the two standardised variables. The
   ten fitted coefficients `c₀…c₉` (in `config.yml`, `polynomial_coefficients`)
   multiply these ten terms, in this order:

   ```
   log_thickness = c₀·1     + c₁·a     + c₂·d
                 + c₃·a²     + c₄·a·d   + c₅·d²
                 + c₆·a³     + c₇·a²·d  + c₈·a·d²  + c₉·d³
   ```

4. **Return to linear space.** The fit is in natural-log space, so the predicted
   sediment thickness (metres) is:

   ```
   sediment_thickness = exp(log_thickness)
   ```

All of the constants — `mean_age`, `mean_distance`, `variance_age`,
`variance_distance`, `max_age`, `max_distance`, and the ten
`polynomial_coefficients` — live in the `sediment_thickness` block of
`config.yml`. The shipped values are the published Dutkiewicz et al. (2017)
relationship. If you retrain the relationship (for a different present-day
age grid or sediment-thickness dataset), replace these constants and the
workflow will use your new fit unchanged.

### Step 4 — Paleobathymetry with sediment isostasy

**Input:** the basement-depth grids (Step 1) and the sediment-thickness grids
(Step 3).
**Output:** `output/Paleobathymetry/paleobathymetry_<t>Ma.nc` — the final
sediment-covered ocean-floor depth (metres, negative down).

Adding a sediment pile of thickness `h` on top of the basement does *not* raise
the seafloor by the full `h`, because the added weight pushes the crust down
(isostasy). Sykes (1996) gives a correction for this. With sediment thickness `h`
expressed in **kilometres**:

```
correction(m) = (0.43422·h − 0.010395·h²) · 1000        (clamped to ≥ 0)
```

The final paleobathymetry is the basement depth plus the sediment thickness,
minus that isostatic correction:

```
paleobathymetry = basement_depth + sediment_thickness − correction
```

Because basement depth is negative (below sea level) and sediment thickness is a
positive amount of infill, the net effect is a shallower (less negative)
seafloor than the bare basement. This streamlined workflow does not add
large igneous provinces (LIPs) or seamounts, so there is no volcanic-height
term.

### Step 5 — pyBacktrack paleobathymetry (merged)

**Input:** the Step 4 paleobathymetry grids (read-only), plus the plate model
(rotations, static polygons, anchor plate) and the present-day seafloor-age
grid.
**Output:** `output/PaleobathymetryPyBacktrack/paleobathymetry_<t>Ma.nc` — the
same-format paleobathymetry grids, now also covering submerged continental
crust and crust that has since been subducted. This is **optional** (default
on; toggle with `run.steps.pybacktrack`) and requires the `pybacktrack`
package. The Step 4 grids in `output/Paleobathymetry/` are left untouched.

Steps 1–4 only reconstruct **oceanic** crust that is defined by the seafloor-age
grids; they leave continental crust and long-subducted ocean crust as `NaN`.
[pyBacktrack](https://github.com/EarthByte/pyBacktrack) (Müller, Cannon,
Williams & Dutkiewicz, 2018) fills that gap: it backtracks/backstrips a uniform
grid of synthetic drill sites on **present-day** crust — including submerged
continental crust — back through time. It cannot, however, generate
paleobathymetry for crust that no longer exists today (already subducted), so
this step **merges in** the Step 4 grids to fill those regions, using
pyBacktrack's own `reconstruct_paleo_bathymetry_grids()` merge support
(pyBacktrack's reconstructed values take precedence on crust that still exists
today; the Step 4 values fill in the rest).

To keep the two paleobathymetry sources aligned, Step 5 reuses the same
settings as Steps 1–4:

- the same plate model, rotations and anchor plate (`plate_model`);
- the shipped present-day seafloor-age grid (`model["age_grid"](0)`) — the
  static polygons needed to assign plate IDs are fetched from the same plate
  model via the Plate Model Manager (or from
  `plate_model.local.static_polygon_files` for a local plate model);
- the same ocean age → depth model (`age_depth.model`) — mapped onto
  pyBacktrack's equivalent constant (`gdh1` → `AGE_TO_DEPTH_MODEL_GDH1`,
  `rhcw18` → `AGE_TO_DEPTH_MODEL_RHCW18`, `crosby09` →
  `AGE_TO_DEPTH_MODEL_CROSBY_2007`; `parsons_sclater` has no pyBacktrack
  equivalent and is not supported by this step);
- the same output grid spacing and time range (`grids.output_spacing`,
  `time.min`/`max`/`step`).

pyBacktrack's own bundled global sediment-thickness, crustal-thickness and
topography grids are used as-is (as recommended when swapping plate models in
the pyBacktrack documentation) — only the plate model, ocean model, present-day
age grid, and the Step 4 merge grids are overridden.

---

## Configuration reference

Everything is driven by `config.yml`. The shipped file is fully commented; this
is a summary of the main blocks.

| Block | Key | Meaning |
|-------|-----|---------|
| `plate_model` | `use_pmm` | `true`: fetch the model from the GPlately Plate Model Manager. `false`: use your own local files (fill in the `local:` block). |
| | `name` | PMM model name (default `zahirovic2022`). |
| | `anchor_plate_id` | Plate held fixed in the reconstruction reference frame (default `0` = the model's absolute frame). Controls how ocean floor is reconstructed through time in Step 2, and the reference frame used in Step 5. |
| | `local.static_polygon_files` | Static polygon file(s) for a local plate model (`use_pmm: false`); required by Step 5 (pyBacktrack) to assign plate IDs. |
| `time` | `min`, `max`, `step` | Reconstruction times (Ma) to compute. |
| | `max_reconstruction_time` | Do not reconstruct ocean floor older than this; `null` = use the plate model's limit. |
| `proximity` | `use_continent_contouring` | `false` (default, off): distance to the COB line segments in `cob_line_segments`. `true`: use dynamically contoured passive margins. |
| | `cob_line_segments` | Path to the passive-margin **COB line-segment** file (default: the dataset shipped in `data/`). Must be polylines along passive margins — **not** COB polygons. |
| | `use_continent_obstacles` | `true` (default): shortest path *around* continents. `false`: straight great-circle distance. |
| | `clamp_mean_distance_km` | Cap on distance to margin (default 3000 km; see Step 2). |
| `grids` | `internal_spacing` | Grid spacing used for the internal distance computation (coarser = faster). |
| | `output_spacing` | Spacing of all final output grids (e.g. `0.1`°). |
| `age_grids` | `source` | `pmm` (age grids that ship with the PMM model) or `local` (your own, via `local_template`). |
| `age_depth` | `model` | `gdh1` (default), `rhcw18`, `parsons_sclater`, or `crosby09` (see Step 1). |
| | `richards_table` | Path to the RHCW18 lookup table (shipped in `data/`; only used by `rhcw18`). |
| `sediment_thickness` | (constants) | The Dutkiewicz et al. (2017) standardisation constants and polynomial coefficients (see Step 3). |
| `features` | `use_lips`, `use_seamounts` | Off by default; reserved for extensions. |
| `run` | `predicting_sediment_thickness_dir` | Where the public sediment-thickness engine lives (default `submodules/…`). |
| | `output_dir`, `num_cpus` | Output location and parallelism. |
| | `mask_continents` | `true` (default): blank continents in every output using the age-grid mask (see below). |
| | `steps` | Turn each of the five steps on/off individually. |

### Plate model — Plate Model Manager or your own files

By default the workflow fetches the **Zahirovic2022** model from the GPlately
Plate Model Manager (PMM), which supplies rotations, topologies, coastlines and
age grids automatically. (The passive-margin COB **line segments** do not come
from the PMM — they are the dataset shipped in `data/`; see Step 2.)

```yaml
plate_model:
  use_pmm: true
  name: zahirovic2022
```

To list the models the PMM offers:

```bash
python -c "from plate_model_manager import PlateModelManager as M; print(M().get_available_model_names())"
```

To use **your own plate model** instead, set `use_pmm: false`, fill in the
`local:` block (rotation, topology and coastline files), set
`age_grids.source: local`, and give an age-grid `local_template` (using
`{time}` for the reconstruction age in Ma). The COB line segments are set
separately under `proximity.cob_line_segments` and are independent of the plate
model.

---

## Inputs and outputs

**Inputs** (all supplied automatically when using the PMM defaults):

- a **plate model**: rotation files and topological features;
- **seafloor-age grids**, one per reconstruction time;
- **passive-margin COB line segments** (shipped in `data/`; or contoured margins);
- **continent obstacles** (coastlines) for the shortest-path distance;
- the shipped **RHCW18 age-depth table** (only for the `rhcw18` model).

**Outputs** (global netCDF grids, one file per time step, metres negative down):

```
output/
├── BasementDepth/               basement_depth_<t>Ma.nc           (Step 1)
├── Distances/                   mean_distance_<spacing>d_<t>.nc   (Step 2)
├── SedimentThickness/           sediment_thickness_<t>Ma.nc       (Step 3)
├── Paleobathymetry/             paleobathymetry_<t>Ma.nc          (Step 4)  ← ocean crust only
├── ContinentContours/           passive_margin_features.gpmlz     (only if continent contouring is on)
└── PaleobathymetryPyBacktrack/  paleobathymetry_<t>Ma.nc          (Step 5, optional) ← ocean + continental crust
```

### Continental masking (on by default)

Every product describes **ocean crust only**. The seafloor-age grids define
where ocean crust exists — cells with no age (`NaN`) are continental crust / no
ocean floor. The workflow propagates that mask to **every** output grid
(basement depth, distance, sediment thickness, paleobathymetry), so continents
are always blank. Basement depth (Step 1) gets `NaN` directly from the age→depth
conversion; the distance and sediment-thickness grids are additionally masked
with the age-grid mask (nearest-neighbour, so the coastline stays crisp); and
paleobathymetry inherits `NaN` from the masked basement + sediment grids. Turn
this off with `run.mask_continents: false` only if you deliberately want values
over continental crust. This masking applies to Steps 1–4 only; the optional
Step 5 (pyBacktrack) output is deliberately **not** masked, since it covers
submerged continental crust as well as ocean crust.

---

## Installation

### 1. Create the environment

```bash
conda env create -f environment.yml
conda activate simple_paleobathymetry
```

This installs `pygplates`, `gplately`, the `plate-model-manager`, `pybacktrack`
(for the optional Step 5), `GMT`, and the usual scientific-Python stack
(`numpy`, `scipy`, `pandas`, `xarray`, `netcdf4`, `joblib`, `pyyaml`).

### 2. Add the sediment-thickness engine

Steps 2 and 3 call two modules — `ocean_basin_proximity.py` (distances) and
`predict_sediment_thickness.py` (the polynomial prediction) — from the public
EarthByte **predicting-sediment-thickness** package. Add it as a git submodule:

```bash
git submodule add https://github.com/EarthByte/predicting-sediment-thickness \
    submodules/predicting-sediment-thickness
git submodule update --init --recursive
```

If you prefer, clone it anywhere and point `run.predicting_sediment_thickness_dir`
in the config at that path.

---

## Frequently asked questions

**Which age-depth model should I use?** `gdh1` is the standard default. Run the
others (`rhcw18`, `parsons_sclater`, `crosby09`) to gauge how sensitive your
results are to that choice, especially over old ocean crust.

**Why is my distance grid capped at 3000 km?** By design — the sediment-thickness
relationship is only calibrated within ~0–3000 km of a margin, so distances are
clamped there (`proximity.clamp_mean_distance_km`). See Step 2.

**Can I run only some steps?** Yes. Set the unwanted steps to `false` under
`run.steps`. For example, once distances exist you can re-run just Step 3 and
Step 4 after changing the sediment-thickness constants.

**Do I need to supply COB features?** No — a default passive-margin COB
line-segment dataset ships in `data/` and is used automatically. Override it via
`proximity.cob_line_segments` only if your plate model has severely altered COB
outlines (the Clennett / Alfonso family), or if you have a more appropriate COB
line-segment set. Do **not** point it at COB polygons (see Step 2).

**Do I need the internet?** Only when using the PMM (`use_pmm: true`), which
downloads and caches the plate model and age grids on first use. With a local
plate model and local age grids the workflow runs offline.

**What about submerged continental crust?** Steps 1–4 do not cover it (see
*Scope* above). The optional Step 5 does, using
[pyBacktrack](https://github.com/EarthByte/pyBacktrack) to reconstruct
paleobathymetry on submerged continental crust and on ocean crust preserved
today, merged with the Step 4 grids to also cover subducted ocean crust.

---

## See also

- **[pyBacktrack](https://github.com/EarthByte/pyBacktrack)** — backtracking /
  backstripping of drill-site and grid stratigraphy to reconstruct paleo-water
  depths for **submerged continental crust** as well as **ocean crust preserved
  today**. Used directly by the optional Step 5 of this workflow, which merges
  pyBacktrack's paleobathymetry with the Step 4 grids.
- **[predicting-sediment-thickness](https://github.com/EarthByte/predicting-sediment-thickness)**
  — the engine used here for distance-to-margin and sediment-thickness prediction
  (Steps 2–3).
- **[continent-contouring](https://github.com/EarthByte/continent-contouring)** —
  dynamically contoured passive margins for deep-time reconstructions (see Step 2).

---

## References

- **Stein, C.A. & Stein, S. (1992).** A model for the global variation in oceanic
  depth and heat flow with lithospheric age. *Nature*, 359, 123–129.
  *(GDH1 age–depth model.)*
- **Parsons, B. & Sclater, J.G. (1977).** An analysis of the variation of ocean
  floor bathymetry and heat flow with age. *Journal of Geophysical Research*,
  82, 803–827. *(Parsons & Sclater age–depth model.)*
- **Crosby, A.G., McKenzie, D. & Sclater, J.G. (2006);** and
  **Crosby, A.G. & McKenzie, D. (2009).** The relationship between depth, age
  and gravity in the oceans; and an analysis of young ocean depth, gravity and
  global residual topography. *Geophysical Journal International.*
  *(Crosby et al. age–depth model.)*
- **Richards, F.D., Hoggard, M.J., Cowton, L.R. & White, N.J. (2018).**
  Reassessing the thermal structure of oceanic lithosphere with revised
  global inventories of basement depths and heat flow measurements.
  *Journal of Geophysical Research: Solid Earth*, 123, 9136–9161.
  *(RHCW18 age–depth model; preferred parameters from the
  [RHCW18_Plate_Model](https://github.com/freddrichards/RHCW18_Plate_Model)
  repository.)*
- **Dutkiewicz, A., Müller, R.D., Wang, X., O'Callaghan, S., Cannon, J. &
  Wright, N.M. (2017).** Predicting sediment thickness on vanished ocean crust
  since 200 Ma. *Geochemistry, Geophysics, Geosystems*, 18, 4586–4603.
  *(Sediment-thickness relationship.)*
- **Sykes, T.J.S. (1996).** A correction for sediment load upon the ocean floor:
  uniform versus varying sediment density estimations. *Marine Geology*, 133,
  35–49. *(Isostatic sediment-load correction.)*
- **Müller, R.D., Cannon, J., Williams, S. & Dutkiewicz, A. (2018).**
  PyBacktrack 1.0: A tool for reconstructing paleobathymetry on oceanic and
  continental crust. *Geochemistry, Geophysics, Geosystems*, 19, 1898–1909,
  doi: [10.1029/2017GC007313](https://doi.org/10.1029/2017GC007313).
  *(pyBacktrack; Step 5.)*
- **Müller, R.D., Seton, M., Zahirovic, S., Williams, S.E., Matthews, K.J.,
  Wright, N.M., Shephard, G.E., Maloney, K.T., Barnett-Moore, N., Hosseinpour,
  M., Bower, D.J. & Cannon, J. (2016).** Ocean basin evolution and global-scale
  plate reorganization events since Pangea breakup. *Annual Review of Earth and
  Planetary Sciences*, 44, 107–138. *(Shipped present-day COB line-segment
  dataset, `data/Global_EarthByte_GPlates_PresentDay_COBs.gpmlz`.)*

---

## License

GNU General Public License, version 2 (see `LICENSE`).
