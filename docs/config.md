# Configuration reference

All simulation settings live in **one YAML file** (JSON also works). Relative
paths in the file are resolved against the directory of the config file
itself.

Generate a fully commented template:

```bash
python -m hpr_sim --write-template configs/my_cavity.yaml
```

## Top-level sections

```yaml
cavity: tesla_9cell        # output file prefix
mesh:    { ... }           # where the cavity geometry comes from
water:   { ... }           # fluid properties
jet:     { ... }           # jet / droplet physics
nozzle:  { ... }           # spray head geometry
motion:  { ... }           # wand kinematics
sim:     { ... }           # cleaning criteria + resolution
output:  { ... }           # what to write and where
```

## `mesh`

| key | default | meaning |
|---|---|---|
| `stl` | — | path to a triangle-mesh STL (mm) |
| `step` | — | a CAD STEP file; auto-converted to STL with FreeCAD (optional) |
| `deflection_mm` | 0.5 | mesh deflection used during STEP→STL conversion |

## `water`

| key | default | meaning |
|---|---|---|
| `rho` | 1000 | water density [kg/m³] |
| `surface_tension` | 0.073 | [N/m] |
| `rho_air` | 1.2 | [kg/m³] |
| `mu_air` | 1.81e-5 | dynamic viscosity of air [Pa·s] |

## `jet`

| key | default | meaning |
|---|---|---|
| `pressure_bar` | 100 | pump pressure [bar] |
| `nozzle_diameter_mm` | 0.5 | hole diameter [mm] |
| `cd` | 0.7 | discharge coefficient |
| `velocity_model` | `droplet` | `droplet` (drag) or `linear` (exponential) |
| `drop_cd_mode` | `schiller-naumann` | drag correlation (see [physics.md](physics.md)) |

## `nozzle`

| key | default | meaning |
|---|---|---|
| `n_holes` | 6 | holes around the head |
| `tilt_fwd_deg` / `tilt_bwd_deg` | 25 / −25 | forward / backward tilt [deg] |
| `n_spray` | 11 | rays per hole (cone sampling) |
| `jet_angle_deg` | 5 | spray half-angle [deg] |
| `head_diameter_mm` | 25 | mechanical head diameter [mm] |

## `motion`

| key | default | meaning |
|---|---|---|
| `x_start_mm` | `auto` | scan start; `auto` detects the cell chain |
| `x_end_mm` | `auto` | scan end |
| `trans_speed_mm_s` | 10 | wand translation speed [mm/s] |
| `rot_speed_rpm` | 24 | head rotation [rpm] |
| `n_recip` | 2 | reciprocation count (passes) |

## `sim`

| key | default | meaning |
|---|---|---|
| `step_trans_mm` | 1.0 | axial pose spacing [mm] |
| `step_rot_deg` | 5 | azimuth pose spacing [deg] |
| `dt_s` | 0.05 | time step between poses [s] |
| `clean_cos_thresh` | 0.5 | minimum `|cos φ|` for effective cleaning |
| `d_min_mm` | 5 | minimum usable stand-off [mm] |
| `d_max_mm` | 200 | maximum usable stand-off [mm] |
| `dmax_kd` | 0 | if >0, `d_max = k_d D` (Guha stand-off) overrides `d_max_mm` |
| `min_impact_bar` | 0 | require wall impact pressure ≥ this [bar] |
| `ray_len_max_mm` | `auto` | ray length; `auto` = `r_eq·√(1+tan²(tilt)) + margin` |

## `output`

| key | default | meaning |
|---|---|---|
| `dir` | `../results` | output directory |
| `vtk` | true | write `*_coverage.vtk` (ParaView) |
| `plots` | true | write PNG plots |
| `json` | true | write `*_summary.json` |

## `auto` geometry detection

When `motion.x_start_mm`, `motion.x_end_mm` or `sim.ray_len_max_mm` are the
string `auto`, the values are inferred from the mesh:

1. radial profile `r_min(x)` of the wall;
2. the cell chain = all points where `r_min` crosses the iris level
   (0.7 × largest equator radius) — multi-cell cavities are handled;
3. scan range = outermost irises ± margin; ray length from the largest cell
   and the nozzle tilt.

`geometry.propose_scan_range()` exposes the same function programmatically.
