# hpr-sim documentation

Welcome to hpr-sim, a general-purpose **high-pressure rinsing (HPR) coverage
simulation** for superconducting RF (SRF) cavities.

```bash
pip install -e .
python -m hpr_sim --help
```

## Contents

| Topic | File |
|---|---|
| Concepts and the physics model | [physics.md](physics.md) |
| YAML configuration reference | [config.md](config.md) |
| Supported environments | [compatibility.md](compatibility.md) |
| Version history | [../CHANGELOG.md](../CHANGELOG.md) |
| Contribution guide | [../CONTRIBUTING.md](../CONTRIBUTING.md) |

## Quick start (10 seconds)

1. Build the example 9-cell TESLA cavity mesh and run the simulation:

   ```bash
   python tools/make_tesla_cavity.py
   python -m hpr_sim -f examples/tesla_9cell.yaml
   ```

2. Open the outputs in `data/outputs/`:

   * `tesla_9cell_coverage3d.png` — effective cleaning coverage cloud;
   * `tesla_9cell_unrolled.png` — axial–azimuth coverage map;
   * `tesla_9cell_coverage.vtk` — per-face data for [ParaView](https://www.paraview.org);
   * `tesla_9cell_summary.json` — metrics.

3. Try a parameter sweep:

   ```bash
   python -m hpr_sim -f examples/tesla_9cell.yaml --sweep threshold,tilt --jobs 4
   ```

## What it computes

For every nozzle pose of the translating, rotating wand, `hpr-sim` casts the
water jets onto the cavity mesh, finds the first hit with a vectorised
Möller–Trumbore intersection, and classifies it as *effectively cleaned* when

* the incidence cosine `|cos φ|` is above the threshold, and
* the jet at that distance still delivers the required impact (droplet-drag
  model), and
* (optional) the wall impact pressure is above `min_impact_bar`.

The headline output is the **effective cleaning coverage** =
cleaned area / reachable area of the RF surface.
