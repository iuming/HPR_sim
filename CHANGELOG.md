# Changelog

All notable changes to **hpr-sim** are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/) and this
project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Standalone open-source repository layout: `docs/`, `examples/`, `tools/`,
  `data/`, CI workflows, contributing guide, code of conduct, changelog.
- **TESLA 1.3 GHz 9-cell example**: parametric cavity builder
  (`tools/make_tesla_cavity.py`) based on the shape published by
  Aune et al. (PRST-AB 3, 092001) — water-tight inner-surface mesh, example
  YAML, output assets.
- Multi-cell aware auto-detection: the scan region now spans the whole cell
  chain (first → last iris) instead of only the largest single cell.
- 2-D strict-1:1 half-cavity coverage top view (`coverage_topview`) and a
  `tools/replot.py` helper to re-render plots from saved VTK data.
- `hpr_sim/export.read_vtk_scalars()` to read the coverage fields back.

### Changed
- Coverage cloud maps: layered vivid colours (green = cleaned,
  amber = reached-not-clean, faint outline), ROI filtering, consistent
  (extent-proportional) 3-D axes.
- Tests write into `data/_test/` (git-ignored) instead of `results/`.

## [0.2.0] — 2026-11-19

### Added
- Public package API: `hpr_sim.{config,config_io,geometry,mesh,meshcheck,
  motion,nozzle,physics,raycast,simulate,analyze,visualize,export,sweep,cli}`.
- YAML/JSON configuration with relative-path resolution and a commented
  template generator (`--write-template`).
- Automatic beam-axis / scan-region / ray-length detection.
- Mesh self-check (water-tightness via deduplicated shared edges, normals,
  degenerate faces, units).
- Droplet-drag jet model with `schiller-naumann`, `morsi-alexander`, `putnam`
  and `constant` drag laws; RCK2 numerical integration.
- VTK export, JSON summaries, parameter sweeps with optional process-pool
  parallelism (graceful serial fallback on restricted hosts).
- `pyproject.toml`, console script `hpr-sim`, pytest suite (24 tests).

## [0.1.0] — 2026-11

- Initial standalone scripts (`run_simulation.py`, `run_sweep.py`) and the
  HPR ray-casting physics developed for the 500 MHz / 1.5 GHz cavity
  vertical-test campaign.
