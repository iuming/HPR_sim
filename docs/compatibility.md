# Compatibility & Version Matrix

Supported environments and feature requirements of `hpr-sim` (v0.2.0).

This matrix reflects what is **unit-tested in CI** (Linux / Windows / macOS,
Python 3.9–3.12) and **runtime-tested locally** (Windows 11, Python 3.13.9,
numpy 2.3.5, matplotlib 3.10.6, PyYAML 6.0.3, FreeCAD 1.1.1 headless).

## 1. Python & OS

| Python | Linux | Windows | macOS | Status |
|---|---|---|---|---|
| 3.9  | ✅ CI | ✅ CI | ✅ CI | supported (lower bound) |
| 3.10 | ✅ CI | ✅ CI | ✅ CI | supported |
| 3.11 | ✅ CI | ✅ CI | ✅ CI | supported |
| 3.12 | ✅ CI | ✅ CI | ✅ CI | supported |
| 3.13 | 🟡 runtime-tested | 🟡 runtime-tested | — | works, not in CI yet |

> The package requires Python **≥ 3.9** (`requires-python = ">=3.9"`) and
> avoids 3.10+ syntax; nothing on the NumPy 2.x line is assumed either (1.22+
> is fine).

## 2. Core dependencies

| Package | Min | Recommended | Tested against | Required for |
|---|---|---|---|---|
| `numpy` | ≥ 1.22 | latest 2.x | 2.3.5 | core simulation |
| `matplotlib` | ≥ 3.5 | latest 3.x | 3.10.6 | plots only (headless-safe) |
| `PyYAML` | ≥ 5.4 | latest | 6.0.3 | YAML configuration (JSON works without it) |

Install with:

```bash
pip install -e .
```

## 3. Optional components

| Component | Version | Used for | Setup |
|---|---|---|---|
| [FreeCAD](https://www.freecad.org) (`freecadcmd`) | ≥ 0.20 (tested 1.1.1) | STEP → STL conversion | on PATH, or `FREECAD_BIN=<path>/freecadcmd.exe` |
| [ParaView](https://www.paraview.org) | ≥ 5.x (any recent) | viewing `*_coverage.vtk` | none (output-only) |
| [Jupyter](https://jupyter.org) | any | `examples/hpr_sim_demo.ipynb` | `pip install jupyter` |
| `pytest` | ≥ 7 | unit tests | `pip install -e .[dev]` |

## 4. Feature-by-feature matrix

| Feature | Requires | Notes |
|---|---|---|
| Core ray-casting simulation | numpy | pure Python + numpy; no compiled deps |
| YAML configuration | PyYAML | JSON configs work without PyYAML |
| Plots & sweeps | numpy, matplotlib | `matplotlib.use('Agg')` headless-safe |
| STEP → STL | FreeCAD (`freecadcmd`) | skipped/errs with a clear message if absent |
| Parallel sweeps (`--jobs N`) | OS process/pipe support | auto-falls back to serial if the pool cannot be created (e.g. restricted sandboxes) |
| VTK export | — | hand-written legacy VTK ASCII, no VTK package needed |
| Demo notebook | jupyter | — |
| Mesh auto-detection & checks | numpy | — |

## 5. Determinism & numerics

* The simulation is **deterministic** (`sim.seed = 0` → identical results on a
  given machine/blas).
* Cross-platform agreement is at the level of the last displayed digit:
  wall-clock times differ, coverage ratios agree to ≳0.1 % (BLAS rounding).
* The internal units are mm / s / bar / m/s; STL is assumed in **mm**
  (`meshcheck` flags wildly out-of-scale inputs).

## 6. Known platform notes

| Platform | Note |
|---|---|
| Windows | parallel sweeps require OS named-pipe/IPC access; restricted shells (e.g. sandboxed CI executors) may deny it → graceful serial fallback |
| macOS | FreeCAD path resolution: set `FREECAD_BIN` if `freecadcmd` is not on PATH |
| Any headless Linux | set `MPLBACKEND=Agg` or rely on the built-in `matplotlib.use('Agg')` fallback (already applied in CLI/notebook) |
