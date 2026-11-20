# -*- coding: utf-8 -*-
"""
YAML/JSON configuration handling.

A configuration file describes one complete simulation case:

.. code-block:: yaml

    cavity: my_cavity          # name used for outputs
    mesh:
      step: ../path/to/model.STEP   # optional STEP -> auto converted
      stl:  ../path/to/generated.stl
      deflection_mm: 0.5            # STEP meshing tolerance
    water:  {rho: 1000, ...}        # see hpr_sim.config.WaterParams
    jet:    {pressure_bar: 100, ...}
    nozzle: {n_holes: 6, ...}
    motion: {x_start_mm: -160, ...} # motion.path auto-detected if omitted
    sim:    {clean_cos_thresh: 0.5, ...}
    output: {dir: ../results}

Relative paths are resolved against the directory of the config file.
"""
from __future__ import annotations

import json
import os

try:
    import yaml
    _HAS_YAML = True
except Exception:  # pragma: no cover
    _HAS_YAML = False

from .config import to_dict


def _resolve(path: str, base_dir: str) -> str:
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(base_dir, path))


def load_config(path: str) -> dict:
    """Load a YAML (or JSON) configuration file.

    Returns the raw config dict with all relative paths resolved to absolute
    paths (keys ``mesh.step``, ``mesh.stl``, ``output.dir``).
    """
    base = os.path.dirname(os.path.abspath(path))
    with open(path, "r", encoding="utf-8") as fh:
        if path.lower().endswith(".json"):
            cfg = json.load(fh)
        else:
            if not _HAS_YAML:
                raise RuntimeError("PyYAML is required for .yaml config files")
            cfg = yaml.safe_load(fh)
    cell = cfg.setdefault("mesh", {})
    for key in ("step", "stl"):
        if cell.get(key):
            cell[key] = _resolve(str(cell[key]), base)
    out = cfg.setdefault("output", {})
    if out.get("dir"):
        out["dir"] = _resolve(str(out["dir"]), base)
    return cfg


def save_config(path: str, cfg: dict) -> str:
    """Write a config dict back to a YAML file."""
    if not _HAS_YAML:
        raise RuntimeError("PyYAML is required for .yaml config files")
    with open(path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, allow_unicode=True)
    return path


def template_config() -> dict:
    """A fully-commented template configuration."""
    from .config import CONFIG_SECTIONS
    cfg = {
        "# info": "hpr-sim configuration template, see hpr_sim/config_io.py",
        "cavity": "my_cavity",
        "mesh": {
            "step": None,
            "stl": "../converted_meshes/my_cavity.stl",
            "deflection_mm": 0.5,
        },
        "output": {"dir": "../results", "vtk": True},
    }
    for name, cls in CONFIG_SECTIONS.items():
        cfg[name] = to_dict(cls())
    return cfg


def write_reference_config(path: str) -> str:
    """Write the template configuration to `path`."""
    return save_config(path, template_config())
