# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""
Export of simulation results to standard formats.

* VTK legacy ASCII (readable by ParaView / Mayavi / PyVista) with one scalar
  per triangle: hits, effective hits, exposure, max intensity, max impact
  pressure and a combined status field (0 = never reached, 1 = reached but
  not cleaned, 2 = effectively cleaned).
* JSON: the metric summary plus the parameters that produced it.
"""
from __future__ import annotations

import json
import os

import numpy as np

from .mesh import TriangleMesh
from .simulate import SimResults


def write_vtk(path: str, mesh: TriangleMesh, res: SimResults):
    """Write per-face coverage data as a VTK unstructured grid."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    n = mesh.n_tri
    status = np.where(res.n_effective > 0, 2, np.where(res.n_hits > 0, 1, 0))
    reach = (res.n_hits > 0).astype(float)
    clean = (res.n_effective > 0).astype(float)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# vtk DataFile Version 3.0\n")
        fh.write("hpr-sim cleaning coverage\n")
        fh.write("ASCII\n")
        fh.write("DATASET UNSTRUCTURED_GRID\n")
        fh.write(f"POINTS {n} double\n")
        for pt in mesh.centroids:
            fh.write(f"{pt[0]:.6g} {pt[1]:.6g} {pt[2]:.6g}\n")
        fh.write(f"CELLS {n} {n * 4}\n")
        for fa in mesh.faces:
            fh.write(f"3 {fa[0]} {fa[1]} {fa[2]}\n")
        fh.write(f"CELL_TYPES {n}\n")
        fh.write(("5\n" * n))
        fh.write(f"CELL_DATA {n}\n")
        names = [("status", status), ("reachable", reach), ("cleaned", clean),
                 ("n_hits", res.n_hits), ("n_effective", res.n_effective),
                 ("exposure_s", res.exposure_s),
                 ("max_intensity", res.max_intensity),
                 ("max_impact_bar", res.max_impact_bar)]
        for name, vals in names:
            fh.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
            for v in vals:
                fh.write(f"{v:.6g}\n")
    return path


def write_json(path: str, summary: dict, params: dict = None):
    """Write the metric summary (and optional parameters) as JSON."""
    obj = {"summary": summary}
    if params:
        obj["parameters"] = _jsonable(params)
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, default=_jsonable)
    return path


def read_vtk_scalars(path: str) -> dict:
    """Read the per-face SCALARS blocks written by :func:`write_vtk`.

    Returns ``{name: numpy array}``.  Only the legacy ASCII VTK files written
    by this package are supported.
    """
    n_faces = 0
    scalars = {}
    current = None
    counter = 0
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            key = parts[0]
            if key == "CELLS":
                n_faces = int(parts[1])
            elif key == "SCALARS":
                current = parts[1]
                scalars[current] = np.zeros(n_faces)
                counter = 0
            elif key == "LOOKUP_TABLE":
                continue
            elif current is not None and key not in ("POINTS", "CELL_TYPES",
                                                     "CELL_DATA", "DATASET"):
                for tok in parts:
                    if counter < n_faces:
                        scalars[current][counter] = float(tok)
                        counter += 1
                    else:
                        break
    return scalars


def _jsonable(x):
    if hasattr(x, "__dict__"):
        return {k: _jsonable(v) for k, v in x.__dict__.items() if not k.startswith("_")}
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    return str(x)
