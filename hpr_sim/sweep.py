# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""
Parameter sweeps of the HPR simulation, optionally parallelised.

A sweep runs a series of simulations with one parameter varied.  Because a
cavity mesh can be large, each worker process loads the mesh once (from a
file path) through the pool initializer, and the parameter configurations
(plain dataclasses) are shipped over the pipe -- small and fast.
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .analyze import summarize
from .mesh import TriangleMesh
from .simulate import run_simulation

SWEEPS = {
    "threshold": dict(attr=("sim", "clean_cos_thresh"),
                      vals=[0.3, 0.4, 0.5, 0.6, 0.7],
                      label="cleaning intensity threshold (cos)"),
    "tilt": dict(attr=("nozzle", "tilt_fwd_deg"),
                 vals=[15, 20, 25, 30, 35], label="hole tilt angle [deg]"),
    "recip": dict(attr=("motion", "n_recip"),
                  vals=[1, 2, 3, 4, 6], label="reciprocation count"),
    "speed": dict(attr=("motion", "trans_speed_mm_s"),
                  vals=[5, 10, 20, 40, 80],
                  label="wand translation speed [mm/s]"),
    "rot": dict(attr=("motion", "rot_speed_rpm"),
                vals=[6, 12, 24, 48, 96], label="head rotation speed [rpm]"),
    "jetangle": dict(attr=("nozzle", "jet_angle_deg"),
                     vals=[3, 5, 8, 12, 16], label="jet spray half-angle [deg]"),
    "pressure": dict(attr=("jet", "pressure_bar"),
                     vals=[70, 85, 100, 130, 160], label="pump pressure [bar]"),
    "impact": dict(attr=("sim", "min_impact_bar"),
                   vals=[0, 1, 2, 3, 5], label="min wall impact pressure [bar]"),
}


def _clone(section):
    """Deep-copy a dataclass instance by reconstructing it."""
    import copy
    return copy.deepcopy(section)


def _set(p, path, value):
    obj = p
    for k in path[:-1]:
        obj = obj[k]
    if path[-1] == "n_recip":
        value = int(round(value))
    setattr(obj, path[-1], value)


_WORKER = {}


def _init(mesh_path, base_cfg):
    _WORKER["mesh"] = TriangleMesh.from_stl(mesh_path)
    _WORKER["cfg"] = base_cfg
    _WORKER["roi"] = None


def _run_one(job):
    attr, fast, value = job
    p = {k: _clone(v) for k, v in _WORKER["cfg"].items()}
    if fast:
        p["sim"].dt_s = 0.2
        p["sim"].step_rot_deg = 20.0
        p["motion"].n_recip = min(p["motion"].n_recip, 1)
    _set(p, attr, value)
    motion = p["motion"]
    cfg = p["sim"]
    roi = (motion.x_start_mm - 30.0, motion.x_end_mm + 30.0)
    res = run_simulation(_WORKER["mesh"], p["water"], p["jet"], p["nozzle"],
                         motion, cfg, axis="x", progress=False)
    s = summarize(_WORKER["mesh"], res, roi=roi)
    s["parameter"] = float(value)
    return s


def sweep_variant(mesh: TriangleMesh, p, attr, value, fast=False):
    """Run a single sweep point in-process (useful for testing / jobs=1)."""
    p2 = {k: _clone(v) for k, v in p.items()}
    if fast:
        p2["sim"].dt_s = 0.2
        p2["sim"].step_rot_deg = 20.0
        p2["motion"].n_recip = min(p2["motion"].n_recip, 1)
    _set(p2, attr, value)
    roi = (p2["motion"].x_start_mm - 30.0, p2["motion"].x_end_mm + 30.0)
    res = run_simulation(mesh, p2["water"], p2["jet"], p2["nozzle"],
                         p2["motion"], p2["sim"], axis="x", progress=False)
    s = summarize(mesh, res, roi=roi)
    s["parameter"] = float(value)
    return s


def run_sweep(mesh_path: str, p, sweep_name: str, vals=None, jobs: int = 1,
              fast: bool = False, progress: bool = True):
    """Run a sweep of the HPR simulation.

    Parameters
    ----------
    mesh_path : STL file (loaded once per worker process)
    p : dict with water/jet/nozzle/motion/sim dataclasses
    sweep_name : key of :data:`SWEEPS`
    vals : values to test (default from SWEEPS)
    jobs : number of parallel workers (ProcessPoolExecutor); if the pool
           cannot be created (e.g. restricted/hostile sandboxes that forbid
           OS pipes), the sweep transparently falls back to in-process serial
           execution.
    fast : reduce resolution for quick exploratory sweeps
    progress : print one line per completed point
    """
    spec = SWEEPS[sweep_name]
    attr = spec["attr"]
    if vals is None:
        vals = list(spec["vals"])
    if jobs <= 1:
        mesh = TriangleMesh.from_stl(mesh_path)
        rows = [sweep_variant(mesh, p, attr, v, fast=fast) for v in vals]
        if progress:
            for v, s in zip(vals, rows):
                print(f"  [{sweep_name}={v}] coverage="
                      f"{s['effective_coverage_ratio'] * 100:.2f}% "
                      f"({s['wall_time_s']:.0f}s)", flush=True)
        return rows, spec

    os.makedirs(os.path.dirname(os.path.abspath(mesh_path)) or ".", exist_ok=True)
    base_cfg = {k: _clone(v) for k, v in p.items()}
    jobs = min(jobs, len(vals))
    try:
        executor = ProcessPoolExecutor(max_workers=jobs, initializer=_init,
                                       initargs=(os.path.abspath(mesh_path),
                                                 base_cfg))
    except Exception as exc:  # pragma: no cover - sandbox dependent
        print(f"  [sweep] parallel pool unavailable ({exc!r}); "
              f"falling back to serial execution", flush=True)
        return run_sweep(mesh_path, p, sweep_name, vals=vals, jobs=1,
                         fast=fast, progress=progress)
    with executor as ex:
        futures = [ex.submit(_run_one, (attr, fast, v)) for v in vals]
        rows = [None] * len(futures)
        for i, f in enumerate(futures):
            s = f.result()
            rows[i] = s
            if progress:
                print(f"  [{sweep_name}={s['parameter']}] coverage="
                      f"{s['effective_coverage_ratio'] * 100:.2f}% "
                      f"({s['wall_time_s']:.0f}s)", flush=True)
    return rows, spec
