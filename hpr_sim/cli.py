# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""
Command-line entry point of the general HPR simulation package.

Usage
-----
    python -m hpr_sim -f examples/tesla_9cell.yaml
    python -m hpr_sim -f examples/tesla_9cell.yaml --sweep threshold,tilt --jobs 4
    python -m hpr_sim --write-template configs/new_cavity.yaml

The same functionality is available as the ``hpr-sim`` console script after
``pip install -e .``.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

from .analyze import print_summary
from .config import parse_config_dict
from .config_io import load_config, template_config, save_config
from .export import write_json, write_vtk
from .geometry import apply_auto
from .mesh import TriangleMesh
from .meshcheck import check_mesh, print_report
from .simulate import run_simulation
from .sweep import run_sweep
from . import visualize


def _ensure_mesh(cfg: dict, args) -> str:
    """Return an existing STL path, converting from STEP if necessary."""
    mesh_section = cfg.get("mesh", {})
    stl = mesh_section.get("stl")
    step = mesh_section.get("step")
    if stl and os.path.exists(stl):
        return stl
    if step and os.path.exists(step):
        print(f"[mesh] converting STEP -> STL ({step})", flush=True)
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
        from step_to_stl import step_to_stl
        out = step_to_stl(step, stl,
                          lin_def=float(mesh_section.get("deflection_mm", 0.5)),
                          ang_def=10.0, quiet=True)
        return out
    raise SystemExit(f"mesh not found: stl={stl} step={step}")


def save_sweep_outputs(cavity: str, out_dir: str, sweep_name: str, rows, spec):
    import csv
    base = os.path.join(out_dir, f"{cavity}_sweep_{sweep_name}")
    with open(base + ".csv", "w", newline="", encoding="utf-8") as fh:
        keys = ["parameter", "coverage", "cleaned_area_m2", "reachable_area_m2",
                "wall_time_s", "mean_exposure_s"]
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows:
            w.writerow([r["parameter"], r["effective_coverage_ratio"],
                        r["cleaned_area_m2"], r["reachable_area_m2"],
                        r["wall_time_s"], r["mean_exposure_cleaned_s"]])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = [r["parameter"] for r in rows]
    y = [r["effective_coverage_ratio"] * 100 for r in rows]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x, y, "o-", color="crimson", lw=2)
    ax.set_xlabel(spec["label"]); ax.set_ylabel("effective coverage [%]")
    ax.set_ylim(0, 100)
    ax.set_title(f"{cavity}: HPR coverage vs {spec['label']}")
    fig.tight_layout(); fig.savefig(base + ".png", dpi=140); plt.close(fig)
    print(f"[sweep] saved {base}.csv, {base}.png", flush=True)


def run_from_config(cfg: dict, args) -> dict:
    cavity = cfg.get("cavity", "case")
    out_dir = cfg.get("output", {}).get("dir", "results")
    os.makedirs(out_dir, exist_ok=True)

    stl = _ensure_mesh(cfg, args)
    print(f"[{cavity}] loading mesh {stl}")
    mesh = TriangleMesh.from_stl(stl, name=cavity)
    print("  " + mesh.summary())

    # --- mesh self-checks ------------------------------------------------
    rep = check_mesh(mesh)
    print_report(rep)

    # --- parameters -------------------------------------------------------
    p = parse_config_dict(cfg)
    axis = p["motion"].axis
    prop = apply_auto(mesh, p["motion"], p["sim"], axis=axis)
    if prop:
        print(f"[auto] scan region x=[{p['motion'].x_start_mm:.1f}, "
              f"{p['motion'].x_end_mm:.1f}] mm, equator@x={prop['equator_mm']:.1f} "
              f"(r={prop['equator_radius_mm']:.1f} mm), ray_len={p['sim'].ray_len_max_mm:.1f} mm",
              flush=True)

    res = run_simulation(mesh, p["water"], p["jet"], p["nozzle"], p["motion"],
                         p["sim"], axis=axis, progress=True)
    roi = (p["motion"].x_start_mm - 30.0, p["motion"].x_end_mm + 30.0)
    summary = print_summary(mesh, res, roi=roi)
    summary["cavity"] = cavity
    summary["roi_mm"] = list(roi)
    summary["mesh_check"] = rep["score"]

    # --- outputs ----------------------------------------------------------
    plot_enabled = cfg.get("output", {}).get("plots", True)
    if plot_enabled:
        base = os.path.join(out_dir, cavity)
        visualize.coverage_cloud(base + "_coverage3d.png", mesh, res,
                                 title=f"{cavity}: effective HPR coverage",
                                 roi_mm=roi)
        visualize.intensity_cloud(base + "_intensity3d.png", mesh, res,
                                  title=f"{cavity}: max cleaning intensity")
        visualize.unwrapped_map(base + "_unrolled.png", mesh, res,
                                title=f"{cavity}: coverage (effective/reachable)")
        visualize.histograms(base + "_hist.png", mesh, res)
        print(f"[{cavity}] plots -> {base}_*.png", flush=True)
    if cfg.get("output", {}).get("vtk", False):
        write_vtk(os.path.join(out_dir, cavity + "_coverage.vtk"), mesh, res)
        print(f"[{cavity}] vtk   -> {os.path.join(out_dir, cavity + '_coverage.vtk')}",
              flush=True)
    if cfg.get("output", {}).get("json", True):
        write_json(os.path.join(out_dir, cavity + "_summary.json"), summary,
                   params={k: v for k, v in p.items()})
        print(f"[{cavity}] json  -> {os.path.join(out_dir, cavity + '_summary.json')}",
              flush=True)

    # --- parameter sweeps ------------------------------------------------
    sweeps = list(args.sweep) if args.sweep else []
    for sw in sweeps:
        if sw not in ("threshold", "tilt", "recip", "speed", "rot",
                      "jetangle", "pressure", "impact"):
            print(f"[sweep] unknown sweep '{sw}', skipping", flush=True)
            continue
        rows, spec = run_sweep(stl, p, sw, jobs=args.jobs, fast=args.fast)
        save_sweep_outputs(cavity, out_dir, sw, rows, spec)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="hpr-sim",
        description="General-purpose high-pressure rinsing (HPR) coverage "
                    "simulation for SRF cavities.")
    ap.add_argument("-f", "--config", help="YAML/JSON configuration file")
    ap.add_argument("--write-template", metavar="PATH",
                    help="write a commented template config and exit")
    ap.add_argument("--sweep", default="",
                    help="comma-separated sweeps to run after the baseline "
                         "(threshold,tilt,recip,speed,rot,jetangle,pressure,impact)")
    ap.add_argument("--jobs", type=int, default=1,
                    help="parallel workers for the sweeps (default 1)")
    ap.add_argument("--fast", action="store_true",
                    help="reduced resolution for quick sweeps")
    args = ap.parse_args(argv)

    if args.write_template:
        cfg = template_config()
        save_config(args.write_template, cfg)
        print(f"template written to {args.write_template}")
        return

    if not args.config:
        ap.error("either -f/--config or --write-template is required")
    cfg = load_config(args.config)
    args.sweep = [s for s in args.sweep.split(",") if s.strip()]
    run_from_config(cfg, args)


if __name__ == "__main__":
    main()
