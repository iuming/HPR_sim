# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Quickstart for the TESLA 9-cell example: build the cavity, run the
simulation with the provided YAML and print the headline result.

    python examples/run_tesla_quickstart.py
"""
from __future__ import annotations

import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from hpr_sim.config_io import load_config
from hpr_sim.cli import run_from_config

STL = os.path.join(ROOT, "data", "tesla_9cell_inner.stl")


def _args():
    ns = argparse.Namespace()
    ns.config = os.path.join(ROOT, "examples", "tesla_9cell.yaml")
    ns.template = None
    ns.sweep = None
    ns.jobs = 1
    ns.fast = False
    return ns


def main():
    # 1) make sure the cavity mesh exists
    if not os.path.exists(STL):
        from tools.make_tesla_cavity import main as _build
        sys.argv = ["make_tesla_cavity.py", "--stl", STL]
        _build()

    # 2) run the whole simulation from the YAML config
    cfg_file = os.path.join(ROOT, "examples", "tesla_9cell.yaml")
    cfg = load_config(cfg_file)
    summary = run_from_config(cfg, _args())
    print("-" * 62)
    print("TESLA 9-cell effective coverage : "
          f"{summary['effective_coverage_ratio'] * 100:.2f}%")
    print("reachable inner area            : "
          f"{summary['reachable_area_m2']:.4f} m^2")
    print("outputs written to              : data/outputs/")


if __name__ == "__main__":
    main()
