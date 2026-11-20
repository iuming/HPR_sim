#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Regenerate the result plots from a previously saved VTK coverage file.

This avoids re-running the (minutes-long) simulation when only the plots are
changed (e.g. colour scheme, viewpoint).

Usage
-----
    python tools/replot.py <mesh.stl> <out.vtk> <png-prefix> [--title ...] [--roi x0,x1]
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from hpr_sim.mesh import TriangleMesh
from hpr_sim.simulate import SimResults
from hpr_sim.export import read_vtk_scalars
from hpr_sim import visualize


def results_from_vtk(mesh, path):
    d = read_vtk_scalars(path)
    res = SimResults(mesh.n_tri)
    name_map = {
        "n_hits": "n_hits", "n_effective": "n_effective",
        "exposure_s": "exposure_s", "max_intensity": "max_intensity",
        "max_impact_bar": "max_impact_bar",
    }
    for k, attr in name_map.items():
        if k in d:
            setattr(res, attr, d[k])
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("vtk")
    ap.add_argument("prefix")
    ap.add_argument("--title", default=None)
    ap.add_argument("--roi", default=None, help="x0,x1 in mm")
    ap.add_argument("--cut", default=None, choices=["x", "y", "z"],
                    help="half-cavity cutaway view along a plane through the centre")
    ap.add_argument("--view", default="3d", choices=["3d", "top", "side"],
                    help="'3d' perspective view, or flat 2-D 'top' (strict 1:1) view of a cutaway")
    ap.add_argument("--out", default=None, help="output prefix override")
    args = ap.parse_args()

    mesh = TriangleMesh.from_stl(args.mesh)
    res = results_from_vtk(mesh, args.vtk)
    roi = None
    if args.roi:
        x0, x1 = [float(x) for x in args.roi.split(",")]
        roi = (x0, x1)
    title = args.title or os.path.basename(args.prefix)
    prefix = args.out or args.prefix
    cut = args.cut or "z"

    if args.view == "top":
        from hpr_sim.visualize import coverage_topview
        coverage_topview(prefix + "_coverage3d.png", mesh, res, cut=cut,
                         roi_mm=roi, title=title + " (top view, 1:1)")
    else:
        visualize.coverage_cloud(prefix + "_coverage3d.png", mesh, res,
                                 title=title, roi_mm=roi, cut=cut)
    visualize.intensity_cloud(prefix + "_intensity3d.png", mesh, res,
                              title=title)
    visualize.unwrapped_map(prefix + "_unrolled.png", mesh, res,
                            title=title)
    visualize.histograms(prefix + "_hist.png", mesh, res)
    print("regenerated:", prefix + "_*.png")


if __name__ == "__main__":
    main()
