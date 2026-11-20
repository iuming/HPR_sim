# -*- coding: utf-8 -*-
"""
Mesh quality checks run before the simulation.

These checks catch the most common problems with imported meshes:
  * open (non-watertight) shells, which make ray casting leak,
  * inconsistent face winding (normals pointing inwards/outwards mixed),
  * degenerate / zero-area or non-finite triangles,
  * wildly different length scales (units sanity check).
"""
from __future__ import annotations

import numpy as np

from .mesh import TriangleMesh

_DEGEN_AREA = 1e-12  # mm^2


def check_mesh(mesh: TriangleMesh, expected_scale_mm: tuple = (1.0, 1e4)) -> dict:
    """Run all checks and return a report dict."""
    rep = {"n_faces": int(mesh.n_tri), "n_points": len(mesh.vertices)}

    # ---- finite vertices --------------------------------------------------
    finite_pts = np.isfinite(mesh.vertices).all(axis=1)
    rep["n_bad_points"] = int((~finite_pts).sum())

    # ---- degenerate / zero-area triangles ---------------------------------
    rep["n_degenerate"] = int(getattr(mesh, "degenerate", np.zeros(mesh.n_tri, bool)).sum())

    # ---- bounding box / units ---------------------------------------------
    bbmin, bbmax = mesh.vertices.min(axis=0), mesh.vertices.max(axis=0)
    dims = bbmax - bbmin
    rep["bbox_mm"] = [float(bbmin[0]), float(bbmin[1]), float(bbmin[2]),
                      float(bbmax[0]), float(bbmax[1]), float(bbmax[2])]
    rep["extent_mm"] = [float(d) for d in dims]
    longest = float(dims.max())
    rep["units_ok"] = expected_scale_mm[0] <= longest <= expected_scale_mm[1]

    # ---- water-tightness (shared edges on deduplicated coordinates) --------
    n = mesh.n_tri
    if n > 0:
        vq = np.round(mesh.vertices, 5)                 # 0.1 um tolerance
        _, inv = np.unique(vq, axis=0, return_inverse=True)
        fid = inv[mesh.faces]                           # (M,3) vertex ids
        e01 = np.sort(fid[:, [0, 1]], axis=1)
        e12 = np.sort(fid[:, [1, 2]], axis=1)
        e20 = np.sort(fid[:, [2, 0]], axis=1)
        edges = np.vstack([e01, e12, e20])
        K = len(np.unique(inv))
        keys = edges[:, 0] * (K + 1) + edges[:, 1]
        vals, counts = np.unique(keys, return_counts=True)
        n_open = int((counts == 1).sum())
        n_manifold = int((counts == 2).sum())
        n_multi = int((counts > 2).sum())
        rep["n_unique_vertices"] = K
        rep["n_edges"] = len(vals)
        rep["n_edges_open"] = n_open
        rep["n_edges_nonmanifold"] = n_multi
        rep["watertight"] = (n_open == 0) and (n_multi == 0)
        rep["watertight_with_nonmanifold"] = (n_open == 0)
    else:
        rep["n_unique_vertices"] = 0
        rep["n_edges"] = 0
        rep["n_edges_open"] = 0
        rep["n_edges_nonmanifold"] = 0
        rep["watertight"] = False
        rep["watertight_with_nonmanifold"] = False

    # ---- normal orientation (outward = consistent with outer surface) -----
    if n > 0:
        vol_center = 0.5 * (bbmin + bbmax)
        outward = np.einsum("ij,ij->i", mesh.centroids - vol_center,
                            mesh.normals_unit)
        rep["n_faces_inward"] = int((outward < 0).sum())
        rep["normals_mixed"] = not (rep["n_faces_inward"] in (0, n))
    else:
        rep["n_faces_inward"] = 0
        rep["normals_mixed"] = False

    rep["score"] = _score(rep)
    return rep


def _score(rep: dict) -> str:
    if rep["n_faces"] == 0:
        return "FAIL: empty mesh"
    if not rep["units_ok"]:
        return "FAIL: bounding box outside expected mm scale"
    if rep["n_bad_points"] or rep["n_degenerate"] > 0.01 * rep["n_faces"]:
        return "FAIL: bad or degenerate triangles"
    if rep.get("n_edges_open"):
        return "WARN: open / non-watertight (ray casting may leak)"
    if rep.get("n_edges_nonmanifold"):
        return "OK (assembly seams: non-manifold shared edges, ray-casting safe)"
    if rep["normals_mixed"]:
        return "WARN: mixed face winding (cleaning uses |cos|, safe)"
    return "OK"


def print_report(rep: dict):
    print(f"  mesh check : {rep['score']}")
    print(f"    faces={rep['n_faces']}, points={rep['n_points']}, "
          f"degenerate={rep['n_degenerate']}, bad-points={rep['n_bad_points']}")
    print(f"    bbox (mm) x:[{rep['bbox_mm'][0]:.1f},{rep['bbox_mm'][3]:.1f}] "
          f"y:[{rep['bbox_mm'][1]:.1f},{rep['bbox_mm'][4]:.1f}] "
          f"z:[{rep['bbox_mm'][2]:.1f},{rep['bbox_mm'][5]:.1f}]")
    print(f"    edges (dedup): total={rep['n_edges']}, open={rep['n_edges_open']}, "
          f"non-manifold={rep['n_edges_nonmanifold']} -> watertight={rep['watertight']}")
    print(f"    inward faces: {rep['n_faces_inward']}/{rep['n_faces']} "
          f"(normals mixed: {rep['normals_mixed']})")
