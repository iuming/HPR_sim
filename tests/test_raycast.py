# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""Tests for mesh I/O, STL round-trip and the batch ray caster."""
import os
import pathlib

import numpy as np
import pytest

from hpr_sim.mesh import TriangleMesh
from hpr_sim.raycast import intersect_bundle_fast, cone_of_dirs


def _unit_cube_mesh(center=(0, 0, 0)):
    """12-triangle closed unit cube centred at `center` (0..1 edge)."""
    c = np.array(center)
    lo, hi = c, c + 1.0
    v = np.array([
        [lo[0], lo[1], lo[2]], [hi[0], lo[1], lo[2]],
        [hi[0], hi[1], lo[2]], [lo[0], hi[1], lo[2]],
        [lo[0], lo[1], hi[2]], [hi[0], lo[1], hi[2]],
        [hi[0], hi[1], hi[2]], [lo[0], hi[1], hi[2]],
    ], dtype=float)
    tri = np.array([
        [0, 1, 2], [0, 2, 3],   # bottom (z=0), outward -z
        [7, 6, 5], [7, 5, 4],   # top (z=1), outward +z
        [0, 4, 5], [0, 5, 1],   # x=0 face
        [3, 2, 6], [3, 6, 7],   # x=1 face
        [0, 3, 7], [0, 7, 4],   # y=0 face
        [1, 5, 6], [1, 6, 2],   # y=1 face
    ], dtype=np.int64)
    return TriangleMesh(v, tri, name="cube")


def test_stl_roundtrip(ws_results):
    mesh = _unit_cube_mesh()
    p = str(pathlib.Path(ws_results) / "_tmp_cube_rnd.stl")
    mesh.write_stl(p)
    try:
        again = TriangleMesh.from_stl(p)
    finally:
        if os.path.exists(p):
            os.remove(p)
    assert again.n_tri == mesh.n_tri
    # STL stores per-face vertices (deduplicated); compare the unique corners
    a = np.unique(np.round(mesh.vertices, 6), axis=0)
    b = np.unique(np.round(again.vertices, 6), axis=0)
    assert a.shape == b.shape
    assert np.allclose(np.sort(a, axis=0), np.sort(b, axis=0), atol=1e-5)


def test_mesh_geometry():
    m = _unit_cube_mesh()
    assert m.total_area == pytest.approx(6.0, rel=1e-6)
    assert m.centroids.shape == (12, 3)
    assert (np.isfinite(m.normals_unit)).all()


def test_mt_ray_hits_triangle():
    # triangle in z=0 plane
    v = np.array([[0, 0, 0], [1, 0, 0.0], [0.0, 1, 0.0]], dtype=float)
    m = TriangleMesh(v, np.array([[0, 1, 2]]), name="tri")
    origin = np.array([0.25, 0.25, 1.0])
    dirs = np.array([[0, 0, -1.0]])
    axis, half = cone_of_dirs(dirs)
    res = intersect_bundle_fast(origin, dirs, m, axis, half, ray_len=10.0)
    assert res["hit"][0]
    assert res["tri"][0] == 0
    assert res["t"][0] == pytest.approx(1.0, rel=1e-6)
    assert res["cos_phi"][0] == pytest.approx(1.0)   # perpendicular


def test_mt_miss():
    v = np.array([[0, 0, 0], [1, 0, 0.0], [0.0, 1, 0.0]], dtype=float)
    m = TriangleMesh(v, np.array([[0, 1, 2]]), name="tri")
    origin = np.array([2.0, 2.0, 1.0])   # misses the triangle
    dirs = np.array([[0, 0, -1.0]])
    axis, half = cone_of_dirs(dirs)
    res = intersect_bundle_fast(origin, dirs, m, axis, half, ray_len=10.0)
    assert not res["hit"][0]


def test_bundle_inside_cube():
    m = _unit_cube_mesh()
    origin = np.array([0.5, 0.5, 0.5])
    dirs = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0.0],
                     [0, 0, 1], [0, 0, -1]], dtype=float)
    axis, half = cone_of_dirs(dirs)
    res = intersect_bundle_fast(origin, dirs, m, axis, half, ray_len=10.0,
                                sub_cone_deg=10.0)
    assert res["hit"].all()
    # each ray should hit a wall at distance 0.5
    assert np.allclose(res["t"], 0.5, atol=1e-8)
