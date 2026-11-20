# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""Tests for the mesh quality checks and the geometry auto-detection."""
import numpy as np
import pytest

from hpr_sim.mesh import TriangleMesh
from hpr_sim.meshcheck import check_mesh
from hpr_sim.geometry import detect_beam_axis, propose_scan_range

from test_raycast import _unit_cube_mesh


def test_watertight_cube():
    rep = check_mesh(_unit_cube_mesh())
    assert rep["watertight"] is True
    assert rep["n_edges_open"] == 0
    assert rep["n_edges_nonmanifold"] == 0
    assert rep["n_degenerate"] == 0


def test_open_box_detected():
    """A cube missing its top face -> open edges."""
    v = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
    tri = np.array([[0, 1, 2], [0, 2, 3],      # bottom
                    [0, 4, 5], [0, 5, 1],      # x=0
                    [3, 2, 6], [3, 6, 7],      # x=1
                    [0, 3, 7], [0, 7, 4],      # y=0
                    [1, 5, 6], [1, 6, 2]], dtype=np.int64)   # y=1 (top missing)
    m = TriangleMesh(v, tri)
    rep = check_mesh(m)
    assert rep["watertight"] is False
    assert rep["n_edges_open"] >= 4


def test_beam_axis_and_scan_proposal():
    # a long tube along x (4 triplets of a ring) - axis should be x
    L = 200.0
    R = 50.0
    theta = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    xs = np.linspace(-L / 2, L / 2, 24)
    points = []
    for x in xs:
        for t in theta:
            points.append([x, R * np.cos(t), R * np.sin(t)])
    # close the ends with discs (rings to centre) - simple fans
    verts = np.array(points, dtype=float)
    m = TriangleMesh(verts, np.array([[0, 1, 2]]), name="tube")
    # (triangle soup is fine for the axis/radial-profile logic)
    assert detect_beam_axis(m) == "x"
    prop = propose_scan_range(m, "x", n_bins=60)
    if prop is not None:
        assert prop["x_start"] < prop["x_end"]
        assert prop["equator_radius_mm"] > 0


def test_check_degenerate_detected():
    v = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=float)
    m = TriangleMesh(v, np.array([[0, 1, 2]]), name="deg")
    rep = check_mesh(m)
    assert rep["n_degenerate"] >= 1
