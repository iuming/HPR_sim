# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""
HPR nozzle-ray generation.

The patent CN121302684A and the FRIB/RAON tools model the HPR nozzle as a set
of rays emitted from a common origin with different direction vectors.  Real
heads have a ring of holes; each hole produces a jet that we sample as a small
cone of rays.  The head rotates around the wand axis, so a full 360-degree
ring of holes is swept as the head turns.
"""
from __future__ import annotations

import numpy as np

from .config import NozzleParams

_AXIS = {"x": 0, "y": 1, "z": 2}


def hole_direction(azimuth: float, tilt_rad: float, axis: str = "x"):
    """Direction of the hole axis at a given azimuth around the wand axis.

    azimuth : angle (rad) around the axis, azimuth=0 -> +y
    tilt_rad : tilt out of the radial plane; + tilt points toward +axis
    """
    r = np.cos(azimuth)
    s = np.sin(azimuth)
    radial = np.zeros(3)
    radial[_AXIS[axis]] = 1.0
    # build orthonormal basis: axis, u1 (at azimuth 0), u2
    if axis == "x":
        u1 = np.array([0.0, np.cos(azimuth), np.sin(azimuth)])   # radial dir
    elif axis == "y":
        u1 = np.array([np.cos(azimuth), 0.0, np.sin(azimuth)])
    else:
        u1 = np.array([np.cos(azimuth), np.sin(azimuth), 0.0])
    # radial component direction along the axis is the tilt
    return np.sin(tilt_rad) * radial + np.cos(tilt_rad) * u1


def _rot_axis_angle(v, k, ang):
    """Rodrigues rotation of v around unit axis k by ang."""
    c, s = np.cos(ang), np.sin(ang)
    return v * c + np.cross(k, v) * s + k * np.dot(k, v) * (1.0 - c)


def generate_rays(
    origin: np.ndarray,
    head_azimuth: float,
    head_diameter_mm: float,
    params: NozzleParams,
    axis: str = "x",
    phase: float = 0.0,
    rng: np.random.Generator = None,
):
    """Generates a bundle of ray directions from a common origin.

    origin : (3,) nozzle position (wand head centre)
    head_azimuth : current head rotation angle (rad)
    phase : extra azimuthal phase per pass (rad)

    Returns (dirs (R,3) normalized, cones ((R,3) hole-axis unit vectors)).
    """
    if rng is None:
        rng = np.random.default_rng(0)

    tilt_fwd = np.deg2rad(params.tilt_fwd_deg)
    tilt_bwd = np.deg2rad(params.tilt_bwd_deg)
    jet = np.deg2rad(params.jet_angle_deg)

    # holes are arranged around the head; half point forward, half backward
    # (typical HPR +25/-25 deg design)
    dirs, hole_ax = [], []
    for i in range(params.n_holes):
        az = 2.0 * np.pi * i / params.n_holes + head_azimuth + phase
        is_fwd = (i % 2) == 0
        tilt = tilt_fwd if is_fwd else tilt_bwd
        haxe = hole_direction(az, tilt, axis)
        # sample a cone of rays around the hole axis
        for k in range(params.n_spray):
            # deterministic low-discrepancy sampling of the cone
            phi = np.pi * (0.5 + k) / params.n_spray     # around axis
            beta = (2.0 * np.pi) * k / params.n_spray    # azimuthal
            # sample half a disk on the sphere cap
            theta = jet * np.sqrt(phi / np.pi)           # polar (concentrated)
            # build a basis around hole axis
            ref = np.array([0.0, 0.0, 1.0])
            if abs(np.dot(haxe, ref)) > 0.99:
                ref = np.array([1.0, 0.0, 0.0])
            u = np.cross(haxe, ref)
            u /= np.linalg.norm(u)
            w = np.cross(haxe, u)
            d = haxe * np.cos(theta) + (u * np.cos(beta) + w * np.sin(beta)) * np.sin(theta)
            dirs.append(d / np.linalg.norm(d))
        hole_ax.append(haxe)
    dirs = np.asarray(dirs)
    return dirs


def ray_origins(origin: np.ndarray, head_diameter_mm: float, params: NozzleParams):
    """Return common ray origin(s).  Patent models use a single shared origin,
    so we return just the head centre."""
    return origin[None, :]
