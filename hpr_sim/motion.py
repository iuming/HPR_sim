# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""
Wand motion: axial translation + rotation + reciprocation.

The wand sits on the beam axis and translates from x_start to x_end and back
(repeated `n_recip` times), while the head rotates continuously.  The motion
is discretised into emission steps bounded by `step_trans_mm` (axial) and
`step_rot_deg` (rotation).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import MotionParams


@dataclass
class Pose:
    x: float            # axial coordinate [mm]
    head_angle: float   # head rotation angle [rad]
    t: float            # time [s]
    pass_idx: int


def build_poses(mp: MotionParams, dt_s: float, step_trans_mm: float,
                step_rot_deg: float) -> list:
    """Discretized wand trajectory.

    Steps are emitted such that the axial step and the rotation step both stay
    under the requested bounds.  Each emission event also gets the physical
    time `dt_s`.
    """
    D = mp.x_end_mm - mp.x_start_mm
    if D <= 0:
        raise ValueError("x_end must be > x_start")

    v = mp.trans_speed_mm_s
    omega = 2.0 * np.pi * mp.rot_speed_rpm / 60.0          # rad/s

    # physical stepping bounded by axial step, rotation step and dt_s
    dt_mech = min(step_trans_mm / max(v, 1e-9),
                  np.deg2rad(step_rot_deg) / max(omega, 1e-9))
    dt = dt_mech
    if dt_s is not None:
        dt = min(dt, dt_s)

    poses = []
    t = 0.0
    phase_sweep = np.deg2rad(mp.phase_shift_pass_deg)
    for p in range(mp.n_recip):
        phase = phase_sweep * p
        # forward
        n_fwd = max(1, int(np.ceil(D / (dt * v))))
        xs_fwd = np.linspace(mp.x_start_mm, mp.x_end_mm, n_fwd + 1)[:-1]
        for x in xs_fwd:
            theta = phase + omega * t
            poses.append(Pose(x, theta, t, p))
            t += dt
        # dwell at the end
        if mp.end_dwell_s > 0:
            for _ in range(int(np.ceil(mp.end_dwell_s / dt))):
                theta = phase + omega * t
                poses.append(Pose(mp.x_end_mm, theta, t, p))
                t += dt
        if p == mp.n_recip - 1:
            break
        # backward
        n_bwd = max(1, int(np.ceil(D / (dt * v))))
        xs_bwd = np.linspace(mp.x_end_mm, mp.x_start_mm, n_bwd)[1:]
        for x in xs_bwd:
            theta = phase + omega * t
            poses.append(Pose(x, theta, t, p))
            t += dt
        if mp.end_dwell_s > 0:
            for _ in range(int(np.ceil(mp.end_dwell_s / dt))):
                theta = phase + omega * t
                poses.append(Pose(mp.x_start_mm, theta, t, p))
                t += dt
    return poses
