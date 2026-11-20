# -*- coding: utf-8 -*-
"""
Main HPR simulation driver.

Casts the nozzle ray bundles from every pose of the wand trajectory against
the cavity mesh, accumulates per-triangle metrics (hit count, effective
cleaning, exposure, intensity, impact) and produces a results dictionary.
"""
from __future__ import annotations

import time

import numpy as np

from .config import (
    JetParams,
    MotionParams,
    NozzleParams,
    SimulationParams,
    WaterParams,
)
from .mesh import TriangleMesh
from .motion import build_poses
from .nozzle import generate_rays
from .physics import (
    cleaning_metric,
    drag_coefficient,
    droplet_diameter_m,
    droplet_reynolds,
    efold_distance_mm,
    impact_pressure,
    intensity,
    breakup_length_m,
    velocity_table,
    weber_number,
)
from .raycast import cone_of_dirs, intersect_bundle_fast

_AXIS_IDX = {"x": 0, "y": 1, "z": 2}


class SimResults:
    """Container for per-triangle simulation results."""

    def __init__(self, n_tri: int):
        self.n_hits = np.zeros(n_tri, np.float64)
        self.n_effective = np.zeros(n_tri, np.float64)
        self.exposure_s = np.zeros(n_tri, np.float64)
        self.max_intensity = np.zeros(n_tri, np.float64)
        self.max_impact_bar = np.zeros(n_tri, np.float64)
        self.min_dist = np.full(n_tri, np.inf)
        self.first_hit_time = np.full(n_tri, np.inf)

    def accumulate_ray_hits(self, tri_ids, effective, dt, runtime_s,
                            intens, dist_mm, p_impact_bar):
        hit, = np.nonzero(np.asarray(tri_ids) >= 0)
        if len(hit) == 0:
            return
        ids = tri_ids[hit]
        eff = effective[hit]
        np.add.at(self.n_hits, ids, 1)
        np.add.at(self.n_effective, ids, eff)
        np.add.at(self.exposure_s, ids, dt * eff)
        # first hit time
        free = np.isinf(self.first_hit_time[ids])
        self.first_hit_time[ids[free]] = runtime_s
        np.maximum.at(self.max_intensity, ids, intens[hit])
        np.maximum.at(self.max_impact_bar, ids, p_impact_bar[hit])
        np.minimum.at(self.min_dist, ids, dist_mm[hit])


def run_simulation(
    mesh: TriangleMesh,
    water: WaterParams,
    jet: JetParams,
    nozzle: NozzleParams,
    motion: MotionParams,
    sim: SimulationParams,
    axis: str = "x",
    progress=True,
) -> SimResults:
    """Run the HPR ray-casting simulation over the mesh.

    Returns a SimResults instance.
    """
    ax = _AXIS_IDX[axis]
    t0 = time.time()

    # wand offset from the axis is zero: wand sits on the beam axis
    origin_off = np.zeros(3)

    v0 = jet.exit_velocity(water)
    flow = jet.jet_flow(water)
    # velocity-vs-distance lookup table (built once per run)
    d_tab, v_tab = velocity_table(v0, jet, water)
    # Guha et al. (2010) stand-off calibration: d_max = kd * nozzle diameter
    d_max_eff = (sim.dmax_kd * jet.nozzle_diameter_mm if sim.dmax_kd > 0
                 else sim.d_max_mm)

    poses = build_poses(motion, sim.dt_s, sim.step_trans_mm, sim.step_rot_deg)
    rng = np.random.default_rng(sim.seed)

    results = SimResults(mesh.n_tri)

    step = 0
    total = len(poses)
    report_every = max(1, total // 20)
    for pose in poses:
        O = origin_off.copy()
        O[ax] = pose.x
        if motion.axis not in (axis,):
            O[_AXIS_IDX[motion.axis]] = 0.0
        origin = np.array(O, dtype=float)

        dirs = generate_rays(origin, pose.head_angle, nozzle.head_diameter_mm,
                             nozzle, axis=axis, phase=0.0, rng=rng)
        if len(dirs) == 0:
            continue

        # pre-filter candidate triangles (cone + sphere conservative) and run
        # the vectorised per-ray Möller-Trumbore intersection
        cone_axis, cone_half = cone_of_dirs(dirs)
        res = intersect_bundle_fast(origin, dirs, mesh, cone_axis, cone_half,
                                    sim.ray_len_max_mm, sub_cone_deg=10.0)

        hit = res["hit"]
        if not hit.any():
            step += 1
            continue

        t = res["t"][hit]
        ids = res["tri"][hit]
        intens = intensity(res["cos_phi"][hit])
        dist_mm = np.clip(t, 0, None)

        # wall impact pressure with the configured velocity-attenuation model
        v_wall = np.interp(dist_mm, d_tab, v_tab)                  # m/s
        p_impact = impact_pressure(v_wall, water.rho) / 1e5        # bar

        eff = cleaning_metric(intens, dist_mm, p_impact, sim, d_max_eff=d_max_eff)

        results.accumulate_ray_hits(ids, eff, sim.dt_s, pose.t,
                                    intens, dist_mm, p_impact)

        step += 1
        if progress and step % report_every == 0:
            msg = (f"  HPR step {step}/{total} (x={pose.x:7.1f} mm, "
                   f"hits/step={int(hit.sum())}) elapsed={time.time()-t0:.1f}s")
            print(msg, flush=True)

    results.n_poses = total
    results.wall_time_s = time.time() - t0
    results.v0_ms = v0
    results.flow_m3s = flow
    # example Weber number at the exit
    results.weber_exit = weber_number(v0, jet.nozzle_diameter_mm, water)
    # droplet-model diagnostics
    if jet.velocity_model == "droplet":
        d_drop = droplet_diameter_m(v0, jet, water) * 1e6
        Lb = breakup_length_m(v0, jet, water) * 1e3
        Re0 = droplet_reynolds(v0, d_drop * 1e-6, water)
        Cd0 = float(drag_coefficient(np.asarray(Re0), jet.drop_cd_mode,
                                     jet.drop_cd))
        ef = efold_distance_mm(v0, jet, water)
        results.d_drop_um = float(d_drop)
        results.breakup_len_mm = float(Lb)
        results.re_exit = float(Re0)
        results.cd_exit = Cd0
        results.efold_mm = float(ef)
        results.drop_cd_mode = jet.drop_cd_mode
        if sim.dmax_kd > 0:
            results.dmax_eff_mm = d_max_eff
        print(f"  HPR done: {total} poses in {results.wall_time_s:.1f}s "
              f"(v0={v0:.1f} m/s, Q/jet={flow*1e6:.2f} L/min, We={results.weber_exit:.0f}), "
              f"droplet:d={d_drop:.0f} um, L_break={Lb:.0f} mm, "
              f"Re0={Re0:.0f}, Cd0={Cd0:.2f}, 1/e={ef:.0f} mm "
              f"(cd={jet.drop_cd_mode}, dmax={d_max_eff:.0f} mm)", flush=True)
    else:
        print(f"  HPR done: {total} poses in {results.wall_time_s:.1f}s "
              f"(v0={v0:.1f} m/s, Q/jet={flow*1e6:.2f} L/min, We={results.weber_exit:.0f}), "
              f"model={jet.velocity_model}", flush=True)
    return results
