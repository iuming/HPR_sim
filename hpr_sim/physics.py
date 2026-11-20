# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""
HPR cleaning physics.

Implements the reference models:

  * cleaning intensity = cos(phi), phi the angle between the jet and the
    surface normal (CN121302684A).  We use the absolute value so winding of
    the imported STL does not matter.
  * jet velocity vs travel distance:
      - "droplet" (default): intact core up to the breakup length, then
        droplets decelerated by aerodynamic drag (WEPS143 / Fermilab style).
        v(d) = v0                    for d <= L_breakup
        v(d) = v0 * exp(-k (d-L_b))  for d >  L_breakup
        with L_breakup = breakup_coef*dN*sqrt(rho_w/rho_air),
             k = 3*rho_air*C_d/(8*rho_w*r_drop),
             r_drop from critical-Weber secondary breakup.
      - "linear": simple placeholder, v(d) = v0*max(0,1-d/v_lin_decay_mm).
  * impact force / impact stagnation pressure from momentum flux
    (WEPS143, FERMILAB-SLIDES-25-0209).
  * Weber number at the wall (atomisation scale, FERMILAB-SLIDES-25-0209).
"""
from __future__ import annotations

import numpy as np

from .config import JetParams, WaterParams


def droplet_diameter_m(v0: float, jp: JetParams, water: WaterParams) -> float:
    """Mean intact droplet diameter [m] from critical-Weber secondary breakup,
    d_drop = We_c * sigma / (rho_air * v0^2)."""
    return (jp.drop_we_critical * water.surface_tension
            / (water.rho_air * max(v0, 1e-9) ** 2))


def breakup_length_m(v0: float, jp: JetParams, water: WaterParams) -> float:
    """Intact (core) jet length [m] before atomization."""
    return (jp.breakup_coef * (jp.nozzle_diameter_mm * 1e-3)
            * (water.rho / water.rho_air) ** 0.5)


def droplet_reynolds(v_ms: float, d_drop_m: float, water: WaterParams):
    """Droplet Reynolds number Re = rho_air * v * d_drop / mu_air."""
    return water.rho_air * np.abs(v_ms) * d_drop_m / water.mu_air


def schiller_naumann_cd(Re):
    """Schiller-Naumann drag correlation for a spherical droplet in air:
    Re < 1000 : Cd = (24/Re) * (1 + 0.15 Re^0.687)
    Re >= 1000: Cd = 0.44
    """
    Re = np.maximum(np.abs(Re), 1e-9)
    return np.where(Re < 1000.0,
                    24.0 / Re * (1.0 + 0.15 * Re ** 0.687),
                    0.44)


# Morsi & Alexander (1972) piecewise coefficients: Cd = a1 + a2/Re + a3/Re^2
# for each Re bin (lower bound, a1, a2, a3).
_MORSI_EDGES = np.array([0.0, 0.1, 1.0, 10.0, 100.0, 1e3, 5e3, 1e4, 5e4])
_MORSI_COEF = np.array([
    (0.0, 24.0, 0.0),             # Re < 0.1 (Stokes)
    (3.69, 22.73, 0.0903),        # 0.1  <= Re < 1
    (1.222, 29.1667, -3.8889),    # 1    <= Re < 10
    (0.6167, 46.5, -116.67),      # 10   <= Re < 100
    (0.3644, 98.33, -2778.0),     # 100  <= Re < 1000
    (0.357, 148.62, -47500.0),    # 1e3  <= Re < 5e3
    (0.46, -490.546, 578700.0),   # 5e3  <= Re < 1e4
    (0.5191, -1662.5, 5416700.0), # 1e4  <= Re < 5e4
    (0.44, 0.0, 0.0),             # Re >= 5e4
])


def morsi_alexander_cd(Re):
    """Morsi & Alexander (1972) drag correlation:
    Cd = a1 + a2/Re + a3/Re^2 with piecewise constants."""
    Re = np.maximum(np.abs(Re), 1e-9)
    row = np.searchsorted(_MORSI_EDGES, Re, side="right") - 1
    row = np.clip(row, 0, len(_MORSI_COEF) - 1)
    a1 = _MORSI_COEF[row, 0]
    a2 = _MORSI_COEF[row, 1]
    a3 = _MORSI_COEF[row, 2]
    return a1 + a2 / Re + a3 / (Re * Re)


def putnam_cd(Re):
    """Putnam (1961) correlation for deforming droplets in sprays:
    Re < 1000 : Cd = (24/Re) * (1 + Re^(2/3)/6)
    Re >= 1000: Cd = 0.44
    """
    Re = np.maximum(np.abs(Re), 1e-9)
    return np.where(Re < 1000.0, 24.0 / Re * (1.0 + Re ** (2.0 / 3.0) / 6.0), 0.44)


def drag_coefficient(Re, mode: str, cd_const: float = 0.5):
    """Dispatch to the requested drag law."""
    if mode == "schiller-naumann":
        return schiller_naumann_cd(Re)
    if mode == "morsi-alexander":
        return morsi_alexander_cd(Re)
    if mode == "putnam":
        return putnam_cd(Re)
    if mode == "constant":
        return np.full_like(np.asarray(Re, dtype=float), cd_const)
    raise ValueError(f"unknown drop_cd_mode: {mode}")


def droplet_decay_coeff(v0: float, jp: JetParams, water: WaterParams) -> float:
    """Constant-`drop_cd` velocity decay coefficient k [1/m] of the
    exponential droplet model (only used when drop_cd_mode == "constant")."""
    r_d = 0.5 * droplet_diameter_m(v0, jp, water)
    return 3.0 * water.rho_air * jp.drop_cd / (8.0 * water.rho * max(r_d, 1e-12))


def _drag_slope(v, jp: JetParams, water: WaterParams, d_drop_m: float) -> float:
    """dv/dx = - (3 rho_air Cd(Re) / (8 rho_w r_d)) * v  [1/s], x in m."""
    r_d = 0.5 * d_drop_m
    if jp.drop_cd_mode == "constant":
        cd = jp.drop_cd
    else:
        Re = droplet_reynolds(v, d_drop_m, water)
        cd = drag_coefficient(Re, jp.drop_cd_mode, jp.drop_cd)
    return 3.0 * water.rho_air * cd / (8.0 * water.rho * max(r_d, 1e-12))


def velocity_table(v0: float, jp: JetParams, water: WaterParams,
                   dmax_mm: float = 600.0, step_mm: float = 0.1):
    """Jet velocity as a function of travel distance for the configured model.

    For the droplet model with a constant C_D this is the closed form
    v0*exp(-k(x-Lb)); for the Re-dependent (Schiller-Naumann) C_D the ODE
    dv/dx = -(3 rho_air Cd(Re(v)) / (8 rho_w r_d)) v is integrated numerically
    with an RK2 (midpoint) scheme from the breakup length.

    Returns (d_mm, v_ms) 1-D numpy arrays.
    """
    d = np.arange(0.0, dmax_mm + step_mm, step_mm)
    if jp.velocity_model == "linear":
        f = np.maximum(0.0, 1.0 - d / max(jp.v_lin_decay_mm, 1e-9))
        return d, v0 * f

    Lb = breakup_length_m(v0, jp, water)
    d_d = droplet_diameter_m(v0, jp, water)
    if jp.drop_cd_mode == "constant":
        k = droplet_decay_coeff(v0, jp, water)
        v = v0 * np.exp(-k * np.maximum(d * 1e-3 - Lb, 0.0))
        return d, v

    # --- Schiller-Naumann: numerical integration (RK2 midpoint) -----------
    v = np.empty_like(d)
    v.fill(v0)
    if dmax_mm * 1e-3 <= Lb:
        return d, v
    start = int(np.ceil(Lb * 1e3 / step_mm))
    x = Lb              # m
    vcur = v0           # m/s
    h = step_mm * 1e-3  # m
    n_after = len(d) - start
    if n_after <= 0:
        return d, v
    for i in range(n_after):
        # slope at current point
        k1 = -_drag_slope(vcur, jp, water, d_d) * vcur
        vmid = vcur + 0.5 * h * k1
        k2 = -_drag_slope(vmid, jp, water, d_d) * vmid
        vcur = vcur + h * k2
        v[start + i] = vcur
        x += h
        if vcur < 1e-3:
            v[start + i:] = 0.0
            break
    return d, v


def efold_distance_mm(v0: float, jp: JetParams, water: WaterParams,
                      dmax_mm: float = 600.0) -> float:
    """Distance (past the breakup length) at which v drops to v0/e."""
    d, v = velocity_table(v0, jp, water, dmax_mm=dmax_mm)
    Lb = breakup_length_m(v0, jp, water) * 1e3
    th = v0 / np.e
    idx = np.nonzero(v <= th)[0]
    if len(idx) == 0:
        return float("nan")
    i = idx[0]
    d0 = d[max(0, i - 1)]
    v0_ = v[max(0, i - 1)]
    v1_ = v[i]
    denom = v1_ - v0_
    frac = (th - v0_) / denom if abs(denom) > 1e-12 else 0.0
    return float(d0 + frac * (d[i] - d0) - Lb)


def jet_velocity_at_distance(v0: float, dist_mm, jp: JetParams,
                             water: WaterParams):
    """Jet velocity [m/s] at a distance [mm] from the nozzle.

    Accepts a scalar or an array of distances.  Uses the velocity table built
    for the configured model + drag mode.
    """
    d_tab, v_tab = velocity_table(v0, jp, water)
    d = np.asarray(dist_mm, dtype=float)
    return np.interp(d, d_tab, v_tab)


def impact_force(flow_m3s: float, vel_ms: float, rho: float) -> float:
    """Impact force [N] = rho * Q * v (momentum flux)."""
    return rho * flow_m3s * vel_ms


def impact_pressure(vel_ms: float, rho: float) -> float:
    """Stagnation pressure [Pa] = 0.5 rho v^2."""
    return 0.5 * rho * vel_ms ** 2


def weber_number(vel_ms: float, nozzle_d_mm: float, water: WaterParams) -> float:
    """Weber number at the wall (air-side, water droplets)."""
    return water.rho_air * vel_ms ** 2 * (nozzle_d_mm * 1e-3) / water.surface_tension


def intensity(hit_cos: np.ndarray) -> np.ndarray:
    """Cleaning intensity per ray: |cos(phi)| in [0,1]."""
    return np.abs(hit_cos)


def cleaning_metric(intens: np.ndarray, dist_mm: np.ndarray,
                    p_impact_bar: np.ndarray, sp, d_max_eff: float = None) -> np.ndarray:
    """Effective cleaning flag:
      intensity above threshold AND usable nozzle-surface distance window AND
      (optionally) minimum wall impact pressure.
    `d_max_eff` overrides the distance limit (e.g. from a Guha-style
    stand-off calibration d_max = kd * nozzle_diameter).
    """
    dmax = d_max_eff if d_max_eff is not None else sp.d_max_mm
    ok = ((intens >= sp.clean_cos_thresh)
          & (dist_mm >= sp.d_min_mm) & (dist_mm <= dmax))
    if getattr(sp, "min_impact_bar", 0.0) > 0:
        ok &= p_impact_bar >= sp.min_impact_bar
    return ok
