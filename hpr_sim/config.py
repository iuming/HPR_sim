# -*- coding: utf-8 -*-
"""
Configuration / parameters for the HPR simulation.

All lengths in millimetres, times in seconds, pressures in bar, velocities in
m/s unless stated otherwise. Defaults follow the reference material:

  * pump pressure      ~100 bar,  Cd ~0.7,  nozzle diameter 0.4-0.6 mm
    (FERMILAB-SLIDES-25-0209)
  * exit velocity      ~100-140 m/s  (Bernoulli with Cd; ~140 m/s at 100 bar
    reported by Fermilab)
  * cleaning intensity defined as cos(phi) between jet and surface normal
    (CN121302684A); effective when above a threshold
  * distance / angle / exposure are the dominant coverage parameters
    (MOPB070 / PIP-II)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List


@dataclass
class WaterParams:
    """Water / process fluid properties."""
    rho: float = 1000.0            # kg/m^3
    surface_tension: float = 0.073 # N/m  (water-air interface)
    rho_air: float = 1.2           # kg/m^3
    mu_air: float = 1.81e-5        # Pa s  (air dynamic viscosity)

    @property
    def rho_mm3(self):             # kg / mm^3
        return self.rho * 1e-9


@dataclass
class JetParams:
    """High pressure jet / nozzle physics.

    Two velocity-attenuation models are available for the jet as it travels:

    ``velocity_model = "linear"``  (previous, simple placeholder)
        v(d) = v0 * max(0, 1 - d / v_lin_decay_mm)

    ``velocity_model = "droplet"`` (refined, WEPS143-style droplet drag)
        - intact core region   0 <= d <= L_breakup : v = v0
          L_breakup = breakup_coef * dN * sqrt(rho_w / rho_air)
        - droplet region       d > L_breakup :
          v(d) = v0 * exp(-k * (d - L_breakup))       [mm/s exponent made from m]
          with  k = 3 * rho_air * drop_cd / (8 * rho_w * r_droplet)
          and droplet diameter from critical-Weber secondary breakup:
          d_droplet = drop_we_critical * sigma / (rho_air * v0^2)
    """
    pressure_bar: float = 100.0      # pump pressure
    nozzle_diameter_mm: float = 0.5  # per-hole diameter
    cd: float = 0.70                 # discharge coefficient
    velocity_model: str = "droplet"  # "linear" | "droplet"
    # --- linear model (kept for comparison / backwards compatibility) ---
    v_lin_decay_mm: float = 250.0    # velocity -> 0 at this distance
    # --- droplet-drag model parameters ---
    # drag coefficient law:
    #   "schiller-naumann"  classic rigid-sphere correlation (Schiller & Naumann 1933)
    #   "morsi-alexander"   piecewise drag from Morsi & Alexander (1972), JFM
    #   "putnam"            deformation-aware correlation (Putnam 1961)
    #   "constant"          fixed `drop_cd`
    drop_cd_mode: str = "schiller-naumann"
    # constant drag coefficient used when drop_cd_mode == "constant"
    drop_cd: float = 0.50            # droplet drag coefficient (sphere, Re~500)
    drop_we_critical: float = 12.0   # critical Weber number for secondary breakup
    breakup_coef: float = 4.5        # L_breakup = coef * dN * sqrt(rho_w/rho_air)
    # if True use pure Bernoulli exit velocity; else let user override
    use_bernoulli: bool = True
    velocity_override_ms: float = 0.0

    @property
    def area_m2(self):
        import math
        return math.pi * (self.nozzle_diameter_mm * 1e-3) ** 2 / 4.0

    def exit_velocity(self, water: WaterParams) -> float:
        """Jet exit velocity in m/s (Bernoulli) or the user override."""
        if not self.use_bernoulli and self.velocity_override_ms > 0:
            return self.velocity_override_ms
        dp = self.pressure_bar * 1e5                  # Pa
        return self.cd * (2.0 * dp / water.rho) ** 0.5

    def jet_flow(self, water: WaterParams) -> float:
        """Volumetric flow per jet in m^3/s."""
        dp = self.pressure_bar * 1e5
        return self.cd * self.area_m2 * (2.0 * dp / water.rho) ** 0.5


@dataclass
class NozzleParams:
    """Wand / nozzle head layout (holes around the head)."""
    n_holes: int = 6                      # holes distributed in azimuth
    tilt_deg: float = 25.0                # tilt of hole axis from the radial plane
    # tilt sign convention: +tilt aims holes toward the +axis direction,
    # -tilt toward the -axis direction; a list allows mixed rings.
    tilt_fwd_deg: float = 25.0
    tilt_bwd_deg: float = -25.0
    n_spray: int = 11                     # rays used to sample each jet cone
    jet_angle_deg: float = 5.0            # spray cone half-angle at each hole
    head_diameter_mm: float = 25.0        # wand head outer diameter


@dataclass
class MotionParams:
    """Wand motion: axial translation with rotation + reciprocation."""
    axis: str = "x"                       # beam axis coordinate
    x_start_mm: float = -160.0
    x_end_mm: float = 160.0
    trans_speed_mm_s: float = 10.0        # wand translation speed
    rot_speed_rpm: float = 24.0           # head rotation speed
    n_recip: int = 2                      # reciprocation count
    # azimuth phase shifted between passes (deg); 0 keeps passes identical
    phase_shift_pass_deg: float = 0.0
    # dwell: keep head at each end for this many seconds before reversing
    end_dwell_s: float = 0.0


@dataclass
class SimulationParams:
    """Numerical stepping and cleaning criteria."""
    step_trans_mm: float = 1.0     # axial movement between emission steps
    step_rot_deg: float = 3.0      # rotation between emission steps
    dt_s: float = 0.1              # physical time between emission steps
    # cleaning criteria (CN121302684A)
    clean_cos_thresh: float = 0.5  # intensity (cos of incidence) threshold
    d_min_mm: float = 5.0          # too close -> damaged / not cleaning
    d_max_mm: float = 200.0        # too far -> jet too weak
    # Guha et al. (2010) stand-off model: if >0, the effective maximum
    # nozzle-wall distance becomes  d_max = dmax_kd * nozzle_diameter_mm
    # (their measured critical stand-off is ~26 D, optimal ~5 D).
    dmax_kd: float = 0.0           # 0 -> use d_max_mm directly
    # optional: effective cleaning additionally requires the wall impact
    # stagnation pressure (0.5 rho v^2 with the droplet/decay model) >= this
    # value in bar.  0.0 disables the criterion.
    min_impact_bar: float = 0.0
    ray_len_max_mm: float = 300.0  # max ray travel considered
    # if >0, rotate nozzle orientation bias by this angle vs radial (optimization)
    wand_bias_deg: float = 0.0

    seed: int = 0                  # for phase jitter (0 = deterministic)


def default_nozzle() -> NozzleParams:
    return NozzleParams()


# ----------------------------------------------------------------------
# Generic dataclass <-> dict conversion (for the YAML/JSON configuration)
# ----------------------------------------------------------------------
from dataclasses import fields, is_dataclass


def to_dict(obj) -> dict:
    """Convert a (nested) dataclass into a plain dict."""
    if is_dataclass(obj):
        return {f.name: to_dict(getattr(obj, f.name)) for f in fields(obj)}
    if isinstance(obj, (list, tuple)):
        return [to_dict(x) for x in obj]
    if isinstance(obj, dict):
        return {k: to_dict(v) for k, v in obj.items()}
    return obj


def from_dict(cls, data: dict):
    """Build a dataclass from a dict, ignoring unknown keys and coercing
    numbers/strings/ints."""
    if not data:
        return cls()
    kwargs = {}
    for f in fields(cls):
        if f.name not in data:
            continue
        val = data[f.name]
        if val is None:
            continue
        # simple scalar coercion for config file robustness
        if isinstance(val, str):
            if val in ("auto", "None", "none"):
                # keep sentinels (e.g. 'auto' for motion.x_start_mm)
                kwargs[f.name] = val
                continue
            try:
                if f.type in ("int", int):
                    val = int(val)
                elif f.type in ("float", float):
                    val = float(val)
            except ValueError:
                pass  # keep the raw string
        kwargs[f.name] = val
    return cls(**kwargs)


CONFIG_SECTIONS = {
    "water": WaterParams,
    "jet": JetParams,
    "nozzle": NozzleParams,
    "motion": MotionParams,
    "sim": SimulationParams,
}


def parse_config_dict(cfg: dict) -> dict:
    """Turn a raw (YAML/JSON) split into the dataclass instances used by the
    simulator.  Returns {"water":..., "jet":..., "nozzle":..., ...}."""
    out = {}
    for key, cls in CONFIG_SECTIONS.items():
        out[key] = from_dict(cls, cfg.get(key) or {})
    return out
