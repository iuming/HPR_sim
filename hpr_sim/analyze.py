# -*- coding: utf-8 -*-
"""
Result post-processing: cleaning coverage and exposure metrics.
"""
from __future__ import annotations

import numpy as np

from .mesh import TriangleMesh
from .simulate import SimResults


def summarize(mesh: TriangleMesh, res: SimResults, roi=None, axis: str = "x") -> dict:
    """Compute the headline metrics of the simulation.

    All surface areas are returned in m^2 (converted from the mm^2 mesh data).

    `roi` : optional (xmin, xmax) axial window.  Only surface inside the window
    is taken as the *relevant* HPR domain; everything outside is ignored both in
    the denominator and numerator (it is usually adjacent hardware reached
    through open apertures, not RF surface).
    """
    ax = {"x": 0, "y": 1, "z": 2}[axis]
    area_mm2 = mesh.areas
    tri_hit = res.n_hits > 0
    tri_clean = res.n_effective > 0
    if roi is not None:
        xc = mesh.centroids[:, ax]
        inside = (xc >= roi[0]) & (xc <= roi[1])
        tri_hit &= inside
        tri_clean &= inside

    a_hit = float(area_mm2[tri_hit].sum())          # mm^2
    a_clean = float(area_mm2[tri_clean].sum())      # mm^2
    a_total = float(mesh.total_area)                # mm^2

    out = {
        "n_tri": int(mesh.n_tri),
        "surface_area_total_m2": a_total * 1e-6,
        "reachable_area_m2": a_hit * 1e-6,
        "cleaned_area_m2": a_clean * 1e-6,
        "effective_coverage_ratio": (a_clean / a_hit) if a_hit > 0 else 0.0,
        "cleaned_frac_of_mesh": (a_clean / a_total) if a_total > 0 else 0.0,
        "reachable_frac_of_mesh": (a_hit / a_total) if a_total > 0 else 0.0,
        "missed_area_m2": (a_total - a_hit) * 1e-6,
        "effective_only_area_m2": (a_hit - a_clean) * 1e-6,
    }

    # exposure on actually (effectively) cleaned area
    exp_clean = res.exposure_s[tri_clean]
    if len(exp_clean):
        out["mean_exposure_cleaned_s"] = float(exp_clean.mean())
        out["median_exposure_cleaned_s"] = float(np.median(exp_clean))
        out["max_exposure_cleaned_s"] = float(exp_clean.max())
    else:
        out.update(mean_exposure_cleaned_s=0.0, median_exposure_cleaned_s=0.0,
                   max_exposure_cleaned_s=0.0)

    hi = res.max_intensity[tri_hit]
    out["mean_max_intensity_hit"] = float(hi.mean()) if len(hi) else 0.0
    d = res.min_dist
    d = d[np.isfinite(d)]
    out["mean_min_dist_mm"] = float(d.mean()) if len(d) else np.nan
    out["mean_impact_bar_hit"] = float(res.max_impact_bar[tri_hit].mean()) if tri_hit.any() else 0.0
    out["wall_time_s"] = res.wall_time_s
    out["n_poses"] = res.n_poses
    out["jet_velocity_ms"] = res.v0_ms
    out["jet_flow_m3s"] = res.flow_m3s
    out["weber_exit"] = res.weber_exit
    for key in ("d_drop_um", "breakup_len_mm", "decay_coeff_1m",
                "re_exit", "cd_exit", "efold_mm", "drop_cd_mode"):
        if hasattr(res, key):
            out[key] = getattr(res, key)
    return out


def print_summary(mesh: TriangleMesh, res: SimResults, roi=None, axis: str = "x"):
    s = summarize(mesh, res, roi=roi, axis=axis)
    sep = "-" * 62
    print(sep)
    rtxt = f"  (ROI x: [{roi[0]:.1f}, {roi[1]:.1f}] mm)" if roi else ""
    print(f"HPR simulation result for {mesh.name}{rtxt}")
    print(sep)
    print(f"  mesh triangles                    : {s['n_tri']}")
    print(f"  total surface area                : {s['surface_area_total_m2']*1e6:10.0f} mm^2  ({s['surface_area_total_m2']:.4f} m^2)")
    print(f"  reachable (ever hit)              : {s['reachable_area_m2']*1e6:10.0f} mm^2  ({s['reachable_frac_of_mesh']*100:5.1f} % of mesh)")
    print(f"  effectively cleaned               : {s['cleaned_area_m2']*1e6:10.0f} mm^2")
    print(f"  -> effective coverage ratio       : {s['effective_coverage_ratio']*100:6.2f} %  (cleaned / reachable)")
    print(f"  missed (never reached) area       : {s['missed_area_m2']*1e6:10.0f} mm^2")
    print(f"  reached but not clean             : {s['effective_only_area_m2']*1e6:10.0f} mm^2")
    print(f"  mean exposure (cleaned, s)     : {s['mean_exposure_cleaned_s']:8.2f} s")
    print(f"  mean max intensity on hits     : {s['mean_max_intensity_hit']:.3f}")
    print(f"  mean impact pressure on hits   : {s['mean_impact_bar_hit']:8.2f} bar")
    print(f"  jet velocity / flow (exit)     : {s['jet_velocity_ms']:.1f} m/s, {s['jet_flow_m3s']*1e6:.2f} L/min per jet")
    print(f"  Weber number at exit           : {s['weber_exit']:.0f}")
    if "d_drop_um" in s:
        mode = s.get("drop_cd_mode", "?")
        if "re_exit" in s:
            print(f"  droplet model ({mode:16s}) : d_drop={s['d_drop_um']:.0f} um, "
                  f"L_breakup={s['breakup_len_mm']:.0f} mm, Re0={s['re_exit']:.0f}, "
                  f"Cd0={s['cd_exit']:.2f}, 1/e={s['efold_mm']:.0f} mm")
        else:
            print(f"  droplet model ({mode:16s}) : d_drop={s['d_drop_um']:.0f} um, "
                  f"L_breakup={s['breakup_len_mm']:.0f} mm, k={s['decay_coeff_1m']:.1f} 1/m")
    print(f"  wall time                      : {s['wall_time_s']:.1f} s")
    print(sep)
    return s
