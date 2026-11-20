# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""
Automatic geometric analysis of a cavity mesh.

For a general-purpose HPR tool the user should not have to hard-code the
beam axis or the scan region for every new cavity.  These helpers analyse the
STL mesh and propose:

  * the beam axis (longest or axis-of-symmetry diameter is searched), and
  * an HPR scan region [x_start, x_end] centred on the largest cell.

Heuristic (radial profile):
  1. r_min(x) = minimum wall radius per axial bin (lowest radius surface seen
     from the axis);
  2. the cell chain is delimited by the axial positions where r_min crosses
     the iris level (``iris_factor`` of the largest equator radius), so
     multi-cell cavities (e.g. 9-cell TESLA) scan from the first to the last
     iris;
  3. the scan range extends a small fraction of the chain span beyond each
     outermost iris (including the beam tubes right outside the cells);
  4. the default ray length follows from the largest cell radius and the
     nozzle tilt: ray ~ r_eq * sqrt(1 + tan(tilt)^2) + margin.

For cavities without a clean radial minimum the fallback is the central
20%-80% of the bounding box along the axis.
"""
from __future__ import annotations

import numpy as np

from .mesh import TriangleMesh

_AXIS = {"x": 0, "y": 1, "z": 2}


def radial_profile(mesh: TriangleMesh, axis: str = "x", n_bins: int = 400):
    """Minimum / median / maximum wall radius vs the axial coordinate."""
    ai = _AXIS[axis]
    c = mesh.centroids
    x = c[:, ai]
    r = np.sqrt((np.delete(c, ai, axis=1) ** 2).sum(axis=1))
    edges = np.linspace(x.min(), x.max(), n_bins + 1)
    idx = np.clip(np.digitize(x, edges) - 1, 0, n_bins - 1)
    rmin = np.full(n_bins, np.inf)
    rmed = np.zeros(n_bins)
    rmax = np.zeros(n_bins)
    cnt = np.zeros(n_bins, int)
    for i in range(n_bins):
        sel = idx == i
        cnt[i] = int(sel.sum())
        if cnt[i] == 0:
            rmin[i] = np.nan
        else:
            rs = r[sel]
            rmin[i] = rs.min()
            rmed[i] = np.median(rs)
            rmax[i] = rs.max()
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, rmin, rmed, rmax, cnt


def _smooth(a, k=5):
    if a.ndim == 0:
        return a
    kernel = np.ones(k) / k
    return np.convolve(np.nan_to_num(a, nan=np.inf), kernel, mode="same")


def detect_beam_axis(mesh: TriangleMesh) -> str:
    """Pick the axis with the longest extent (proxy for the beam axis)."""
    bbmin = mesh.vertices.min(axis=0)
    bbmax = mesh.vertices.max(axis=0)
    dims = bbmax - bbmin
    return "xyz"[int(np.argmax(dims))]


def propose_scan_range(mesh: TriangleMesh, axis: str = "x",
                       iris_factor: float = 0.7, n_bins: int = 400,
                       tilt_deg: float = 25.0, margin_mm: float = 30.0):
    """Propose [x_start, x_end], the equator position and a default ray
    length for the HPR wand scan.

    Multi-cell cavities (e.g. a 9-cell TESLA) are handled automatically: the
    scan region spans from the first to the last iris of the whole cell chain,
    not just the single largest cell.

    Returns a dict or ``None`` if no reliable profile could be found.
    """
    centers, rmin, rmed, rmax, cnt = radial_profile(mesh, axis, n_bins)
    ok = np.isfinite(rmin) & (cnt > 0)
    if ok.sum() < 10:
        return None
    # fill empty axial bins by interpolation so the profile is continuous
    x_valid = centers[ok]
    for arr in (rmin, rmed, rmax):
        arr[~ok] = np.interp(centers[~ok], x_valid, arr[ok])
    sm = _smooth(rmin)          # fully populated, no NaNs

    # cell equator = maximum of the (smoothed) minimum-radius profile
    xi = int(np.nanargmax(sm))
    req = sm[xi]
    if not np.isfinite(req) or req <= 0:
        return None
    thr = iris_factor * req

    # --- scan region: first->last iris of the cell chain -------------------
    # crossing points where the radial profile rises/falls through the iris
    # threshold delimit every cell bulge (tube -> iris -> cell -> iris -> ...)
    above = sm >= thr
    cross = np.flatnonzero(np.diff(above.astype(int)) != 0)
    if len(cross) >= 2:
        Li, Ri = int(cross[0]), int(cross[-1])
    else:
        # fallback: walk out from the largest equator (single-cell cavities
        # where the beam tube is wider than the iris threshold)
        Li = xi
        while Li > 0 and (sm[Li - 1] >= thr):
            Li -= 1
        Ri = xi
        while Ri < n_bins - 1 and (sm[Ri + 1] >= thr):
            Ri += 1

    span = centers[Ri] - centers[Li]
    if span <= 0 or not np.isfinite(span):
        return None

    x_iris_l, x_iris_r = centers[Li], centers[Ri]
    x_start = x_iris_l - 0.1 * span - margin_mm
    x_end = x_iris_r + 0.1 * span + margin_mm

    import math
    ray_len = req * math.sqrt(1.0 + math.tan(math.radians(abs(tilt_deg))) ** 2) + margin_mm
    return dict(x_start=float(x_start), x_end=float(x_end),
                equator_mm=float(centers[xi]), equator_radius_mm=float(req),
                iris_left_mm=float(x_iris_l), iris_right_mm=float(x_iris_r),
                ray_len_mm=float(ray_len))


def apply_auto(mesh: TriangleMesh, motion, sim, axis: str = "x"):
    """Fill ``motion.x_start_mm`` / ``x_end_mm`` and ``sim.ray_len_max_mm``
    with the automatic proposal when their value is ``'auto'`` (a string)."""
    auto_any = (motion.x_start_mm == "auto" or motion.x_end_mm == "auto"
                or sim.ray_len_max_mm == "auto")
    if not auto_any:
        return None
    prop = propose_scan_range(mesh, axis, tilt_deg=abs(getattr(motion, "tilt_fwd_deg", 25.0)))
    if prop is None:
        bbmin = mesh.vertices.min(axis=0)[_AXIS[axis]]
        bbmax = mesh.vertices.max(axis=0)[_AXIS[axis]]
        lo, hi = bbmin + 0.2 * (bbmax - bbmin), bbmin + 0.8 * (bbmax - bbmin)
        prop = dict(x_start=float(lo), x_end=float(hi),
                    ray_len_mm=0.5 * (hi - lo))
    if motion.x_start_mm == "auto":
        motion.x_start_mm = prop["x_start"]
    if motion.x_end_mm == "auto":
        motion.x_end_mm = prop["x_end"]
    if sim.ray_len_max_mm == "auto":
        sim.ray_len_max_mm = prop["ray_len_mm"]
    return prop
