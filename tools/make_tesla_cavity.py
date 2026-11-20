#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parametric builder of a TESLA 1.3 GHz 9-cell cavity (inner niobium surface).

The half-cell shape follows the TESLA cavity definition published in
B. Aune et al., "The Superconducting TESLA Cavities", Phys. Rev. ST Accel.
Beams 3, 092001 (2000) - Table "Half-cell shape parameters" (all in mm):

    midcup :  R_equat = 103.3,  R_iris = 35.0,  R_arc = 42.0,
              iris ellipse half-axes a = 12, b = 19,   l = 57.7
    endcup :  R_equat = 103.3,  R_iris = 39.0,  R_arc = 40.3,
              iris ellipse half-axes a = 10, b = 13.5, l = 56.0

Construction of one half cell (z: 0 = equator plane -> l = iris plane):

    * iris ellipse  : centre (l - a, R_iris); the wall leaves the iris with
      a vertical tangent (right-most point of the ellipse) and sweeps to the
      bottom point (r = R_iris + b, horizontal tangent);
    * circular arc  : centre (0, R_equat - R_arc), radius R_arc, starts at the
      equator with horizontal tangent;
    * straight band : the common tangent between the two arcs.

The full 9-cell cavity = 2 end cells + 7 mid cells + beam tubes (r = 39 mm).
The inner wall is a single-valued profile r(z) revolved around the axis and
closed with end caps -> one water-tight surface, ready for ray casting.

Usage
-----
    python tools/make_tesla_cavity.py [--stl data/tesla_9cell_inner.stl]
                                      [--dz 1.0] [--nphi 48] [--plot profile.png]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from hpr_sim.mesh import TriangleMesh  # noqa: E402

# ---------------------------------------------------------------------------
# TESLA half-cell parameters (Aune et al., PRST-AB 3, 092001 (2000))
# ---------------------------------------------------------------------------
MIDCUP = dict(l=57.7, R_iris=35.0, R_equat=103.3, R_arc=42.0, a=12.0, b=19.0)
ENDCUP = dict(l=56.0, R_iris=39.0, R_equat=103.3, R_arc=40.3, a=10.0, b=13.5)
BEAM_TUBE_R = 39.0
BEAM_TUBE_L = 100.0
N_MIDCELLS = 7


def _common_tangent(p: dict, kmax: float = 2.0):
    """Slope/intercept of the straight section tangent to both the equator
    circular arc and the iris ellipse, together with the two tangent points
    (z_t1, r_t1 on the circle; z_t2, r_t2 on the ellipse)."""
    zc, rc = 0.0, p["R_equat"] - p["R_arc"]
    Ra = p["R_arc"]
    ce_z, ce_r = p["l"] - p["a"], p["R_iris"]
    a, b = p["a"], p["b"]

    # line y = k z + q ; circle tangency gives q as a function of k
    def ellipse_disc(k, q):
        A = 1 / a ** 2 + k * k / b ** 2
        B = 2 * (k * (q - ce_r) / b ** 2 - ce_z / a ** 2)
        C = ce_z ** 2 / a ** 2 + (q - ce_r) ** 2 / b ** 2 - 1
        return B * B - 4 * A * C

    def f(k):
        q = rc - k * zc - Ra * np.sqrt(1 + k * k)
        return ellipse_disc(k, q)

    # bisection on the slope k in (0, kmax]
    lo, hi = 1e-4, kmax
    f_lo, f_hi = f(lo), f(hi)
    if f_lo * f_hi > 0:                 # no sign change -> fallback slope
        k = 0.5
    else:
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            f_mid = f(mid)
            if f_lo * f_mid <= 0:
                hi, f_hi = mid, f_mid
            else:
                lo, f_lo = mid, f_mid
        k = 0.5 * (lo + hi)
    q = rc - k * zc - Ra * np.sqrt(1 + k * k)

    # tangent point on the circle (lower-right branch: tangent slope == k)
    z_t1 = zc + Ra * (k / np.sqrt(1 + k * k))
    r_t1 = rc + Ra * (1.0 / np.sqrt(1 + k * k))
    # tangent point on the ellipse (larger-z branch)
    A = 1 / a ** 2 + k * k / b ** 2
    B = 2 * (k * (q - ce_r) / b ** 2 - ce_z / a ** 2)
    C = ce_z ** 2 / a ** 2 + (q - ce_r) ** 2 / b ** 2 - 1
    z_t2 = (-B - np.sqrt(max(B * B - 4 * A * C, 0.0))) / (2 * A)
    r_t2 = k * z_t2 + q
    return k, (float(z_t1), float(r_t1)), (float(z_t2), float(r_t2))


def half_cell_r(p: dict, n: int = 400):
    """r(z) samples of one half cell; z in [0, l] from equator to iris."""
    _, t1, t2 = _common_tangent(p)
    z_t1, r_t1 = t1
    z_t2, r_t2 = t2
    zc, rc = 0.0, p["R_equat"] - p["R_arc"]
    Ra, ce_z, ce_r = p["R_arc"], p["l"] - p["a"], p["R_iris"]
    a, b = p["a"], p["b"]

    z = np.linspace(0.0, p["l"], n)
    r = np.empty_like(z)
    for i, zi in enumerate(z):
        if zi <= z_t1:
            r[i] = rc + np.sqrt(max(Ra * Ra - zi * zi, 0.0))
        elif zi <= z_t2:
            k = (r_t2 - r_t1) / max(z_t2 - z_t1, 1e-9)
            r[i] = r_t1 + k * (zi - z_t1)
        else:
            r[i] = ce_r + b * np.sqrt(max(1 - ((zi - ce_z) / a) ** 2, 0.0))
    # enforce monotonic decrease (numerical safety)
    r = np.maximum.accumulate(r[::-1])[::-1]
    return z, r


def full_profile(dz: float = 1.0):
    """Return (z, r) of the whole 9-cell inner surface, z=0 at the centre."""
    U = []

    def add_half(pdict, rising: bool):
        zz, rr = half_cell_r(pdict, n=max(2, int(np.ceil(pdict["l"] / dz))))
        if not rising:                      # equator -> iris (descending)
            zz = _z[0] + zz
        else:                               # iris -> equator (rising)
            zz = _z[0] + pdict["l"] - zz
        U.extend(zip(zz, rr))
        _z[0] += pdict["l"]

    # most-negative point: beam tube + cup half + 8 mid halves (down to the
    # equator plane of the central cell 5, which becomes z = 0)
    _z = [-(BEAM_TUBE_L + ENDCUP["l"] + (N_MIDCELLS + 1) * MIDCUP["l"])]
    z_start_ = _z[0]

    # negative-side beam tube (iris radius = 39)
    U.append((_z[0], BEAM_TUBE_R))
    n_tube = int(np.ceil(BEAM_TUBE_L / dz))
    for i in range(1, n_tube + 1):
        U.append((_z[0] + i * BEAM_TUBE_L / n_tube, BEAM_TUBE_R))
    _z[0] += BEAM_TUBE_L

    # the 18 half cells (2 end cups + 16 mid cups), see module docstring
    add_half(ENDCUP, rising=True)          # h1  : cup, iris0 -> eq1
    add_half(MIDCUP, rising=False)         # h2  : eq1 -> iris1
    for _ in range(N_MIDCELLS):            # h3..h16 : cells 2..8
        add_half(MIDCUP, rising=True)
        add_half(MIDCUP, rising=False)
    add_half(MIDCUP, rising=True)          # h17 : iris8 -> eq9
    add_half(ENDCUP, rising=False)         # h18 : eq9 -> iris9 (beam tube)

    # positive-side beam tube
    for i in range(1, n_tube + 1):
        U.append((_z[0] + i * BEAM_TUBE_L / n_tube, BEAM_TUBE_R))
    _z[0] += BEAM_TUBE_L

    U = np.array(sorted(U, key=lambda x: x[0]))
    keep = np.abs(np.diff(U[:, 0])) > 1e-6
    U = np.concatenate([U[:1], U[1:][keep]])
    return U[:, 0], U[:, 1], z_start_


def revolve_stl(z, r, nphi: int = 48):
    """Revolve a single-valued profile r(z) into a closed (water-tight)
    triangle mesh and return it."""
    z = np.asarray(z, float)
    r = np.asarray(r, float)
    n = nphi
    th = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    verts = np.empty((len(z) * n, 3))
    for i in range(len(z)):
        slice_ = verts[i * n:(i + 1) * n]
        slice_[:, 0] = z[i]
        slice_[:, 1] = r[i] * np.cos(th)
        slice_[:, 2] = r[i] * np.sin(th)

    faces = []
    for i in range(len(z) - 1):
        for j in range(n):
            j2 = (j + 1) % n
            a0 = i * n + j
            a1 = i * n + j2
            a2 = (i + 1) * n + j2
            a3 = (i + 1) * n + j
            faces.append((a0, a1, a2))
            faces.append((a0, a2, a3))
    # end caps (fan to the axis)
    cap0 = len(z) * n
    for j in range(n):
        faces.append((j, (j + 1) % n, cap0))
    cap1 = cap0 + 1
    last = (len(z) - 1) * n
    for j in range(n):
        faces.append((last + j, cap1, last + (j + 1) % n))

    verts = np.vstack([verts, [[z[0], 0, 0]], [[z[-1], 0, 0]]])
    return TriangleMesh(verts, np.asarray(faces, dtype=np.int64),
                        name="tesla_9cell_inner")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stl", default=os.path.join("data", "tesla_9cell_inner.stl"))
    ap.add_argument("--dz", type=float, default=1.5)
    ap.add_argument("--nphi", type=int, default=48)
    ap.add_argument("--plot", default="")
    args = ap.parse_args()

    z, r, z_start = full_profile(args.dz)
    mesh = revolve_stl(z, r, args.nphi)
    os.makedirs(os.path.dirname(os.path.abspath(args.stl)) or ".", exist_ok=True)
    mesh.write_stl(args.stl)
    print(f"written {args.stl}: {mesh.n_tri} triangles, "
          f"z=[{z.min():.1f},{z.max():.1f}] mm, r_max={r.max():.1f} mm")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.plot(z, r, lw=1.0, color="tab:blue")
        ax.set_aspect("equal")
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("r [mm]")
        ax.set_title("TESLA 9-cell inner profile (1.3 GHz)")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"profile plot -> {args.plot}")


if __name__ == "__main__":
    main()
