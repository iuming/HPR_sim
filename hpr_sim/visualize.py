# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""
Visualisation of HPR simulation results.
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .mesh import TriangleMesh
from .simulate import SimResults


def _scatter(ax, mesh: TriangleMesh, values, vmin, vmax, cmap, alpha=1.0,
             s=1.2, nr=40000):
    """Scatter triangle centroids coloured by a per-triangle value."""
    vals = np.asarray(values, dtype=float)
    stride = max(1, int(len(vals) // nr))
    v = vals[::stride]
    c = mesh.centroids[::stride]
    m = np.isfinite(v)
    sc = ax.scatter(c[m, 0], c[m, 1], c[m, 2], c=v[m], cmap=cmap,
                    vmin=vmin, vmax=vmax, s=s, alpha=alpha, linewidths=0)
    return sc


def _scatter_rgb(ax, mesh: TriangleMesh, mask, color, alpha=1.0, s=2.0,
                 nr=100000):
    """Scatter the (masked) centroids with a single vivid RGB colour."""
    idx = np.nonzero(mask)[0]
    if len(idx) == 0:
        return
    stride = max(1, int(len(idx) // nr))
    idx = idx[::stride]
    c = mesh.centroids[idx]
    ax.scatter(c[:, 0], c[:, 1], c[:, 2], color=color, s=s, alpha=alpha,
               linewidths=0, depthshade=False)


def _set_true_aspect(ax, pts: np.ndarray = None):
    """Scale the 3-D axes consistently so 1 mm has the same screen length on
    every axis (box aspect follows the displayed data extents).

    ``pts`` = the points actually visible (e.g. a half-cutaway subset); the
    limits and the box aspect are then made proportional to that subset, so
    there is no wasted empty space and the proportions match the picture.
    """
    pts = np.asarray(pts) if pts is not None else None
    if pts is None or len(pts) == 0:
        return
    bbmin = pts.min(axis=0)
    bbmax = pts.max(axis=0)
    span = np.maximum(bbmax - bbmin, 1e-6)
    try:
        ax.set_box_aspect((float(span[0]), float(span[1]), float(span[2])))
    except Exception:
        pass
    mid = 0.5 * (bbmin + bbmax)
    half = 0.5 * span + 0.02 * span.max()   # small symmetric margin
    for i, name in enumerate(("set_xlim", "set_ylim", "set_zlim")):
        try:
            getattr(ax, name)(mid[i] - half[i], mid[i] + half[i])
        except Exception:
            pass


def coverage_cloud(out_path: str, mesh: TriangleMesh, res: SimResults,
                   title: str = "HPR cleaning coverage", show_missed=True,
                   roi_mm=None, elev: float = 10.0, azim: float = -60.0,
                   context_nr: int = 6000, wireframe: bool = True,
                   cut: str = None):
    """3D coverage cloud map with vivid layered colours.

    Layout:

      * **effectively-cleaned faces** -> bright green, opaque, large markers;
      * **reached but not cleaned**   -> amber;
      * **background**                -> the *whole* cavity only as a very
        sparse, faint grey outline (≈``context_nr`` points, alpha≈0.06) so it
        does not occlude the inner surface, but keeps the cavity shape legible.

    ``roi_mm`` (x_start, x_end) restricts the *reached/cleaned* layers to the
    wand-scan area (the background outline always shows the full cavity).

    ``cut``: if ``'x'``/``'y'``/``'z'``, display only the **minus half** of
    the cavity cut by the plane through its centre along that axis and look
    at it from above -- a classic half-cavity cutaway that shows the inner RF
    surface unobstructed.
    """
    reached = res.n_hits > 0
    # reached faces that lie inside the scan region
    in_roi = np.ones(mesh.n_tri, dtype=bool)
    if roi_mm is not None:
        xc = mesh.centroids[:, 0]
        in_roi = (xc >= roi_mm[0]) & (xc <= roi_mm[1])
    cleaned = (res.n_effective > 0) & in_roi
    not_clean = reached & (~(res.n_effective > 0)) & in_roi

    keep = np.ones(mesh.n_tri, dtype=bool)
    if cut in ("x", "y", "z"):
        ai = {"x": 0, "y": 1, "z": 2}[cut]
        keep = mesh.centroids[:, ai] <= 0.0

    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111, projection="3d")

    # 1) faint full-cavity outline for orientation (never occludes)
    if wireframe:
        _scatter_rgb(ax, mesh, keep, "#8894a6", alpha=0.08, s=1.0,
                     nr=context_nr)
        if roi_mm is not None:
            for xb in (roi_mm[0], roi_mm[1]):
                ax.plot([xb, xb], [mesh.centroids[:, 1].min(),
                                   mesh.centroids[:, 1].max()],
                        [0, 0], color="#666666", lw=1, alpha=0.5)
    # 2) reached but not effectively cleaned
    _scatter_rgb(ax, mesh, not_clean & keep, "#ff9900", alpha=0.7, s=2.4,
                 nr=80000)
    # 3) effectively cleaned - vivid green, drawn last (on top)
    _scatter_rgb(ax, mesh, cleaned & keep, "#00e676", alpha=1.0, s=3.8,
                 nr=120000)

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    if cut in ("x", "y", "z"):
        # look straight down the cut axis
        if cut == "x":
            elev, azim = 0.0, 90.0
        elif cut == "y":
            elev, azim = 0.0, 0.0
        else:
            elev, azim = 89.0, -60.0
    ax.view_init(elev=elev, azim=azim)
    # make the axes proportional to what is actually shown (cutaway subset
    # for the half views)
    _set_true_aspect(ax, mesh.centroids[keep] if cut in ("x", "y", "z")
                     else mesh.vertices)

    # legend
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", ls="", color="#00e676", ms=9,
               label="effectively cleaned"),
        Line2D([0], [0], marker="o", ls="", color="#ff9900", ms=8,
               label="reached, not cleaned"),
        Line2D([0], [0], marker="o", ls="", color="#8894a6", ms=7,
               alpha=0.4, label="cavity surface (outline)"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=9, framealpha=0.85)

    fig.tight_layout(pad=0.4)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def intensity_cloud(out_path: str, mesh: TriangleMesh, res: SimResults,
                    title: str = "HPR max cleaning intensity (cos phi)"):
    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111, projection="3d")
    hit = res.n_hits > 0
    views = _scatter(ax, mesh, np.where(hit, res.max_intensity, -1), -1, 1,
                     "hot", alpha=0.95, s=2.0, nr=120000)
    fig.colorbar(views, ax=ax, shrink=0.6, label="max intensity")
    ax.set_title(title)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    ax.view_init(elev=12.0, azim=-58.0)
    _set_true_aspect(ax, mesh.vertices)
    fig.tight_layout(pad=0.4)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def unwrapped_map(out_path: str, mesh: TriangleMesh, res: SimResults,
                  axis: str = "x", cell_mm: float = 3.0, n_phi: int = 180,
                  title: str = "Effective coverage, unrolled (axial x azimuth)"):
    """Unrolled (axial-angle) coverage map for roughly axisymmetric cavities."""
    ai = {"x": 0, "y": 1, "z": 2}[axis]
    c = mesh.centroids
    xa = c[:, ai]
    u1 = (ai + 1) % 3
    u2 = (ai + 2) % 3
    phi = np.arctan2(c[:, u2], c[:, u1])          # [-pi, pi]
    xmin, xmax = xa.min(), xa.max()
    nx = max(2, int(np.ceil((xmax - xmin) / cell_mm)))
    bins_x = np.linspace(xmin, xmax, nx + 1)
    bins_p = np.linspace(-np.pi, np.pi, n_phi + 1)

    # effective coverage fraction per bin, weighted by area
    eff = (res.n_effective > 0).astype(float)
    hit = (res.n_hits > 0).astype(float)
    area = mesh.areas
    X = (np.digitize(xa, bins_x) - 1)
    P = (np.digitize(phi, bins_p) - 1)
    valid = (X >= 0) & (X < nx) & (P >= 0) & (P < n_phi)
    weight = np.bincount(X[valid] * n_phi + P[valid], weights=area[valid],
                         minlength=nx * n_phi).reshape(nx, n_phi)
    weff = np.bincount(X[valid] * n_phi + P[valid], weights=area[valid] * eff[valid],
                       minlength=nx * n_phi).reshape(nx, n_phi)
    whit = np.bincount(X[valid] * n_phi + P[valid], weights=area[valid] * hit[valid],
                       minlength=nx * n_phi).reshape(nx, n_phi)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_eff = np.where(whit > 0, weff / np.maximum(whit, 1e-12), np.nan)
        frac_reach = np.where(weight > 0, whit / np.maximum(weight, 1e-12), np.nan)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax0 = axes[0]
    im0 = ax0.imshow(frac_eff.T, origin="lower", aspect="auto",
                     extent=[xmin, xmax, -np.pi, np.pi],
                     cmap="RdYlGn", vmin=0, vmax=1)
    ax0.set_ylabel("azimuth [rad]")
    ax0.set_title(title + "\n(effective / reachable per bin)")
    fig.colorbar(im0, ax=ax0)
    ax1 = axes[1]
    im1 = ax1.imshow(frac_reach.T, origin="lower", aspect="auto",
                     extent=[xmin, xmax, -np.pi, np.pi],
                     cmap="Blues", vmin=0, vmax=1)
    ax1.set_xlabel(f"{axis} [mm]")
    ax1.set_ylabel("azimuth [rad]")
    ax1.set_title("reachable fraction of all surface per bin")
    fig.colorbar(im1, ax=ax1)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def coverage_topview(out_path: str, mesh: TriangleMesh, res: SimResults,
                     cut: str = "z", half: int = -1, roi_mm=None,
                     title: str = "HPR coverage (top view, 1:1)"):
    """Flat 2-D top view of a half-cavity cutaway with **strictly equal** x/y
    spacing (``aspect='equal'``), so the axes are exactly to scale.

    ``cut`` : axis through whose mid-plane the cavity is cut
    ``half`` : -1 keeps centroid coord <= 0, +1 keeps >= 0 on the cut axis.
    The two remaining axes are drawn as a 2-D scatter map (green = effectively
    cleaned, orange = reached but not cleaned, grey = untouched context).
    """
    ai = {"x": 0, "y": 1, "z": 2}[cut]
    other = [i for i in range(3) if i != ai]
    names = ["x [mm]", "y [mm]", "z [mm]"]
    keep = mesh.centroids[:, ai]
    keep = keep <= 0.0 if half < 0 else keep >= 0.0

    reached = res.n_hits > 0
    in_roi = np.ones(mesh.n_tri, dtype=bool)
    if roi_mm is not None:
        xc = mesh.centroids[:, 0]
        in_roi = (xc >= roi_mm[0]) & (xc <= roi_mm[1])
    cleaned = (res.n_effective > 0) & in_roi
    not_clean = reached & (~(res.n_effective > 0)) & in_roi

    fig, ax = plt.subplots(figsize=(9.2, 7.0))

    def plot2(mask, color, alpha, s, nr):
        idx = np.nonzero(mask & keep)[0]
        if len(idx) == 0:
            return
        stride = max(1, int(len(idx) // nr))
        c = mesh.centroids[idx[::stride]]
        ax.scatter(c[:, other[0]], c[:, other[1]], color=color, s=s,
                   alpha=alpha, linewidths=0)

    plot2(np.ones(mesh.n_tri, dtype=bool), "#8894a6", 0.18, 0.8, 60000)
    plot2(not_clean, "#ff9900", 0.75, 3.0, 80000)
    plot2(cleaned, "#00e676", 1.0, 5.0, 120000)

    ax.set_xlabel(names[other[0]])
    ax.set_ylabel(names[other[1]])
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")   # strict 1:1 on both axes
    ax.grid(alpha=0.25)

    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0], [0], marker="o", ls="", color="#00e676", ms=8,
               label="effectively cleaned"),
        Line2D([0], [0], marker="o", ls="", color="#ff9900", ms=7,
               label="reached, not cleaned"),
        Line2D([0], [0], marker="o", ls="", color="#8894a6", ms=6, alpha=0.5,
               label="cut surface (context)"),
    ], loc="upper right", fontsize=9, framealpha=0.85)

    fig.tight_layout(pad=0.5)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def histograms(out_path: str, mesh: TriangleMesh, res: SimResults):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    hit = res.n_hits > 0
    ax = axes[0]
    ax.hist(res.exposure_s[res.n_effective > 0], bins=40,
            color="steelblue", edgecolor="white")
    ax.set_xlabel("exposure time [s]")
    ax.set_ylabel("faces")
    ax.set_title("Exposure on cleaned area")
    ax = axes[1]
    ax.hist(res.max_intensity[hit], bins=40, color="coral", edgecolor="white")
    ax.set_xlabel("max intensity (cos phi)")
    ax.set_title("Max intensity distribution (hit faces)")
    ax = axes[2]
    ax.hist(res.max_impact_bar[hit], bins=40, color="seagreen", edgecolor="white")
    ax.set_xlabel("max impact pressure [bar]")
    ax.set_title("Impact pressure on hits")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path
