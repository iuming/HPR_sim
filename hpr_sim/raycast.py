# -*- coding: utf-8 -*-
"""
Vectorised Möller-Trumbore ray-triangle intersection and conservative
candidate pre-filtering.

The nozzle is modelled as a bundle of rays sharing a single origin (as in
CN121302684A where the spray head is a set of direction vectors from a common
starting point O), so for every emission step we have one origin and a small
set of ray directions.

Candidate pre-filter: only triangles whose centroid lies inside the cone+sphere
volume spanned by the ray set are tested.  The cone is the union cone of all ray
directions; the sphere radius is the maximum ray length (plus a per-triangle
bbox margin) so that no hit can be missed.
"""
from __future__ import annotations

import numpy as np

from .mesh import TriangleMesh

_EPS = 1e-9


def candidate_mask(
    mesh: TriangleMesh,
    origin: np.ndarray,
    cone_axis: np.ndarray,
    cone_half_angle: float,
    ray_len: float,
) -> np.ndarray:
    """Boolean mask over mesh triangles that may intersect the ray bundle.

    Parameters
    ----------
    mesh : TriangleMesh
    origin : (3,) array, common ray origin
    cone_axis : (3,) unit vector, axis of the union cone
    cone_half_angle : float (radians)
    ray_len : float, maximum length along any ray to consider (mm)
    """
    c = mesh.centroids
    d = c - origin
    dist = np.linalg.norm(d, axis=1)
    # conservative sphere test (per-triangle bbox margin)
    tri_r = 0.5 * np.linalg.norm(mesh.bbox_max - mesh.bbox_min, axis=1)
    rmax = ray_len + tri_r
    ok_sphere = dist <= rmax

    # conservative cone test
    avec = d / np.maximum(dist, 1e-9)[:, None]
    cosang = np.abs(avec @ cone_axis)
    # angular margin from triangle size
    sin_marg = np.minimum(tri_r / np.maximum(dist, 1e-6), 1.0)
    marg = np.arcsin(sin_marg)
    cos_marg = np.cos(cone_half_angle + marg)
    ok_cone = cosang >= cos_marg
    return ok_sphere & ok_cone


def ray_triangle_intersect(
    origin: np.ndarray,
    dirs: np.ndarray,
    mesh: TriangleMesh,
    tri_idx: np.ndarray,
    tmin: float = 1e-6,
    tmax: float = np.inf,
) -> dict:
    """Intersect batched rays (all sharing `origin`) against a sub-set of
    triangles.

    Parameters
    ----------
    origin : (3,) float
    dirs : (R,3) unit ray directions
    mesh : TriangleMesh
    tri_idx : (C,) int64 indices of candidate triangles
    tmin, tmax : float, accepted ray-parameter window

    Returns
    -------
    dict with keys (all arrays with one entry per ray, -1 / NaN when missed):
        hit      (R,) bool
        tri      (R,) int   - hit triangle index into mesh.faces
        t        (R,) float - travel distance in mm
        u,v      (R,) float - barycentric coordinates
        cos_phi  (R,) float - |dot(ray, face normal)| (cleaning intensity)
    """
    C = len(tri_idx)
    if C == 0:
        R = len(dirs)
        return dict(hit=np.zeros(R, bool), tri=np.full(R, -1, np.int64),
                    t=np.full(R, np.nan), u=np.full(R, np.nan),
                    v=np.full(R, np.nan), cos_phi=np.full(R, np.nan))

    v0 = mesh.tri_v[tri_idx, 0]            # (C,3)
    e1 = mesh.e1[tri_idx]
    e2 = mesh.e2[tri_idx]
    nrm = mesh.normals_unit[tri_idx]

    R = len(dirs)
    hit_tri = np.full(R, -1, np.int64)
    hit_t = np.full(R, np.inf)
    hit_u = np.zeros(R)
    hit_v = np.zeros(R)
    hit_cos = np.zeros(R)
    hit_mask = np.zeros(R, bool)

    # Precompute det for all candidates once (independent of ray)
    # We loop over rays and vectorise over candidates.
    for r in range(R):
        d = dirs[r]
        p = np.cross(d[None, :], e2)            # (C,3)
        det = np.einsum("ij,ij->i", e1, p)
        ok_det = np.abs(det) > _EPS
        if not ok_det.any():
            continue
        inv = 1.0 / det
        tv = origin - v0
        uu = np.einsum("ij,ij->i", tv, p) * inv
        ok_u = (uu >= 0.0) & (uu <= 1.0)
        q = np.cross(tv, e1)
        vv = np.einsum("ij,ij->i", d[None, :], q) * inv
        ok_v = (vv >= 0.0) & (uu + vv <= 1.0)
        tt = np.einsum("ij,ij->i", e2, q) * inv
        ok_t = (tt > tmin) & (tt < tmax)
        cand = ok_det & ok_u & ok_v & ok_t
        if not cand.any():
            continue
        idx = np.nonzero(cand)[0]
        k = idx[np.argmin(tt[idx])]
        if tt[k] < hit_t[r]:
            hit_t[r] = tt[k]
            hit_tri[r] = tri_idx[k]
            hit_u[r] = uu[k]
            hit_v[r] = vv[k]
            hit_cos[r] = abs(np.dot(d, nrm[k]))
            hit_mask[r] = True

    return dict(hit=hit_mask, tri=hit_tri, t=hit_t, u=hit_u, v=hit_v,
                cos_phi=hit_cos)


def intersect_bundle_fast(
    origin: np.ndarray,
    dirs: np.ndarray,
    mesh: TriangleMesh,
    cone_axis: np.ndarray,
    cone_half: float,
    ray_len: float,
    tmin: float = 1e-6,
    sub_cone_deg: float = 10.0,
) -> dict:
    """Fully-vectorised intersection of a whole ray bundle against the mesh.

    Candidate filtering is done in two stages:

      1. sphere stage: triangles whose centroid lies within `ray_len` (+ a
         per-triangle bbox margin) of the common origin -- the union cone is
         usually ~180 deg for a ring of tilted holes, so a sphere filter is
         just as tight and much cheaper;
      2. ray stage: per ray, keep only triangles whose centroid lies within a
         small cone (``sub_cone_deg`` + triangle angular margin) around that
         ray's direction.

    The surviving (ray, triangle) pairs are then intersected in ONE batch
    (Möller-Trumbore), and the closest hit per ray is kept.

    Returns the same dict as :func:`ray_triangle_intersect`.
    """
    R = len(dirs)
    hit_tri = np.full(R, -1, np.int64)
    hit_t = np.full(R, np.inf)
    hit_u = np.zeros(R)
    hit_v = np.zeros(R)
    hit_cos = np.zeros(R)
    hit_mask = np.zeros(R, bool)

    # ---- stage 1: sphere around the common origin -------------------------
    dc_all = mesh.centroids - origin
    dist_all = np.linalg.norm(dc_all, axis=1)
    tri_r = 0.5 * np.linalg.norm(mesh.bbox_max - mesh.bbox_min, axis=1)
    cand = dist_all <= (ray_len + tri_r)
    tri_idx = np.nonzero(cand)[0]
    if len(tri_idx) == 0:
        return dict(hit=hit_mask, tri=hit_tri, t=hit_t, u=hit_u, v=hit_v,
                    cos_phi=hit_cos)

    # ---- stage 2: per-ray narrow cone --------------------------------------
    dc = dc_all[tri_idx]
    dist = dist_all[tri_idx]
    avec = dc / np.maximum(dist, 1e-9)[:, None]
    r_tri = tri_r[tri_idx]
    marg = np.arcsin(np.minimum(r_tri / np.maximum(dist, 1e-6), 1.0))
    cos_lim = np.cos(np.deg2rad(sub_cone_deg) + marg)          # (C1,)

    cosang = avec @ dirs.T                                     # (C1, R)
    ok = cosang >= cos_lim[:, None]
    ci, ri = np.nonzero(ok)                                    # surviving pairs
    # note: ok[row, col] -> ci are candidate rows, ri are ray columns
    if len(ri) == 0:
        return dict(hit=hit_mask, tri=hit_tri, t=hit_t, u=hit_u, v=hit_v,
                    cos_phi=hit_cos)

    gid = tri_idx[ci]                                          # global face ids

    # ---- batched Möller-Trumbore -------------------------------------------
    with np.errstate(divide="ignore", invalid="ignore"):
        v0 = mesh.tri_v[gid, 0]
        e1 = mesh.e1[gid]
        e2 = mesh.e2[gid]
        d = dirs[ri]
        p = np.cross(d, e2)
        det = np.einsum("ij,ij->i", e1, p)
        good_det = np.abs(det) > _EPS
        if not good_det.any():
            return dict(hit=hit_mask, tri=hit_tri, t=hit_t, u=hit_u, v=hit_v,
                        cos_phi=hit_cos)
        inv = 1.0 / det
        tv = origin - v0
        uu = np.einsum("ij,ij->i", tv, p) * inv
        ok_u = (uu >= 0.0) & (uu <= 1.0)
        q = np.cross(tv, e1)
        vv = np.einsum("ij,ij->i", d, q) * inv
        ok_v = (vv >= 0.0) & (uu + vv <= 1.0)
        tt = np.einsum("ij,ij->i", e2, q) * inv
        ok_t = (tt > tmin) & (tt < ray_len)
        valid = good_det & ok_u & ok_v & ok_t
    if not valid.any():
        return dict(hit=hit_mask, tri=hit_tri, t=hit_t, u=hit_u, v=hit_v,
                    cos_phi=hit_cos)

    ri_v = ri[valid]
    gid_v = gid[valid]
    t_v = tt[valid]
    u_v = uu[valid]
    v_v = vv[valid]
    cos_v = np.abs(np.einsum("ij,ij->i", d[valid], mesh.normals_unit[gid_v]))

    # closest hit per ray
    order = np.lexsort((t_v, ri_v))
    r_s, g_s, t_s, u_s, v_s, c_s = (ri_v[order], gid_v[order], t_v[order],
                                    u_v[order], v_v[order], cos_v[order])
    ur, first = np.unique(r_s, return_index=True)
    hit_mask[ur] = True
    hit_tri[ur] = g_s[first]
    hit_t[ur] = t_s[first]
    hit_u[ur] = u_s[first]
    hit_v[ur] = v_s[first]
    hit_cos[ur] = c_s[first]

    return dict(hit=hit_mask, tri=hit_tri, t=hit_t, u=hit_u, v=hit_v,
                cos_phi=hit_cos)


def cone_of_dirs(dirs: np.ndarray):
    """Union cone (axis, half-angle) covering a set of ray directions."""
    d = dirs / np.linalg.norm(dirs, axis=1, keepdims=True)
    # find direction maximally centred: direction of the sum; for symmetric
    # ray sets the sum can vanish, then fall back to the first ray
    axis = d.sum(axis=0)
    nrm = np.linalg.norm(axis)
    if nrm < 1e-9:
        axis = d[0].copy()
    else:
        axis = axis / nrm
    cosang = d @ axis
    cosang = np.clip(cosang, -1, 1)
    half = np.arccos(cosang.min())
    return axis, float(half)
