# -*- coding: utf-8 -*-
"""
Mesh loading and preprocessing for HPR ray-casting simulation.

Reads STL files (binary or ASCII), computes face normals / areas / centroids,
and provides a lightweight AABB-assisted triangle lookup used by the ray
caster.
"""
from __future__ import annotations

import os
import struct

import numpy as np


class TriangleMesh:
    """Triangle mesh (vertices + triangle indices) with per-face geometry."""

    def __init__(self, vertices: np.ndarray, faces: np.ndarray, name: str = "mesh"):
        self.name = name
        self.vertices = np.asarray(vertices, dtype=np.float64)  # (N,3)
        self.faces = np.asarray(faces, dtype=np.int64)          # (M,3)
        self._compute_geometry()

    # ------------------------------------------------------------------
    # geometry builders
    # ------------------------------------------------------------------
    def _compute_geometry(self):
        v = self.vertices
        f = self.faces
        self.tri_v = v[f]                       # (M,3,3)
        self.e1 = self.tri_v[:, 1] - self.tri_v[:, 0]
        self.e2 = self.tri_v[:, 2] - self.tri_v[:, 0]
        self.normals = np.cross(self.e1, self.e2)               # (M,3), unnormalised
        nlen = np.linalg.norm(self.normals, axis=1)
        self.degenerate = np.asarray(nlen < 1e-12, dtype=bool)  # zero-area faces
        # guard against degenerate triangles before normalising
        nlen[nlen < 1e-30] = 1.0
        self.areas = 0.5 * nlen
        self.normals_unit = self.normals / nlen[:, None]
        self.centroids = self.tri_v.mean(axis=1)                # (M,3)
        # per-face bounding boxes
        self.bbox_min = self.tri_v.min(axis=1)                  # (M,3)
        self.bbox_max = self.tri_v.max(axis=1)                  # (M,3)
        self.n_tri = len(f)
        self.total_area = float(self.areas.sum())

    def submesh(self, mask) -> "TriangleMesh":
        return TriangleMesh(self.vertices, self.faces[mask], name=self.name + "_sub")

    # ------------------------------------------------------------------
    # STL I/O
    # ------------------------------------------------------------------
    @classmethod
    def from_stl(cls, path: str, name: str = None) -> "TriangleMesh":
        with open(path, "rb") as fh:
            data = fh.read()
        if name is None:
            name = os.path.splitext(os.path.basename(path))[0]
        if data[:5] == b"solid":
            return cls(*_parse_ascii_stl(data), name=name)
        return cls(*_parse_binary_stl(data), name=name)

    def write_stl(self, path: str):
        """Write a binary STL file (normals included)."""
        m = self.n_tri
        with open(path, "wb") as fh:
            fh.write(b"\0" * 80)
            fh.write(struct.pack("<I", m))
            flat = np.empty((m, 3, 3))
            flat[:] = self.tri_v
            # normals: use unit normals
            nrm = self.normals_unit
            for i in range(m):
                fh.write(struct.pack("<3f", *nrm[i]))
                for j in range(3):
                    fh.write(struct.pack("<3f", *flat[i, j]))
                fh.write(struct.pack("<H", 0))

    def summary(self) -> str:
        bbmin = self.vertices.min(axis=0)
        bbmax = self.vertices.max(axis=0)
        return (
            f"{self.name}: {self.n_tri} triangles, area={self.total_area/1e6:.3f} m^2, "
            f"bbox=[{bbmin[0]:.1f},{bbmin[1]:.1f},{bbmin[2]:.1f}] -> "
            f"[{bbmax[0]:.1f},{bbmax[1]:.1f},{bbmax[2]:.1f}] mm"
        )


def _parse_binary_stl(data: bytes):
    n = struct.unpack_from("<I", data, 80)[0]
    off = 84
    expected = 84 + n * 50
    if len(data) < expected:
        raise ValueError(f"binary STL truncated: {len(data)} < {expected}")
    count = 0
    tris = []
    for i in range(n):
        v0 = struct.unpack_from("<3f", data, off + 12)
        v1 = struct.unpack_from("<3f", data, off + 24)
        v2 = struct.unpack_from("<3f", data, off + 36)
        if v0 != v1 and v1 != v2 and v0 != v2:
            tris.append((v0, v1, v2))
            count += 1
        off += 50
    verts = np.asarray(tris, dtype=np.float64).reshape(-1, 3)   # (count*3, 3)
    faces = np.arange(count * 3, dtype=np.int64).reshape(count, 3)
    return verts, faces


def _parse_ascii_stl(data: bytes):
    text = data.decode("ascii", errors="ignore")
    coords = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("vertex"):
            coords.append([float(x) for x in s.split()[1:4]])
    coords = np.asarray(coords, dtype=np.float64)
    n = len(coords) // 3
    faces = np.arange(n * 3, dtype=np.int64).reshape(n, 3)
    return coords, faces


def merge_meshes(meshes, name="merged"):
    verts = []
    faces = []
    base = 0
    for m in meshes:
        faces.append(m.faces + base)
        verts.append(m.vertices)
        base += len(m.vertices)
    return TriangleMesh(np.vstack(verts), np.vstack(faces), name=name)
