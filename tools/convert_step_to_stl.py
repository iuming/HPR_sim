#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert a STEP (or other CAD) assembly/solid into an STL triangle mesh using
FreeCAD in command-line (headless) mode.

FreeCAD's freecadcmd does not forward arbitrary trailing CLI arguments to the
script (it tries to open them as documents), so all parameters are passed
through environment variables:

    STEP2STL_SRC  : input  STEP file (absolute path)
    STEP2STL_DST  : output STL file (absolute path)
    STEP2STL_LIN_DEF : linear deflection in mm (optional, default 0.5)
    STEP2STL_ANG_DEF : angular deflection in degrees (optional, default 10)

Example (PowerShell):
    $env:STEP2STL_SRC  = "D:\\...\\model.STEP"
    $env:STEP2STL_DST  = "D:\\...\\model.stl"
    freecadcmd.exe convert_step_to_stl.py
"""
import os
import sys
import time

LOG = os.environ.get("STEP2STL_LOG")

# Early marker so we can see the script executing even before main()
with open(os.environ.get("STEP2STL_MARKER", os.devnull), "w") as _mf:
    _mf.write("started\n")


def log(msg):
    print(msg, flush=True)
    if LOG:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(str(msg) + "\n")


def main():
    try:
        _main()
    except SystemExit:
        raise
    except BaseException:
        import traceback

        tb = traceback.format_exc()
        log(tb)
        sys.exit(3)


def _main():
    import FreeCAD as App
    import Part
    import MeshPart
    import Mesh
    import Import
    src = os.environ.get("STEP2STL_SRC")
    dst = os.environ.get("STEP2STL_DST")
    if not src or not dst:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    lin_def = float(os.environ.get("STEP2STL_LIN_DEF", "0.5"))
    ang_def = float(os.environ.get("STEP2STL_ANG_DEF", "10.0"))

    t0 = time.time()
    log(f"Opening {src} ...")
    # Import.open is more robust than App.openDocument for CAD files run from
    # the command line (openDocument can fail with 'Invalid project file').
    Import.open(src)
    doc = App.getDocument("Unnamed")
    shapes = []
    for obj in doc.Objects:
        if hasattr(obj, "Shape") and not obj.Shape.isNull():
            try:
                if obj.Shape.Solids:
                    shapes.extend(obj.Shape.Solids)
                else:
                    shapes.append(obj.Shape)
            except Exception:
                pass

    if not shapes:
        log(f"ERROR: no solids found in {src}")
        sys.exit(2)

    log(f"Loaded {len(shapes)} solid(s) from {src}")

    import math

    angular = math.radians(ang_def)
    meshes = []
    total_faces = 0
    for i, shp in enumerate(shapes):
        bb = shp.BoundBox
        dims = [bb.XLength, bb.YLength, bb.ZLength]
        # Skip degenerate construction objects (datum planes / infinite solids)
        # that SolidWorks sometimes exports: huge or empty bounding boxes.
        if max(dims) > 1e6 or not shp.Faces:
            log(f"  skip solid {i} (degenerate bbox {dims[0]:.1f} x {dims[1]:.1f} x {dims[2]:.1f})")
            continue
        log(
            f"  mesh solid {i}: {dims[0]:.2f} x {dims[1]:.2f} x {dims[2]:.2f} mm "
            f"(points[{len(shp.Vertexes)}], faces[{len(shp.Faces)}])"
        )
        try:
            m = MeshPart.meshFromShape(
                Shape=shp,
                LinearDeflection=lin_def,
                AngularDeflection=angular,
                Relative=False,
            )
        except Exception as e:
            log(f"  solid {i} meshing failed: {e}; skipped")
            continue
        if m is not None and len(m.Facets):
            meshes.append(m)
            total_faces += len(m.Facets)

    if not meshes:
        log("ERROR: no valid meshes produced")
        sys.exit(2)

    # Register each mesh as a Mesh::Feature in the document: FreeCAD's
    # Mesh.export only accepts document objects, not raw MeshObject instances.
    features = []
    for i, m in enumerate(meshes):
        feat = doc.addObject("Mesh::Feature", f"MeshSolid{i}")
        feat.Mesh = m
        features.append(feat)

    os.makedirs(os.path.dirname(dst), exist_ok=True)
    Mesh.export(features, dst)
    log(
        f"Exported STL -> {dst}: {total_faces} triangles across {len(features)} solid(s), "
        f"deflection={lin_def} mm / {ang_def} deg, {time.time()-t0:.1f} s"
    )


if __name__ == "__main__":
    main()
else:
    # FreeCAD's freecadcmd does not set __name__ == "__main__" when running a
    # script, so call main() explicitly in that case too.
    main()
