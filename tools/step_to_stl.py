# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""
Python wrapper around the FreeCAD headless STEP -> STL conversion
(`tools/convert_step_to_stl.py`).

Locates `freecadcmd.exe`, runs the conversion, and returns the STL path.
"""
from __future__ import annotations

import os
import shutil
import subprocess

_TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "convert_step_to_stl.py")


def find_freecadcmd():
    """Locate freecadcmd.exe."""
    for p in (os.environ.get("FREECAD_BIN"),
              r"D:\Software\FreeCAD 1.1\bin\freecadcmd.exe"):
        if p and os.path.isfile(p):
            return p
    p = shutil.which("freecadcmd")
    if p:
        return p
    raise FileNotFoundError(
        "freecadcmd.exe not found; set FREECAD_BIN or install FreeCAD")


def step_to_stl(src_step: str, dst_stl: str, lin_def: float = 0.5,
                ang_def: float = 10.0, log_path: str = None, quiet: bool = True):
    """Convert a STEP file to an STL mesh via FreeCAD.

    Returns the STL path.  Raises RuntimeError on failure.
    """
    src_step = os.path.abspath(src_step)
    dst_stl = os.path.abspath(dst_stl)
    if not os.path.exists(src_step):
        raise FileNotFoundError(src_step)

    os.makedirs(os.path.dirname(dst_stl), exist_ok=True)
    exe = find_freecadcmd()
    env = dict(os.environ)
    env["STEP2STL_SRC"] = src_step
    env["STEP2STL_DST"] = dst_stl
    env["STEP2STL_LIN_DEF"] = str(lin_def)
    env["STEP2STL_ANG_DEF"] = str(ang_def)
    if log_path:
        env["STEP2STL_LOG"] = os.path.abspath(log_path)

    cfg_dir = os.path.dirname(dst_stl)
    cmd = [exe, "-u", os.path.join(cfg_dir, "fc_user.cfg"),
           "-s", os.path.join(cfg_dir, "fc_system.cfg"), _TOOL]
    r = subprocess.run(cmd, env=env, cwd=os.path.dirname(exe),
                       capture_output=True, text=True, timeout=1800)
    if r.returncode != 0 or not os.path.exists(dst_stl):
        tail = (r.stdout or "")[-2000:] + (r.stderr or "")[-1000:]
        if log_path and os.path.exists(log_path):
            with open(log_path, "r", encoding="utf-8", errors="replace") as fh:
                tail += "\n---log---\n" + fh.read()[-2000:]
        raise RuntimeError(f"STEP->STL conversion failed for {src_step}:\n{tail}")
    if not quiet:
        print(f"STEP->STL: {src_step} -> {dst_stl}")
    return dst_stl
