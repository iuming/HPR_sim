# -*- coding: utf-8 -*-
# Author: Ming Liu (刘铭)
#         National Synchrotron Radiation Laboratory (NSRL),
#         University of Science and Technology of China (USTC)
# Contact: mliu1018@ustc.edu.cn
# SPDX-License-Identifier: MIT

"""HPR (High-Pressure Rinsing) simulation of superconducting RF cavities.

A ray-casting based coverage / cleaning-efficiency simulator following the
methods used in published SRF HPR practice:

  * CN121302684A - Möller-Trumbore ray-triangle intersection, nozzle as a
    bundle of rays from a common origin, cleaning intensity =
    cos(angle between jet and surface normal) > threshold, effective cleaning
    coverage ratio as the key metric.
  * FRIB (STL-based ray casting, cleaning intensity, missed-coverage
    detection), RAON (barycentric intersection, translating + rotating nozzle
    wand) and PIP-II simulations (exposure time, nozzle-surface distance,
    angle of incidence as the dominant cleaning parameters).

Quick start
-----------
.. code-block:: console

    python tools/make_tesla_cavity.py                           # example mesh
    python -m hpr_sim --write-template configs/new_cavity.yaml  # scaffold
    python -m hpr_sim -f examples/tesla_9cell.yaml              # run a case
    python -m hpr_sim -f examples/tesla_9cell.yaml --sweep threshold,tilt --jobs 4
"""

__version__ = "0.2.0"
