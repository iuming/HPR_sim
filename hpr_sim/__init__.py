# -*- coding: utf-8 -*-
"""HPR (High-Pressure Rinsing) simulation of superconducting RF cavities.

A ray-casting based coverage / cleaning-efficiency simulator following the
methods from the reference material in `References/`:

  * CN121302684A (Dongguan/CSNS patent) - Möller-Trumbore ray-triangle
    intersection, nozzle as a bundle of rays from a common origin, cleaning
    intensity = cos(angle between jet and surface normal) > threshold,
    effective cleaning coverage ratio as the key metric.
  * document.pdf (FRIB) - STL based ray casting, cleaning intensity, missed
    coverage detection.
  * MOP51 / WEPS143 (RAON) - barycentric intersection, coverage by a
    translating + rotating nozzle wand.
  * MOPB070 / PIP-II - exposure time, nozzle-surface distance, angle of
    incidence as the dominant cleaning parameters.

Quick start
-----------
.. code-block:: console

    python -m hpr_sim --write-template configs/new_cavity.yaml   # scaffold
    python -m hpr_sim -f configs/main_cavity.yaml                # run a case
    python -m hpr_sim -f configs/main_cavity.yaml --sweep threshold,tilt --jobs 4
"""

__version__ = "0.2.0"
