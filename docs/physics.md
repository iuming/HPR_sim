# Physics model

This page summarises the models implemented in `hpr_sim` and the literature
they come from. A full Chinese survey with all references is kept in the
original project notes (`References/HPR文献调研_液滴阻力与射流衰减.md`).

## 1. Cleaning criterion (CN121302684A; FRIB / PIP-II / RAON practice)

The core idea used by modern SRF HPR simulations is **ray casting**: the
nozzle head is treated as a bundle of rays (one per hole × spray direction),
and the cleaning intensity at the wall is

```
I = |cos φ|
```

with φ the angle between a ray and the surface normal. A wall point is
*effectively cleaned* only if

```
I >= clean_cos_thresh            (default 0.5)
d_min <= d <= d_max              (usable stand-off window)
(optional) p_impact >= min_impact_bar
```

The patent CN121302684A (Simulation system for the high-pressure rinsing of
SRF cavities) and the FRIB/PIP-II/RAON simulations all express coverage as
the fraction of the RF surface cleaned this way.

## 2. Jet model

### Exit conditions (Bernoulli)

```
v0 = C_d sqrt(2 ΔP / ρw)
Q  = C_d A_nozzle v0            (per hole)
We = ρa v0² d_drop / σ          (Weber number at the exit)
```

With the defaults (100 bar, d=0.5 mm, C_d=0.7) this gives
`v0 ≈ 99 m/s`, `Q/jet ≈ 19.4 L/min`, `We ≈ 81`.

### Intact core and breakup

The compact jet keeps its initial velocity until the breakup length
(dimensionless breakup coefficient `k_b ≈ 12`, e.g. Sallam–Dai–Faeth 2002):

```
L_breakup = k_b D sqrt(ρw/ρa)
```

For the reference jet `L_breakup ≈ 65 mm`. Up to this distance the jet is
treated as a solid column hitting the wall; beyond it, the spray is modelled
as droplets of initial size (Weber-based estimate, e.g. Kooij et al.)

```
d_drop = We_c σ / (ρa v0²)
```

with `We_c ≈ 12`, giving `d_drop ≈ 74 µm`.

### Droplet deceleration

Loss of axial velocity due to air drag is integrated along the flight:

```
dv/dx = - 3 ρa C_D(Re) v / (8 ρw r_d)
```

with the **Reynolds-number-dependent drag coefficient**; four correlations
are shipped (`jet.drop_cd_mode`):

| mode | formula | reference |
|---|---|---|
| `schiller-naumann` (default) | `C_D = 24/Re (1 + 0.15 Re^0.687)` | Schiller & Naumann (1933) |
| `morsi-alexander` | piecewise polynomial | Morsi & Alexander, *JFM* 55 (1972) 355 |
| `putnam` | `24/Re (1 + Re^2/6)^...` approximation | Putnam (1961) |
| `constant` | fixed `C_D` | — |

The 1/e velocity distance of the default jet is ≈ 121 mm, far less than the
`linear` (exponential) model's simple decay coefficient — a real difference
that matters for which surfaces pass the distance window.

### Wall impact pressure

At the wall the *stagnation* pressure of the local jet/droplet velocity is
reported per face (`mean_impact_bar_hit`, and required via
`sim.min_impact_bar`):

```
p_impact = ½ ρw v_loc²
```

Because the intact core is short, far surfaces (e.g. the far side of each
cell equator) are reached with considerably reduced force — this is exactly
what the `impact` sweep exposes.

## 3. Ground-truth benchmarks

* **Guha et al., J. Mater. Process. Technol. 211 (2011) 610** — measured jet
  stand-off experiments on high-pressure water; supports a usable distance
  window and motivates `dmax_kd` (stand-off `d_max = k_d D`; optimal ≈ 5D,
  critical ≈ 26D).
* **Aune et al., PRST-AB 3 (2000) 092001** — the TESLA cavity shape used in
  `tools/make_tesla_cavity.py` (equator 103.3 mm, mid iris 35 mm, end iris
  39 mm, half-cell 57.7 / 56 mm).

## 4. Mesh & ray casting

* Möller–Trumbore intersection, batched in NumPy.
* Two-stage candidate culling: sphere around the camera → a 10° cone per ray
  → per-triangle angular allowance.
* Normals are not required to be consistent (`meshcheck` warns, rays use
  `|cos φ|`); a purely water-tight mesh is still preferred because the
  nearest-hit logic assumes one closed surface.

## 5. Limitations

* Droplets are monodisperse (no size distribution), spherical, and
  non-deforming; no secondary breakup or coalescence.
* No CFD of the near-wall flow, and no surface-particle adhesion model —
  `min_impact_bar` is a corrosion-science-motivated threshold (validation
  pending).
* The mesh is assumed to be in millimetres and to represent the **inner**
  niobium surface.
