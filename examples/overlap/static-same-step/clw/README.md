# Circular Laser Weld: Graded Overlapping Schwarz Pair

Overlapping Schwarz decomposition of a circular laser weld specimen with J2
plasticity, ported from the graded Schwarz pair of the FY26 ASC P&EM L2
milestone CLW study. This is the least-DOF configuration of that study:
fine-domain base-material rings graded to t/4 of the specimen thickness
(t = 0.0604 in), decomposition at 0.30 in center distance from the weld with
0.050 in overlap half-thickness. Coarse plate 14,553 hexes / 20,184 nodes,
fine weld domain 24,804 hexes / 31,672 nodes — 39,357 hexes total against
61,978 for the monolithic reference mesh. All units SI (meshes in meters).

Physics matching the study's Sierra/SM deck: 304L steel, finite-deformation
J2 plasticity (E = 200 GPa, nu = 0.25, rho = 7800 kg/m^3, yield 1 GPa, linear
hardening 20 GPa), quasistatic, 50 steps to t = 1, displacement-controlled
pull of 5 mm in y on the plate top surface with a cosine ramp, x/y symmetry
and a z-fixed line on the coarse plate, overlapping Schwarz coupling with
minimum 1 iteration and relative tolerance 1e-3 (the study's operating
point). Both stops of the smoke case converge in one Schwarz iteration.

## Files

- `clw-geo.jou` — shared laser weld geometry.
- `clw-schwarz-norma.jou` — Cubit journal generating both meshes, including
  the `schwarz_coarse_ss`/`schwarz_fine_ss` side sets the `Schwarz overlap`
  boundary condition couples through. Built with Cubit 17.08:
  `cubit -batch -nographics -nojournal -input clw-schwarz-norma.jou`
- `clw-coarse.g`, `clw-fine.g` — the generated meshes.
- `clw.yaml`, `clw-coarse.yaml`, `clw-fine.yaml` — the full 50-step case.
- `clw-smoke.yaml` — 2-step version for timing and setup checks.
- `nh-*.yaml` — neohookean (elastic) variant; much faster, useful for
  Schwarz-behavior studies that do not need plasticity.

## Running

```
JULIA_NUM_THREADS=8 julia --project=/path/to/Norma.jl /path/to/Norma.jl/src/Norma.jl clw.yaml
```

`use line search: true` in the solver blocks is required: the first Newton
iterate of a step yields spuriously in the boundary layer next to the loaded
surface, and raw full Newton bounces off the yield surface. The case also
needs the J2 state-commit fix (commit b15ace97, 2026-08-28); before it the
first step never converges.

Measured on 8 threads (Linux workstation, 2026-08-28): J2 about 38 s per step
after JIT warm-up (roughly 30-35 min for the 50-step case); neohookean about
1.5 s per step.
