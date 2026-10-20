# LCM reference data

Reference ("gold") outputs of the crystal plasticity and anisotropic
hyperelasticity tests of [LCM](https://github.com/sandialabs/LCM), the
Albany-based code from which Norma's anisotropic elasticity and crystal
plasticity models are ported. The Norma tests compare against these files, so
they need neither an LCM build nor `exodiff`.

The data come from LCM commit `3e3af13` (2026-10-16), which includes the
crystal plasticity fixes of
[sandialabs/LCM#127](https://github.com/sandialabs/LCM/issues/127). The first
line of each file names its source: the gold file and the last LCM commit
that changed it, or the output of an LCM run and the LCM commit of the run.
LCM is distributed under the BSD license reproduced in `LICENSE-LCM.txt`.

## Regenerating

```bash
julia --project=. test/lcm-reference/extract.jl <path to an LCM checkout> [--runs <folder of LCM runs>]
```

LCM ships no gold file for its polycrystal, RubiksCube. The `RubiksCube`
cases are therefore taken from runs of LCM `3e3af13` (Trilinos `373b9b89`,
serial GCC 13.4 release build), with `--runs` pointing at the folder that holds
the runs' outputs and material files:

- `RubiksCube`: `RubiksCube_Implicit.yaml` as LCM ships it.
- `RubiksCube_NoAverage`: `inputs/RubiksCube_Implicit_NoAvg.yaml`, the same
  input with `Weighted Volume Average J` and `Volume Average Pressure` set to
  false in every block of its materials.
- `RubiksCubeTet10`: `inputs/RubiksCubeTet10.yaml`, the latter on
  `meshes/RubiksCubeTet10.g` with `Use Composite Tet 10` in every block.

The materials of each run are `materials/<case>.yaml`, which the inputs call
`<input>_Materials.yaml`. These cases keep only the steps and elements their
test reads (see [Cases](#Cases)); `--every-step` writes all of them instead,
with their nodal files.

## Format

Each case has up to three comma-separated files. Lines starting with `#` are
comments.

- `<case>.csv` holds the element variables, one row per
  `(step, block, element, point)`. When the integration points of an element
  agree to a relative 1e-9, they are collapsed into one row with `point = 0`.
  The comment line gives the largest relative spread between points.
- `<case>-nodal.csv` holds the nodal coordinates and nodal variables, one row
  per `(step, node)`, for the two cases whose test reads them,
  `RubiksCube_NoAverage` and `RubiksCubeTet10`.
- `<case>-times.csv` gives the time of every LCM step, for the cases that
  keep only some steps; a comment line in `<case>.csv` lists those it keeps.
- `materials/<case>.yaml` is the LCM material file of the case (the input
  file itself for the Material Point Simulator cases).
- `meshes/` holds the LCM meshes the FEM tests run (`single-crystal-lcm.jl`,
  `lcm-rubiks-cube.jl`),
  with the block and node set names LCM gives them by default, `block_<id>`
  and `nodelist_<id>`.

The element files of the FEM cases also carry `F_nodal_k`, the deformation
gradient of each element from the affine fit of its nodal displacements,
which is exact for the homogeneous deformations of these cases, and a
comment line with the Exodus block names, which the materials files refer to.
`lcm_crystal_plasticity.jl` translates the LCM materials to Norma's input and
replays every case at the material point along LCM's deformation and times.

Tensors are numbered row-major as 3 × 3 tensors in every file: `Cauchy_Stress_k`
is σ[(k - 1) ÷ 3 + 1, (k - 1) % 3 + 1]. The same holds for `F`, `Fp`,
`Velocity_Gradient_Plastic` (Lp) and `Strain`. The 2D cases carry only the
in-plane components 1, 2, 4 and 5. Per-slip-system variables are
`gamma_<s>`, `gamma_dot_<s>`, `tau_<s>` (resolved shear stress) and
`tau_hard_<s>` (slip resistance, or dislocation density for the Dislocation
Density hardening law).

The material point cases (`CP-uniaxial`, `AHD-*`) write the step index as the
time; their time step is the `Step Size` of the LCM input (0.001 for
`CP-uniaxial`, 0.01 for `AHD-*`).

## Cases

The mesh of every FEM case is a unit cube of one or two HEX8 elements, except
the 2D cases (one QUAD4 element, plane strain) and RubiksCube (27 HEX8
elements, one per grain).

| Case | What it exercises |
|---|---|
| `ElasticShear2d`, `ElasticShear3d` | Cubic elasticity, lattice aligned with the axes, simple shear |
| `ElasticShearRotated3d` | Cubic elasticity, lattice from Bunge angles (0.12, 0.03, 0.28) rad, simple shear |
| `PlasticShear2d` | Two in-plane slip systems, linear power law, no hardening |
| `SingleSlip_Explicit`, `SingleSlip_Implicit` | Single slip, two blocks with different orientations, power law, no hardening |
| `SingleSlipHard_Explicit`, `SingleSlipHard_Implicit` | Single slip with linear hardening |
| `SingleSlipSaturation` | Single slip with saturation hardening, slip-hardness residual |
| `MultiSlipPlaneHard_Implicit` | FCC, 12 systems, Linear Minus Recovery hardening, adaptive steps |
| `MultiSlipPlaneHard_Implicit_Active_Sets` | The same with only the 8 active systems, fixed steps |
| `MinisolverStep_*` | The multi-slip problem with LCM's local solvers (Newton, Newton with line search, trust region; conjugate gradient and regularized line search are not registered LCM tests). No test reads `MinisolverStep_ConjugateGradient`: LCM writes only every fifth of its steps, which a material point cannot replay. It is kept as LCM's reference for its conjugate-gradient local solver |
| `ThermallyActivatedSlip` | BCC, 12 systems, thermally activated flow, saturation hardening |
| `MultiFamily` | Three slip families with different hardening laws |
| `DislocationDensityHardening` | FCC rotated 45° about z, dislocation-density hardening, constant true strain rate |
| `DislocationDensityTwoFamilies` | The same 12 systems as a second slip family, after an inactive first family; reproduces `DislocationDensityHardening` |
| `QuadSlipDislocationDensityTraction` | Loading along [011], traction controlled |
| `RubiksCube`, `RubiksCube_NoAverage` | 27 FCC grains, Linear Minus Recovery hardening, pulled along x on rollers with LCM's adaptive steps; with and without LCM's volume averaging of J and pressure. Not homogeneous: the 8 points of every element are kept, in LCM's order, and LCM numbers the output nodes its own way. `RubiksCube` keeps the last step of the 9 elements on the loaded face x = 0.5, which give the axial reaction; `RubiksCube_NoAverage` keeps every 4th of its 31 steps and the last, which include the steps of its largest deviations from Norma |
| `RubiksCubeTet10` | RubiksCube with every cube cut into 6 TETRA10 (`make_rubiks_cube_tet10.jl`), LCM's composite tet at cubature degree 2, at LCM's 74 steps; every 15th step and the last are kept |
| `CP-uniaxial` | Material point: FCC, single slip, uniaxial strain |
| `AHD-uniaxial`, `AHD-shear`, `AHD-hydrostatic` | Material point: fiber-reinforced hyperelastic model with damage |

## Composite tet 10

`composite-tet/` holds the references of `test/composite-tet10-lcm.jl` for
LCM's composite tetrahedron (`Use Composite Tet 10: true`). They are written by:

```bash
julia --project=. test/lcm-reference/extract_composite_tet.jl <LCM checkout> <runs>
```

| Case | What it exercises |
|---|---|
| `ProjectTet10` | LCM's gold: a cube of 6 composite tets stretched by 1 %, linear elastic, 4 points, with LCM's nodal projection of the linear field of its test (`proj_nodal_linear`). Converged to LCM's residual tolerance of 1e-10 only. Its projected Cauchy stress is not kept: unlike its point stress, it is not uniform |
| `ClampedCt_Explicit` | LCM's ClampedSDBC1DCt test: a bar of 6000 composite tets with a Gaussian initial displacement, 5 points, explicit Newmark with the lumped analytic mass, 150 steps. The nodes on the edge x = y = 0 at steps 0, 50, 100 and 150 |
| `ClampedCt_Implicit` | The same with implicit Newmark (β = 1/4) and the consistent analytic mass, 15 steps. The same nodes at steps 0, 1, 3, 5, 10 and 15. LCM starts from a zero acceleration |
| `Cylinder` | A cylinder of 387 curved composite tets (`inputs/make_cylinder.py`, gmsh 4.15), Neohookean, clamped ends, stretched by 30 % in 10 steps, 4 points. The displacements at 3 steps, and σ and F at every point at the last |

- **ProjectTet10.** It comes from the gold file in the LCM checkout. The node
  sets of its mesh get LCM's default names, `nodelist_<id>`.
- **Other cases.** They come from runs of LCM `3e3af13` on the inputs in
  `composite-tet/inputs`. The clamped inputs are LCM's own, written at every
  step.
- **Clamped mesh.** `meshes/clamped-ct.g` is LCM's `clamped_ct.g` with a node
  set of every node, `NodeSetAll`, for the initial condition.
- **Tensor layout.** LCM writes the component (i, j) of a tensor at point q of
  a tetrahedron as `<field>_<j>_<k>`, with k = 3 (q − 1) + i on two digits:
  each point's components are in column-major order. The files here are
  row-major, like the others.

## History

Before LCM commit `42ba4a5`, LCM rotated the elasticity tensor by Qᵀ while it
rotated the slip systems by Q, where the columns of Q are the basis vectors
([sandialabs/LCM#127](https://github.com/sandialabs/LCM/issues/127)). Of the
cases here, only `ElasticShearRotated3d` has an orientation for which the two
differ; its gold file was regenerated with the fix, and the files here are
from after it.
