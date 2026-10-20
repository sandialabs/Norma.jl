# LCM reference data

Reference ("gold") outputs of the elastic crystal plasticity tests and the
anisotropic hyperelasticity tests of [LCM](https://github.com/sandialabs/LCM),
the Albany-based code from which Norma's anisotropic elasticity and
anisotropic hyperelastic damage models are ported. The Norma tests compare
against these files, so they need neither an LCM build nor `exodiff`.

The data come from LCM commit `3e3af13` (2026-10-16), which includes the
crystal plasticity fixes of
[sandialabs/LCM#127](https://github.com/sandialabs/LCM/issues/127). The first
line of each file names its source: the gold file and the last LCM commit
that changed it. LCM is distributed under the BSD license reproduced in
`LICENSE-LCM.txt`.

## Regenerating

```bash
julia --project=. test/lcm-reference/extract.jl <path to an LCM checkout>
```

## Format

Each case has a comma-separated file. Lines starting with `#` are comments.

- `<case>.csv` holds the element variables, one row per
  `(step, block, element, point)`. When the integration points of an element
  agree to a relative 1e-9, they are collapsed into one row with `point = 0`.
  The comment line gives the largest relative spread between points.
- `materials/<case>.yaml` is the LCM material file of the case (the input
  file itself for the Material Point Simulator cases).
- `meshes/` holds LCM meshes that the tests read, with the block and node
  set names LCM gives them by default, `block_<id>` and `nodelist_<id>`.

The element files of the FEM cases also carry `F_nodal_k`, the deformation
gradient of each element from the affine fit of its nodal displacements,
which is exact for the homogeneous deformations of these cases, and a
comment line with the Exodus block names, which the materials files refer to.

Tensors are numbered row-major as 3 × 3 tensors in every file: `Cauchy_Stress_k`
is σ[(k - 1) ÷ 3 + 1, (k - 1) % 3 + 1]. The same holds for `F` and `Strain`.
The 2D case carries only the in-plane components 1, 2, 4 and 5.

The material point cases (`AHD-*`) write the step index as the time; their
time step is the `Step Size` of the LCM input, 0.01.

## Cases

The mesh of every FEM case is a unit cube of one HEX8 element, except the 2D
case (one QUAD4 element, plane strain).

| Case | What it exercises |
|---|---|
| `ElasticShear2d`, `ElasticShear3d` | Cubic elasticity, lattice aligned with the axes, simple shear |
| `ElasticShearRotated3d` | Cubic elasticity, lattice from Bunge angles (0.12, 0.03, 0.28) rad, simple shear |
| `AHD-uniaxial`, `AHD-shear`, `AHD-hydrostatic` | Material point: fiber-reinforced hyperelastic model with damage |

## History

Before LCM commit `42ba4a5`, LCM rotated the elasticity tensor by Qᵀ while it
rotated the slip systems by Q, where the columns of Q are the basis vectors
([sandialabs/LCM#127](https://github.com/sandialabs/LCM/issues/127)). Of the
cases here, only `ElasticShearRotated3d` has an orientation for which the two
differ; its gold file was regenerated with the fix, and the files here are
from after it.
