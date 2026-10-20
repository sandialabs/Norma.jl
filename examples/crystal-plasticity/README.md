# Crystal plasticity polycrystals

Two polycrystals whose grains carry their own lattice orientation. The
orientation of every element is an element variable of the mesh, read with

```yaml
orientation:
  mesh variables: [phi1, Phi, phi2]
```

where `phi1`, `Phi` and `phi2` are Bunge Euler angles in degrees (see
[Lattice orientation](../../docs/src/reference/materials.md#lattice-orientation)).
`test/polycrystal.jl` runs both and checks them against the Voigt and Reuss
bounds in elasticity, and the RubiksCube with one element per grain against
the Taylor and Sachs estimates in plasticity; the plastic runs of the Voronoi
polycrystal are the examples below, which no test runs. It also runs
`compare_estimates.jl` and the converter's voxel path.
`test/periodic.jl` runs them as periodic cells, with the
[`Periodic`](../../docs/src/reference/boundary-conditions.md#periodic)
boundary condition.

## `rubiks-cube/`: LCM's RubiksCube

The cube [−0.5, 0.5]³ split into 3 × 3 × 3 grains with the 27 random lattice
orientations and the FCC crystal of LCM's `tests/LCM/CrystalPlasticity/RubiksCube`,
pulled along x to 0.5 % strain at 1/s on rollers (`rubiks-cube.yaml`).

`make_rubiks_cube.jl` builds the meshes from an LCM checkout, with n × n × n
HEX8 elements per grain:

```bash
julia --project=../../.. make_rubiks_cube.jl /path/to/LCM 1 2
```

`rubiks-cube-1.g` has one element per grain, like LCM's `RubiksCube.g`, and
`rubiks-cube-2.g` has eight. Every grain is a block named `block_<k>` after
LCM's blocks; the faces are node sets `nsx-`, …, `nsz+` and side sets `ssx-`,
…, `ssz+`.

## `polycrystal/`: Voronoi polycrystals from MicroStructPy

A periodic Voronoi tessellation of the unit cube into 25 grains, generated
with [MicroStructPy](https://docs.microstructpy.org) as tetrahedra
(`polycrystal-tets.xml`) or as voxels (`polycrystal-voxels.xml`). The
periodic mesher and the `periodic_margin` setting of these inputs are not in
a release of MicroStructPy: release 1.5.9 on PyPI stops at `periodic_margin`.
The meshes were made with commit `46c3179b` (2026-10-13) of the `develop`
branch of [kip-hart/MicroStructPy](https://github.com/kip-hart/MicroStructPy),
which still reports version 1.5.9, with its dependency pyvoro-rimoli 1.4.1;
that commit reproduces both meshes byte for byte:

```bash
pip install pyvoro-rimoli==1.4.1 git+https://github.com/kip-hart/MicroStructPy@46c3179b7ebf1de31e9b9355d33b6871f77d2eb2
microstructpy polycrystal-tets.xml
julia --project=../../.. microstructpy_to_exodus.jl polycrystal-tets/trimesh.txt polycrystal-tet4.g
microstructpy polycrystal-voxels.xml
julia --project=../../.. microstructpy_to_exodus.jl polycrystal-voxels/rastermesh.txt polycrystal-hex8.g
```

`microstructpy_to_exodus.jl` draws a lattice orientation for every grain,
uniformly over SO(3) (`--seed N`, default 1), and writes it with the grain
number as element variables. It writes one block named `polycrystal`, or with
`--blocks` one block per grain, and with `--tet10` quadratic tetrahedra. The
node and side sets are those of the RubiksCube meshes, and a node set holds
each corner of the cube, `nsx-y-z-` to `nsx+y+z+`. The converter checks that
the mesh conforms across its periodic faces, and the node set of the upper
face of a periodic axis lists the images of the nodes of the lower face in
order.

`polycrystal-tension.yaml` pulls the voxel polycrystal like the RubiksCube,
on rollers, and `polycrystal-tension-tet4.yaml` the tetrahedral one.
`polycrystal-periodic.yaml` loads the voxel polycrystal as a periodic cell
instead, in uniaxial stress along x: F̄ₓₓ is prescribed, and every other
component of the average first Piola–Kirchhoff stress is zero (see
[`Periodic`](../../docs/src/reference/boundary-conditions.md#periodic)). With
`input mesh file: polycrystal-tet4.g` it runs the tetrahedral mesh. All of
them take LCM's steps of 5 × 10⁻⁴, halved when they fail (see
[Time steps](#time-steps)). At 0.5 % strain, the axial stresses compare with
the estimates of the same 25 grains as follows:

| Mesh and conditions | Steps halved | Sachs | FEM | Taylor |
|---|---|---|---|---|
| voxels, rollers (`polycrystal-tension.yaml`) | 5 | 276.2 | 341.2 | 373.2 |
| voxels, periodic (`polycrystal-periodic.yaml`) | 2 | 276.2 | 346.2 | 372.1 |
| tetrahedra, rollers (`polycrystal-tension-tet4.yaml`) | 6 | 276.2 | 343.0 | 373.5 |
| tetrahedra, periodic | 2 | 276.2 | 347.9 | 372.3 |

The Taylor estimate frees the same components of F̄ as the conditions: F̄_yy
and F̄_zz on rollers, and the five stress-free components of the periodic
cell.

## Time steps

A quasi-static step starts from the previous solution with the new boundary
values. Without a predictor, its first Newton iterate puts the whole increment
of the loaded face on the layer of elements next to it. With crystal
plasticity that strain must stay within the range of the flow rule (a slip
rate below the maximum of its family), or the step fails and adaptive time
stepping halves it: for the RubiksCube crystal, once it flows, the increment
Δu/h of that strain, with h the element size, must stay below about 5 × 10⁻⁴.
The steps would then shrink with the elements next to the faces: on the
tetrahedral mesh the strain that the first iterate puts on them is the
increment over a size as small as 0.013 (the largest |Σ ∂N_a/∂x| of the face
nodes, 78), so the steps would have to stay below about 6.4 × 10⁻⁶. The
examples therefore set
`predictor: tangent` (see
[Predictor](../../docs/src/reference/time-integrators.md#predictor)), which
spreads the increment of the face over the body, and take LCM's steps of
5 × 10⁻⁴, halved when they fail. The steps halved on the way to 0.5 % strain,
with each predictor:

| Mesh | `constant` | `linear` | `tangent` |
|---|---|---|---|
| `rubiks-cube-1.g`, 27 HEX8 | 9 | 1 | 1 |
| `rubiks-cube-2.g`, 216 HEX8 | 49 | 5 | 4 |
| `polycrystal-hex8.g`, 1000 HEX8 | 81 | 7 | 5 |

`linear` also halves the first step on all but the coarsest mesh, as it has
no increment to extrapolate there. With `tangent`, the tetrahedral
polycrystal reaches 0.5 % strain in its 10 stops with 6 steps halved. The
larger steps change the response little: the RubiksCube of `rubiks-cube-2.g`
carries 335.0 MPa at 0.5 % strain, against 335.2 MPa with steps of
6.25 × 10⁻⁵ and no predictor.

With a predictor, the steps halved are those at the onset of flow, between
strains of about 0.05 and 0.25 %. The predictor is not their cause: the local
integration of the flow rule starts from the slip rates of the previous step,
which are still small, and at some integration points its first iterate
exceeds the maximum slip rate. Of the RubiksCube of `rubiks-cube-2.g`, even
the converged displacement at 0.15 % strain of a run with steps of
1.25 × 10⁻⁴ fails to evaluate, at 89 of its 1728 integration points, with one
step of 5 × 10⁻⁴ from the state at 0.1 %.

A periodic cell needs no predictor. The
[periodic condition](../../docs/src/reference/boundary-conditions.md#periodic)
applies the increment of F̄ to the whole cell before the first iterate, which
is the solution of a homogeneous cell. The voxel and tetrahedral cells of
`polycrystal-periodic.yaml` each halve 2 of their steps of 5 × 10⁻⁴, both at
the onset of flow. With `predictor: tangent`, which then predicts in place of
the condition, they halve the same steps and reach the same stresses.

`compare_estimates.jl` runs an example and lists its axial stress next to the
Sachs and Taylor estimates of the same grains, on rollers or for a periodic
cell. It integrates the estimates in steps of at most 5 × 10⁻⁵ strain between
the stops of the run:

```bash
cd rubiks-cube
julia --project=../../.. ../compare_estimates.jl rubiks-cube.yaml
cd ../polycrystal
julia --project=../../.. ../compare_estimates.jl polycrystal-periodic.yaml
```
