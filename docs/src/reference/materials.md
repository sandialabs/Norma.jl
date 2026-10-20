# Materials

Materials are declared inside the [`model`](model.md) block. The `material`
subsection maps each mesh element block to a named material definition, and each
named definition selects a constitutive model and its parameters:

```yaml
model:
  type: solid mechanics
  material:
    blocks:
      cantilever: steel      # element block "cantilever" uses material "steel"
    steel:                   # the material definition named above
      model: linear elastic
      elastic modulus: 6.895e+09
      Poisson's ratio: 0.25
      density: 2768.0
```

The names on the left of `blocks` must match the element block names in the
Exodus mesh, and the number of entries must equal the number of element blocks.
The name on the right is an arbitrary label that must itself appear as a key
under `material`, holding that material's parameters.

## Elastic constants

Every isotropic constitutive model draws its elastic response from exactly
**two** elastic constants (the anisotropic models take the constants of their
symmetry class instead; see [Anisotropic elasticity](#Anisotropic-elasticity)).
Supply any valid pair from the table below; Norma computes the remaining
constants. Supplying zero, one, or an unsupported combination aborts.

| Key | Symbol | Constant |
|---|---|---|
| `elastic modulus` | E | Young's modulus |
| `Poisson's ratio` | ν | Poisson's ratio |
| `bulk modulus` | κ | bulk modulus |
| `Lamé's first constant` | λ | first Lamé parameter |
| `shear modulus` | μ | shear modulus |

Accepted pairs: `elastic modulus` with any one of the other four;
`Poisson's ratio` with `bulk modulus`, `Lamé's first constant`, or
`shear modulus`; `bulk modulus` with `Lamé's first constant` or `shear modulus`;
or `Lamé's first constant` with `shear modulus`. The key spellings are literal,
including the apostrophe in `Poisson's ratio` and the accented, apostrophized
`Lamé's first constant`.

`density` (ρ) is optional and defaults to `0.0`. It is not used by static
analyses but is required in practice by the dynamic time integrators, which need
a mass matrix.

## Constitutive models

The `model` key selects the constitutive law. The accepted strings are literal
and case-sensitive:

| `model` value | Model | Kinematics |
|---|---|---|
| `linear elastic` | Linear elastic | infinitesimal strain |
| `Saint-Venant Kirchhoff` | Saint-Venant–Kirchhoff | finite strain |
| `neohookean` | Neohookean | finite strain |
| `r-neohookean` | Reciprocal Neohookean | finite strain |
| `seth-hill` | Seth-Hill generalized hyperelastic | finite strain |
| `hencky` | Hencky (logarithmic strain) hyperelastic | finite strain |
| `j2 plasticity` | J2 (von Mises) plasticity | finite strain |
| `anisotropic linear elastic` | Anisotropic Hookean elasticity | infinitesimal strain |
| `anisotropic Saint-Venant Kirchhoff` | Anisotropic Saint-Venant–Kirchhoff | finite strain |

### Linear elastic

Small-strain Hookean elasticity.

```yaml
model: linear elastic
elastic modulus: 6.895e+09
Poisson's ratio: 0.25
density: 2768.0
```

Parameters: two elastic constants, optional `density`.

### Saint-Venant–Kirchhoff

Finite-strain extension of Hookean elasticity (linear relation between the
Green–Lagrange strain and the second Piola–Kirchhoff stress).

```yaml
model: Saint-Venant Kirchhoff
elastic modulus: 6.895e+09
Poisson's ratio: 0.25
density: 2768.0
```

Parameters: two elastic constants, optional `density`.

### Neohookean

Compressible neo-Hookean hyperelasticity.

```yaml
model: neohookean
elastic modulus: 6.895e+09
Poisson's ratio: 0.25
density: 2768.0
```

Parameters: two elastic constants, optional `density`.

### Seth-Hill

Generalized hyperelastic model with tunable strain-measure exponents `m` and
`n` (both required integers).

```yaml
model: seth-hill
elastic modulus: 6.895e+09
Poisson's ratio: 0.25
density: 2768.0
m: 2
n: 2
```

Parameters: two elastic constants, optional `density`, required integers `m` and
`n`.

### Hencky

Hyperelasticity quadratic in the logarithmic strain,

```math
\psi = \tfrac{\kappa}{2}\,(\operatorname{tr}\mathbf{E})^2
      + \mu\,\operatorname{dev}\mathbf{E} : \operatorname{dev}\mathbf{E},
\qquad \mathbf{E} = \tfrac{1}{2}\log\mathbf{C}.
```

Under uniaxial stress the lateral log stretch is exactly ``-\nu`` times the
axial one and the Kirchhoff stress is exactly ``E \log\lambda`` at any
stretch, which makes this model useful for testing strong nonlinearities
against closed-form answers.

```yaml
model: hencky
elastic modulus: 1.0e+09
Poisson's ratio: 0.25
density: 1000.0
```

Parameters: two elastic constants, optional `density`.

### J2 plasticity

Finite-deformation J2 plasticity with a multiplicative elastic–plastic split,
radial-return integration, and linear isotropic hardening.

```yaml
model: j2 plasticity
elastic modulus: 70.0e+09
Poisson's ratio: 0.25
density: 1000.0
yield stress: 250.0e+06
hardening modulus: 0.7e+09
```

| Key | Required | Default | Meaning |
|---|---|---|---|
| two elastic constants | yes | — | elastic response |
| `density` | no | `0.0` | mass density |
| `yield stress` | no | `0.0` | initial yield stress σ_y |
| `hardening modulus` | no | `0.0` | linear isotropic hardening modulus H |

## Anisotropic elasticity

The anisotropic models take their elastic constants in the frame of the
material's lattice (for a crystal, the cube axes; for a hexagonal crystal,
lattice axis 3 is the c axis), and an `orientation` that places the lattice in
the sample. `anisotropic linear elastic` is the small-strain law σ = ℂ : ε;
`anisotropic Saint-Venant Kirchhoff` is its finite-strain counterpart, linear
between the Green–Lagrange strain and the second Piola–Kirchhoff stress, and
the elastic law of LCM's crystal plasticity model.

```yaml
model: anisotropic Saint-Venant Kirchhoff
symmetry: cubic
C11: 168.4e+09
C12: 121.4e+09
C44: 75.4e+09
density: 8960.0
orientation:
  Bunge Euler angles: [30.0, 45.0, 0.0]   # degrees unless `angle unit: radians`
```

The constants are the components of the Voigt stiffness, in the order 11, 22,
33, 23, 13, 12 with engineering shear strains, so that `C44` relates σ₂₃ to
2ε₂₃. The `symmetry` key selects which of them are given:

| `symmetry` | Constants | Notes |
|---|---|---|
| `isotropic` | two elastic constants | as for the isotropic models |
| `cubic` | `C11`, `C12`, `C44` | |
| `hexagonal` | `C11`, `C12`, `C13`, `C33`, `C44` | C₆₆ = (C₁₁ − C₁₂)/2; `transversely isotropic` is a synonym |
| `tetragonal` | `C11`, `C12`, `C13`, `C33`, `C44`, `C66` | |
| `orthotropic` | `C11`, `C22`, `C33`, `C12`, `C13`, `C23`, `C44`, `C55`, `C66` | |
| `triclinic` | any of the 21 `Cij`, i ≤ j | the others are zero |

The constants must give a positive definite stiffness; otherwise the input is
rejected. For explicit dynamics the stable time step uses the largest wave
modulus ρc² over all propagation directions in place of λ + 2μ.

### Lattice orientation

The orientation Q holds the lattice axes, expressed in the sample frame, as its
columns. A lattice direction v has sample components Q v, and the elasticity
tensor in the sample frame is

```math
\mathbb{C}_{ijkl} = Q_{ip}\,Q_{jq}\,Q_{kr}\,Q_{ls}\,\mathbb{C}^{L}_{pqrs}.
```

The `orientation` block takes one of:

| Key | Value |
|---|---|
| `basis vectors` | the three lattice axes in the sample frame, `[[...], [...], [...]]` (LCM's `Basis Vector 1/2/3`); each is normalized, and the triad must be orthogonal and right-handed |
| `Bunge Euler angles` | `[φ1, Φ, φ2]`, giving Q = R_z(φ1) R_x(Φ) R_z(φ2) |
| `mesh variables` | names of element variables of the input mesh that give each element its own orientation: three names for Bunge Euler angles, or nine for Q row-major (Q₁₁, Q₁₂, …, Q₃₃) |

`angle unit` (`degrees`, the default, or `radians`) applies to the Bunge Euler
angles, also when they come from the mesh, and `time index` (default 1)
selects the step of the mesh variables. Without an `orientation`, the lattice
axes coincide with the sample axes.

The `basis vectors` follow the convention of LCM's crystal plasticity model.
LCM versions before commit `42ba4a5` rotated the elasticity tensor by Qᵀ
instead ([sandialabs/LCM#127](https://github.com/sandialabs/LCM/issues/127));
Norma reproduces purely elastic results of those versions with the inverse
orientation, the Bunge angles (−φ2, −Φ, −φ1).

## Material point driver

`Norma.material_point_history(material, deformation, times; orientation,
free, stress)` evaluates a constitutive model at one point without a mesh,
under F = `deformation(t)` at the given times, with the lattice orientation
Q of an oriented material. The components of F listed in `free` (a 3 × 3
Boolean matrix) are solved for so that the same components of the first
Piola–Kirchhoff stress follow `stress(t)` (zero by default).

## Canonical examples

- Linear elastic (single domain): `examples/single/static-solid`
- J2 plasticity: `examples/materials/j2/cube.yaml`
- Seth-Hill: `examples/ems/awful-cube/awful-cube.yaml`
- Anisotropic elasticity: `examples/single/static-solid/cube/anisotropic/cube-shear.yaml`
