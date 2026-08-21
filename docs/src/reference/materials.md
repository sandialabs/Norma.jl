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

Every constitutive model draws its elastic response from exactly **two** elastic
constants. Supply any valid pair from the table below; Norma computes the
remaining constants. Supplying zero, one, or an unsupported combination aborts.

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
| `seth-hill` | Seth-Hill generalized hyperelastic | finite strain |
| `hencky` | Hencky (logarithmic strain) hyperelastic | finite strain |
| `j2 plasticity` | J2 (von Mises) plasticity | finite strain |

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

## Canonical examples

- Linear elastic (single domain): `examples/single/static-solid`
- J2 plasticity: `examples/materials/j2/cube.yaml`
- Seth-Hill: `examples/ems/awful-cube/awful-cube.yaml`
