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
| `anisotropic hyperelastic damage` | Fiber-reinforced neo-Hookean with damage | finite strain |
| `crystal plasticity` | Rate-dependent crystal plasticity | finite strain |

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

## Anisotropic hyperelastic damage

A compressible neo-Hookean matrix reinforced by any number of fiber families,
each phase degraded by a damage that grows with the largest strain energy it
has stored. This is LCM's AHD model:

```math
\begin{aligned}
W_m &= \tfrac{\kappa}{2}\left[\tfrac{1}{2}(J^2 - 1) - \log J\right]
      + \tfrac{\mu}{2}\left(J^{-2/3}\operatorname{tr}\mathbf{C} - 3\right), \\
W_f &= \frac{k}{q}\left[e^{q (I_4 - 1)^2} - 1\right],
      \qquad I_4 = \mathbf{M} \cdot \mathbf{C}\,\mathbf{M}, \\
d &= d_{\max}\left[1 - e^{-\alpha / s}\right],
      \qquad \alpha = \max_{\text{history}} W, \\
\mathbf{P} &= \phi_m (1 - d_m)\,\frac{\partial W_m}{\partial \mathbf{F}}
      + \sum_f \phi_f (1 - d_f)\,\frac{\partial W_f}{\partial \mathbf{F}}.
\end{aligned}
```

Without damage, the matrix is the `neohookean` model. As in LCM, the fibers
resist compression as well as tension. A fiber stretch at which
e^{q (I_4 − 1)^2} would overflow, which a long Newton trial step can reach,
fails the evaluation, and the step is cut.

```yaml
model: anisotropic hyperelastic damage
elastic modulus: 2.0e+05          # matrix
Poisson's ratio: 0.3
matrix volume fraction: 0.5
matrix maximum damage: 0.1
matrix damage saturation: 1.0e+04
fibers:
  - direction: [1.0, 0.0, 0.0]    # in the material frame
    k: 1.0e+04
    q: 100.0
    volume fraction: 0.5
    maximum damage: 1.0
    damage saturation: 1.0e+04
orientation:                      # optional; this one turns the fiber to y
  Bunge Euler angles: [90.0, 0.0, 0.0]
```

| Key | Default | Meaning |
|---|---|---|
| two elastic constants | — | matrix stiffness |
| `matrix volume fraction` | `1.0` | φ_m |
| `matrix maximum damage` | `1.0` | d_max of the matrix |
| `matrix damage saturation` | none | s of the matrix; without it the matrix does not damage |
| `fibers` | none | list of fiber families |
| `direction` | — | fiber direction M in the material frame (normalized) |
| `k`, `q` | `1.0` | fiber stiffness and exponential stiffening |
| `volume fraction` | `0.0` | φ_f |
| `maximum damage` | `1.0` | d_max of the fiber family |
| `damage saturation` | none | s of the fiber family; without it the family does not damage |
| `orientation` | identity | as for the [anisotropic elastic models](#Lattice-orientation) |

The stable explicit time step uses the stiffness of the undeformed material,
φ_m (λ + 2μ) + Σ_f 8 φ_f k. The exponential fibers stiffen fast with their
stretch: for q = 100 the largest wave modulus is 3 times the estimate at 5 %
fiber stretch and over 200 times at 10 %, so an explicit run with stretched
fibers needs a CFL that allows for it.

The internal variables are the largest energy and the damage of the matrix
and of each fiber family, `matrix_energy`, `fiber_1_energy`, …,
`matrix_damage`, `fiber_1_damage`, …. LCM keeps the energy of the previous
step instead of the largest over the history, so its damage recovers on
unloading; the two agree under monotonic loading.

## Crystal plasticity

Finite-deformation, rate-dependent crystal plasticity, ported from LCM's
CrystalPlasticity model. The deformation gradient splits as F = Fe Fp. A slip
system α has a unit direction s and normal n in the lattice frame; with the
lattice orientation Q its Schmid tensor is P_α = (Q s_α) ⊗ (Q n_α), and over
a step the plastic part follows the exponential map

```math
\mathbf{F}^p_{n+1} = \exp\Big(\sum_\alpha \Delta\gamma_\alpha \mathbf{P}_\alpha\Big)\,\mathbf{F}^p_n .
```

The elastic law is the [anisotropic Saint-Venant–Kirchhoff](#Anisotropic-elasticity)
law on Fe, and the resolved shear stress is the Mandel stress on the system,
τ_α = (Feᵀ Fe S) : P_α. Each system belongs to a slip family, whose flow rule
gives the slip rate from τ_α and the slip resistance g_α, and whose hardening
law gives the hardening state from the slip rates.

```yaml
model: crystal plasticity
symmetry: cubic
C11: 204600.0
C12: 137700.0
C44: 126200.0
orientation:
  Bunge Euler angles: [30.0, 45.0, 0.0]
integration scheme: implicit       # implicit (default) or explicit
residual: slip                     # slip (default) or slip and hardness
slip families:
  - slip systems: fcc
    flow rule:
      type: power law
      reference slip rate: 1.0
      rate exponent: 20.0
    hardening law:
      type: linear minus recovery
      hardening modulus: 355.0
      recovery modulus: 2.9
      initial hardening state: 122.0
```

The elastic constants and the orientation are those of the anisotropic
elastic models. `slip systems` is either the name of a built-in set, in the
lattice frame of LCM's tables, or a list of systems, each with a `direction`,
a `normal` and an optional `initial hardening state`:

| Set | Systems |
|---|---|
| `fcc` | 12 {111}⟨110⟩ |
| `bcc 110`, `bcc 112`, `bcc 123` | 12, 12 and 24 ⟨111⟩ systems on {110}, {112} and {123} |
| `hcp basal`, `hcp prismatic`, `hcp pyramidal a` | 3, 3 and 6 ⟨a⟩ systems |
| `hcp pyramidal c+a first order`, `hcp pyramidal c+a second order` | 12 and 6 ⟨c+a⟩ systems |

The hexagonal sets take the family key `c/a ratio` (default √(8/3)), with
lattice axis 1 along a₁ and axis 3 along c. `maximum slip rate` (default
ln(max Float64), as in LCM) bounds the slip rate of the family's systems;
a step that exceeds it is cut.

Every family needs a `flow rule` and a `hardening law`, each with its `type`
(`type: none` for a family that does not slip or harden); the keys of each
block are checked against those of its type. Flow rules (`type` and keys),
with x = τ/g:

| `type` | Slip rate | Keys |
|---|---|---|
| `power law` | γ̇₀ \|x\|^(m−1) x | `reference slip rate` γ̇₀, `rate exponent` m |
| `power law with drag` | the power law in series with γ̇₀ x/d | the power-law keys and `drag coefficient` d |
| `thermal activation` | γ̇₀ exp(−F₀ (1 − ξ^p)^q) sign τ, ξ = (\|τ\| − g)/s_t | `reference slip rate`, `activation energy` F₀, `thermal resistance` s_t, `p exponent`, `q exponent` |
| `none` | 0 | |

Hardening laws (`type` and keys); every law takes the `initial hardening
state`:

| `type` | Law | Keys |
|---|---|---|
| `linear minus recovery` | ġ = (H − R (g − g₀)) Σ\|γ̇\|, integrated exactly | `hardening modulus` H, `recovery modulus` R (default 0) |
| `saturation` | ġ_α = h₀ 2 Σ_β q_αβ \|γ̇_β\| (g_s − g_α)/(g_s − g₀), g_s = g_s0 (Σ\|γ̇\|/γ̇_s)^ω | `hardening rate` h₀, `initial saturation stress` g_s0, `saturation exponent` ω (default 0), `reference slip rate` γ̇_s (default 1) |
| `dislocation density` | ρ̇_α = max(k₁ √ρ^F_α − k₂ ρ_α, 0) \|γ̇_α\|, g_α = c μ b √ρ^P_α | `geometric factor` c, `generation factor` k₁, `annihilation factor` k₂, `shear modulus` μ, `burgers vector magnitude` b |
| `none` | g constant | |

The implicit scheme solves for the slips (or, with `residual: slip and
hardness`, for the slips and the hardening states) by Newton's method, with
the tolerances of an optional `local solver` block (`absolute tolerance`,
default 1e-12; `relative tolerance`, default 1e-14; `minimum iterations`,
default 0; `maximum iterations`, default 50). The explicit scheme evaluates
the resolved shear stress on Fp_n and the hardening with the slip rates of
the previous step; it is only stable while Δt (∂γ̇/∂τ) μ < 2. The implicit
scheme starts its Newton iteration from the slip rates of the previous step
(`slip predictor: rate`, the default) or from zero rates (`slip predictor:
none`). The explicit scheme needs `rate`: with `none` its hardening would
never change, as in LCM, and the input is rejected. With `residual: slip`,
as with the explicit scheme, the hardening is integrated by forward Euler,
as in LCM. For saturation hardening this overshoots the saturation stress
when a step slips more than about (g_s − g₀)/(2 h₀) in all, and a step that
crosses it is a local failure; `residual: slip and hardness` integrates the
hardening implicitly and takes such steps. A local failure cuts the global
step under [adaptive time stepping](time-integrators.md), and stops the run
without it. The first Newton iterate of a quasi-static step puts the whole
increment of a prescribed displacement on the elements next to the boundary,
so the time step must keep that strain within the range of the flow rule. The
tangent is the exact algorithmic tangent, which is not symmetric: use `linear
solver: direct` with the Hessian minimizer.

The internal variables are `Fp_ij`, then `gamma_α`, `gamma_dot_α`,
`tau_hard_α` (the hardening state: the slip resistance, or the dislocation
density) and `tau_α` for every system, then `eqps`.

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
- Anisotropic hyperelastic damage: `examples/materials/ahd/cube.yaml`
