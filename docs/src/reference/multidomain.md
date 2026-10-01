# Multidomain and Schwarz coupling

A multidomain simulation couples two or more single-domain models with the
Schwarz alternating method. It has two parts: a **top-level controller file**
(`type: multi`) that lists the subdomains and governs the coupling, and one
**subdomain input file** per domain — each an ordinary single-domain input with
Schwarz coupling boundary conditions added.

```yaml
type: multi
domains: ["cuboid-1.yaml", "cuboid-2.yaml"]
initial time: 0.0
final time: 1.0
time step: 1.0e-3
minimum iterations: 1
maximum iterations: 16
relative tolerance: 1.0e-12
absolute tolerance: 1.0e-08
```

## Controller file (`type: multi`)

### Domains and time window

| Key | Required | Meaning |
|---|---|---|
| `domains` | yes | list of subdomain input files, each a full single-domain input |
| `initial time` | yes | coupled simulation start time |
| `final time` | yes | coupled simulation end time |
| `time step` | yes | controller (Schwarz coupling) time step |

The controller supplies the time window to every subdomain, overriding any
`initial time`/`final time` in the subdomain files. Each subdomain keeps its own
time integrator and its own `time step`, which is clamped to at most the
controller step; a subdomain with a smaller step subcycles within each
controller step.

The controller file does **not** contain `model`, `solver`, or `time integrator`
blocks — those live in the subdomain files.

### Schwarz iteration control

Each controller step performs Schwarz iterations until the interface converges.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `minimum iterations` | yes | — | minimum Schwarz iterations per step |
| `maximum iterations` | yes | — | maximum Schwarz iterations per step |
| `absolute tolerance` | yes | — | absolute interface convergence tolerance |
| `relative tolerance` | yes | — | relative interface convergence tolerance |
| `constraint absolute tolerance` | no | `0` | absolute tolerance of the constrained Dirichlet–Neumann criterion on the root mean square interface jump, in the units of the constrained quantity (m/s or m); the controller's `absolute tolerance` does not apply to constrained pairs |
| `unconverged step action` | no | `warn` | what to do with a step that exhausts `maximum iterations` without meeting either tolerance: `warn` and continue, or `abort` |
| `stalled interface jump action` | no | `warn` | what to do when the interface jump of a paired impedance condition or the interface residual of a constrained Dirichlet–Neumann pair stops decreasing above `relative tolerance`: `warn` and accept the iterate, or `abort` |

A step that reaches `maximum iterations` without converging reports the errors
it stopped at against both tolerances. Under the default it says so and the run
continues, which is the historical behavior made visible; `abort` stops the run
instead, for cases where an unconverged interface must not be carried forward.

### Relaxation and acceleration

| Key | Required | Default | Meaning |
|---|---|---|---|
| `relaxation` | no | fixed | `aitken recursive` (Irons–Tuck) or `aitken secant` adaptive relaxation; omit for a fixed factor |
| `relaxation parameter` | no | `1.0` | relaxation factor θ; the constant factor under fixed relaxation, and under either Aitken method the factor used wherever an adaptive one is not yet available |
| `aitken N0 parameter` | no | `1` | Schwarz iteration, counting from zero, at which the adaptive factor takes over from `relaxation parameter` (Aitken methods only) |
| `interface predictor` | no | `false` | extrapolate the interface state at the start of each step |
| `naive stabilized` | no | `false` | naive interface stabilization |

Relaxation is not tied to a particular transmission condition: it applies to
whatever datum a coupling boundary condition transmits. For `Schwarz overlap`
and `Schwarz DN nonoverlap` the relaxed quantity is the interface displacement;
for `Schwarz impedance nonoverlap` and `Schwarz RR nonoverlap` it is the
interface force right-hand side. Both Aitken forms work
with all of them, and `relaxation parameter` and `aitken N0 parameter` mean the
same thing in each. On the impedance and Robin conditions the acceleration is
substantial: on the cantilever benchmark either Aitken form converges in about
a tenth of the Schwarz iterations that a fixed factor needs.

`relaxation parameter` is not ignored when `relaxation` names an Aitken method.
It is the factor applied for Schwarz iterations below `aitken N0 parameter`,
and the fallback whenever an adaptive factor cannot be formed, which includes
the first iterations of every step before two iterates exist to compare.

Aitken acceleration is applied only to stops with a single substep, that is
when the controller `time step` equals the relaxed subdomain's time step. A
windowed stop couples all of its time slots in one Schwarz iteration, where the adaptive
factors were measured to diverge or to lose to a fixed factor, so such stops
use `relaxation parameter` throughout. This is automatic and needs no input.

If a relaxation factor ever becomes small enough to leave an interface iterate
unchanged, the Schwarz iteration carries no information: every subdomain re-solves against
the data it already had and returns the solution it already had. The
displacement-based convergence test cannot distinguish that from convergence,
so such an iteration is refused as evidence of convergence and the run logs
`Relaxation factor near zero froze an interface iterate`. Seeing that message
repeatedly means the coupling is not advancing; check `relaxation parameter`.

### Output

| Key | Required | Default | Meaning |
|---|---|---|---|
| `Exodus output interval` | no | the controller `time step` | applied uniformly to all subdomains |
| `CSV output interval` | no | `0.0` (disabled) | applied uniformly to all subdomains |
| `blended energy output` | no | `false` | write an Arlequin-blended kinetic/stored/total-energy CSV each stop, removing the double count in overlapping regions |

The blending weights depend only on the reference configuration and are
computed once, with grid searches over the partner elements and the Schwarz
side sets, so enabling the energy output costs a few seconds of setup and a
small fraction of a step per stop.

The file `<name>-energy.csv` has one row per stop. Its columns are, in order:

| Column | Meaning |
|---|---|
| `time` | time of the stop |
| `stored_energy`, `kinetic_energy`, `total_energy` | blended physical energy E1: strain energy, kinetic energy (the staggered form ½ vᵀ M_L v − (Δt²/8) aᵀ M_L a for central difference subdomains), and their sum |
| `e2_total` | E2, the sum over subdomains of the energy of the system differentiated in time |
| `e2_<subdomain>` | E2 of one subdomain, one column per subdomain in the order of `domains` |
| `staggered_kinetic_interface`, `staggered_kinetic_interior` | staggered kinetic energy of the central difference subdomains, summed over the rows of the nodes in any Schwarz side set and over the remaining rows; `NaN` without central difference subdomains |
| `displacement_jump_<D>_<N>`, `velocity_jump_<D>_<N>` | one-sided interface jumps of each Dirichlet–Neumann pair at the end of the stop, D the Dirichlet and N the Neumann subdomain (defined under `Schwarz DN nonoverlap`); for an adjoint-paired impedance pair, D is the subdomain listed first |
| `displacement_jump_rms_<D>_<N>`, `velocity_jump_rms_<D>_<N>` | the same jumps as root mean square values over the interface, ‖q_D − Π_D q_N‖_{W_D} / √\|Γ_D\|, in m and m/s |

E2 is the Newmark discrete energy of the differentiated equation of motion
M ȧ + K v = ḟ (Prakash and Hjelmstad 2004, Eqs. (56)–(58) and (71)),
E2 = ½ aᵀ A a + ½ vᵀ K v with A = M + (Δt²/2)(2β − γ) K, evaluated as
½ aᵀ M a + (Δt²/2)(2β − γ) SE(a) + SE(v), where SE(x) is the strain energy
with the nodal field x in place of the displacement. This is ½ aᵀ M a + SE(v)
for Newmark with β = 1/4, γ = 1/2 (consistent mass) and
½ aᵀ M_L a − (Δt²/4) SE(a) + SE(v) for central difference (lumped mass). The
identity SE(x) = ½ xᵀ K x holds for the small-strain `linear elastic`
material, so E2 is written only for subdomains whose materials are all linear
elastic and whose integrator is Newmark without HHT-α or central difference,
and is `NaN` otherwise. On a linear problem without loads, continuity of the
interface velocity conserves E2 and continuity of the interface displacement
with β = 1/4, γ = 1/2 conserves E1; an undecomposed run conserves both.

## Schwarz coupling boundary conditions (subdomain files)

Inside each subdomain's `boundary conditions` block, coupling to a partner
subdomain is expressed with one of the Schwarz condition types below. All share
these keys:

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source` | yes | — | name of the coupled subdomain (its input-file basename) |
| `side set` | yes | — | this subdomain's coupling surface |
| `source side set` | for non-overlapping | `""` | the partner's coupling surface |
| `source block` | for overlapping | `""` | the partner element block searched for the overlap |
| `search tolerance` | no | `1.0e-6` | geometric search tolerance for locating partner points |

### `Schwarz overlap`

Dirichlet-overlap coupling: the subdomain reads its partner's displacement in
the shared overlap region.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source block` | yes | — | partner element block covering the overlap |
| `weak` | no | `false` | weak (integrated) rather than pointwise transfer |
| `compute overlap L2 relative error` | no | `""` | report the overlap L2 relative error of `disp`, `velo`, or `acce` (drives the overlap-error mesh-swap criterion) |

### `Schwarz DN nonoverlap`

Non-overlapping Dirichlet–Neumann coupling across a shared interface. The two
sides must take opposite roles. The Dirichlet side receives the partner's
projected interface *displacement* (not the projected current position): on a
curved interface the two sides discretize the geometry as different facet
polyhedra, so a position transfer would inject the chordal mismatch between
the two trace meshes — of the order of the coarser side's facet sagitta — as
a spurious scalloped interface displacement at the coarse-facet frequency. On
flat interfaces the two forms coincide, because the L2 projection reproduces
linear functions. Regression: test 123
(`schwarz-nonoverlap-static-inclusion-curved-interface.jl`, a stiffer
circular inclusion in a square matrix with a 36:20 non-conformal interface
discretization).

Note that for a stiff inclusion fully embedded in a softer matrix the DN
fixed-point map has gain greater than one, so *fixed* relaxation diverges for
every relaxation parameter; use Aitken relaxation
(`relaxation: aitken secant` or `aitken recursive` on the controller).

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner interface surface |
| `default BC type` | no | `Dirichlet` | this side's role: `Dirichlet` or `Neumann` (the two sides must be opposite) |
| `swap BC types` | no | `false` | swap the Dirichlet/Neumann roles between Schwarz iterations |
| `constrained` | no | `false` | constrained exchange: the Dirichlet side imposes one projected quantity and derives the other two kinematic fields from its own Newmark relations; the Neumann side receives the d'Alembert reaction of the Dirichlet side; must be set on both sides |
| `constraint` | no | `velocity` | quantity imposed by the constrained exchange: `velocity` or `displacement`; one value per pair, which may be given on either side or on both (then equal) |

**Interface residuals.** At every Schwarz iteration the log reports, for each
Dirichlet–Neumann pair with Dirichlet side D and Neumann side N, the one-sided
jumps ‖q_D − Π_D q_N‖_{W_D} / ‖q_D‖_{W_D} of the velocity and of the
displacement, with Π_D the Dirichlet projector and W_D the boundary mass
matrix of the Dirichlet interface, ‖x‖²_W = Σ_c x_cᵀ W x_c over the three
components, and the force residual ‖Π_Dᵀ r_D − f_N‖_{W_N⁻¹} / ‖Π_Dᵀ r_D‖_{W_N⁻¹},
with r_D = −(M a + f_int − f_body − f_boundary) on the interface rows of D and
f_N the interface force applied on N in the same iteration. For the
unconstrained exchange these values are reported only.

**Constrained exchange.** With `constrained: true` both transfer operators are
built from one cross mass matrix B, integrated over the facets of the side
with more interface nodes, Π_D = W_D⁻¹ B, and the Neumann side's force transfer
is Π_Dᵀ. The Neumann side receives Π_Dᵀ r_D, the d'Alembert reaction of the
Dirichlet side, so the two interface rows add to the row of the undecomposed
problem. The Dirichlet side computes from its state (u_n, v_n, a_n) at the
start of the step u_pre = u_n + Δt v_n + (½ − β) Δt² a_n and
v_pre = v_n + (1 − γ) Δt a_n (β = 0 for central difference) and imposes, with
`constraint: velocity`, v = Π_D v_N, a = (v − v_pre)/(γ Δt),
u = u_pre + β Δt² a, and with `constraint: displacement`, u = Π_D u_N,
a = (u − u_pre)/(β Δt²), v = v_pre + γ Δt a. Relaxation acts on the
constrained quantity only. A pair converges when the jump of its constrained
quantity and the force residual are both at or below `relative tolerance`, or
when the root mean square jump over the interface is at or below
`constraint absolute tolerance` (a top-level key in the units of the
constrained quantity, m/s for the velocity and m for the displacement
constraint; default 0, which leaves the relative test alone); when every
Schwarz coupling is constrained this replaces the displacement criterion,
otherwise both must hold. A residual that stays above the tolerance while the
displacement update has converged and decreased by less than 5% since the
previous such iteration is handled by `stalled interface jump action`.

The velocity constraint admits different time steps on the two sides, as in
Gravouil and Combescure (2001): the side with the finer step receives the
partner's velocity (Dirichlet side) or reaction (Neumann side) interpolated
linearly in time from the partner's substep history, the relaxation state is
kept per substep time, and the jump and force residual of the stopping rule
are the largest over the substeps of the stop. With different steps make the
side with the coarser step the Dirichlet side where the integrators allow it:
on the cantilever with a 4:1 step ratio this kept E2 within 4e-12 and E1
within 0.11% over 10 ms (explicit pair, and implicit coarse side with explicit
fine side), while with the finer side as the Dirichlet side E2 grew by up to
7e-8 at the 1e-12 tolerance and E1 varied by up to 47%. An explicit coarse
side with an implicit fine side has no good choice: as the Dirichlet side the
explicit member diverges, and as the Neumann side E1 varied by up to 17%. The
run aborts with
`swap BC types`, with HHT-α, with reduced order models, and with
`constraint: displacement` unless both sides are Newmark with the same
β > 0, γ, and time step.

The initial acceleration of a constrained pair is found by a Schwarz
iteration at t = 0: the Dirichlet side imposes the projected partner
acceleration, the Neumann side receives the d'Alembert reaction, and each side
recomputes its initial acceleration with the interface rows included. The
Dirichlet datum is relaxed with the fixed factor `relaxation parameter`, also
when `relaxation` names an Aitken method, and the iteration stops when the
acceleration jump and the force residual are at or below `relative
tolerance`. Without it the first step changes the conserved energy by up to 1%
on the cantilever, and between two central difference subdomains the
difference of the two interface accelerations alternates in sign at every
step.

The choice of the Dirichlet side decides whether the iteration converges. Its
linearized gain is G = −Π_D H_N Π_Dᵀ S_D, with H_N = C M̃_N⁻¹ Cᵀ the interface
flexibility of the Neumann side (M̃ = M + β Δt² K, lumped for central
difference) and S_D the interface dynamic stiffness of the Dirichlet side;
fixed relaxation θ converges when every eigenvalue g of G satisfies
|1 − θ + θ g| < 1. On the cantilever with equal steps of 1e-6 s:

- implicit and explicit subdomains: make the implicit subdomain the
  Dirichlet side. With the explicit side as Dirichlet the eigenvalues lie in
  [−13.3, −1.6] on conforming meshes and every θ ≥ 0.5 diverges; with the
  implicit side they lie in [−0.61, −0.08].
- two subdomains with the same integrator on conforming meshes: the
  eigenvalues are near −1, so θ = 0.5 converges in one to eight iterations and
  θ = 1 does not converge.
- nonconforming meshes: make the finer side the Dirichlet side.

### `Schwarz impedance nonoverlap`

Non-overlapping impedance (absorbing) coupling — the default and recommended
non-overlapping method. The interface transmits the partner traction plus a
dashpot term, making the interface energy exchange dissipative. See
`docs/notes/schwarz-coupling` for the theory.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner interface surface |
| `robin parameter` | no | `0.0` | Robin coefficient α; must be identical on both sides under `adjoint pairing` (the default), where it affects the convergence rate and not the converged solution when the interface jump closes (conforming or nested meshes); on nonconforming meshes the accepted iterate carries a residual jump and the effect of α on the converged solution has not been measured; per-side values are allowed with `adjoint pairing: false` |
| `impedance scale` | no | `1.0` | scalar scaling of the dashpot impedance; must be > 0 |
| `adjoint pairing` | no | `true` | use the adjoint-paired shared cross-mass transfer (recommended); `false` restores the legacy per-side transfer |

### `Schwarz RR nonoverlap`

The classical Robin-Robin coupling `traction + α·displacement = data`: the
Robin spring is the only coupling term and there is no dashpot (`impedance
scale` is rejected under this keyword). The condition is not absorbing, so in
elastodynamics it can pump energy at the interface (issue #176; on the
cantilever and nested-cylinders benchmarks of `docs/notes/schwarz-coupling`
nearly every dynamic Robin-Robin run ends by element inversion after
exponential energy growth, on overlap and nonoverlap decompositions alike); it is intended
for quasi-statics and for comparison against the classical Robin-Robin
literature, and the run warns when it is used with a dynamic time integrator.
For dynamics prefer `Schwarz impedance nonoverlap`.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner interface surface |
| `robin parameter` | yes | — | Robin coefficient α (positive); the two sides may use different values under `adjoint pairing: false` (the default) |
| `adjoint pairing` | no | `false` | `true` uses the adjoint-paired shared cross-mass transfer, which makes the Robin spring a conservative interface spring and requires one shared α per interface |

### `Schwarz impedance overlap`

Overlapping impedance coupling with recovered or consistent partner tractions.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source block` | yes | — | partner element block covering the overlap |
| `robin parameter` | no | `0.0` | Robin coefficient α |
| `impedance scale` | no | `1.0` | dashpot impedance scaling (scalar, or a P/S-split schedule) |
| `partner traction` | no | `auto` | `auto`, `consistent traction`, or `recovered stress` |
| `transfer` | no | `variational` | partner-field transfer: `pointwise` or `variational` |
| `transfer quadrature subdivisions` | no | `1` | quadrature refinement for variational transfer (integer ≥ 1) |
| `representable dashpot` | no | `false` | restrict the dashpot to the representable subspace |
| `content aware absorption` | no | `false` | content-aware absorption variant |

Requesting `Schwarz impedance overlap` forces consistent nodal stress recovery
on for the coupled model.

### `Schwarz contact`

Frictionless or tied contact enforced through Schwarz coupling.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner contact surface |
| `friction type` | yes | — | `frictionless` or `tied` |
| `swap BC types` | no | `false` | swap the interface roles between iterations |

## Mesh swapping (`swaps`)

Both single- and multi-domain runs can replace a mesh mid-run when a criterion
fires. `swaps` is a list; each entry names a `replacement` file and a
`criterion` (in a multidomain run it also names the `subsim` to replace). Swaps
cannot be combined with `restart`.

| Key | Required | Meaning |
|---|---|---|
| `subsim` | multidomain | subdomain to replace |
| `replacement` | yes | replacement input file |
| `criterion` | yes | trigger (see below) |

Criterion `type` values: `time` (with `t_swap`); `stress recovery` (with
`tolerance`, default `1.0e-2`, and `direction` `refine`/`coarsen`);
`elastic to plastic transition` (with `tolerance`, default `0.05`);
`overlap l2 relative error` (with `tolerance`, default `1.0e-6`, and `direction`
`refine`/`coarsen`).

## Restart (`restart`)

A run can resume from an Exodus snapshot by adding a `restart` block. Its one
key is `index`, the snapshot time-step index (negative values count back from
the end, so `-1` is the last snapshot). Restart is incompatible with `swaps`,
`initial conditions`, `j2 plasticity`, `Schwarz contact`, and mesh smoothing.

```yaml
restart:
  index: -1
```

## Canonical examples

- Overlapping Schwarz: `examples/overlap/`
- Overlapping Schwarz with J2 plasticity (circular laser weld, graded pair):
  `examples/overlap/static-same-step/clw/`
- Non-overlapping impedance (same and subcycled steps):
  `examples/nonoverlap/dynamic-same-step/`,
  `examples/nonoverlap/dynamic-different-steps/`
- Contact: `examples/contact/`
- Adaptive mesh swapping: `examples/adaptive-time-stepping/`,
  `examples/ahead/`

See `docs/notes/schwarz-coupling` for the theory and stability analysis of the
impedance coupling.
