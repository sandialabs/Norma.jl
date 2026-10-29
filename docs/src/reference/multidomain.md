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
| `absolute tolerance` | unless every coupling is a matching pair | — | absolute tolerance of the update criterion, on the change, between Schwarz iterations, of u + Δt v over all degrees of freedom of all subdomains (Δt the controller step); a stop converges when this or the relative test holds. When every Schwarz coupling is a matching pair it is not used, and it is ignored with a warning if given |
| `relative tolerance` | unless every coupling is a matching pair | — | relative tolerance of the update criterion, on the same change divided by the norm of u + Δt v over all degrees of freedom of all subdomains, so the test does not depend on the coordinate origin. Matching pairs use the matching criterion instead (see `matching relative tolerance`); when every Schwarz coupling is a matching pair it is not used, and it is ignored with a warning if given |
| `matching relative tolerance` | when a pair sets `matching` | — | relative tolerance of the matching criterion, on the one-sided jump of the matched quantity and on the force residual (see `Schwarz DN nonoverlap`); the controller's `relative tolerance` does not apply to matching pairs |
| `matching absolute tolerance` | no | `0` | absolute tolerance of the matching criterion on the root mean square interface jump, in the units of the matched quantity (velocity or displacement); the controller's `absolute tolerance` does not apply to matching pairs |
| `unconverged step action` | no | `abort` | what to do with a step that exhausts `maximum iterations` without converging: `abort` the run, or `warn` and continue with the last iterate |
| `stalled interface jump action` | no | `warn` | what to do when the interface residual of a matching Dirichlet–Neumann pair stops decreasing above `matching relative tolerance`: `warn` and accept the iterate, or `abort` |

A step that reaches `maximum iterations` without converging reports the
measures it stopped at against the tolerances that apply: the change of
u + Δt v (the update criterion) against `absolute tolerance` and `relative tolerance`, the interface
residual of the matching pairs against `matching relative tolerance`, or both.
Under the default `abort` the run stops; `warn` continues with the last
iterate.

### Linear solver of the subdomains

The time of a Schwarz run is that of its subdomain solves: on one thread the
coupling layer itself (the transfer operators, the partner histories, the
relaxation) took 1.5% of the matrix-inclusion split and 7.6% of the
nonconforming Dirichlet-Neumann cantilever, and the construction of the
projectors 23% of the short inclusion run. Subdomains are often small and
solved many times, where the conjugate gradient repeats many iterations on a
small matrix and a factorization takes a few milliseconds. Measured on one
thread of a Ryzen 9 9900X (October 24, 2026), wall seconds of the whole run
by the `linear solver` and `preconditioner` of the subdomain solver blocks
([solvers](solvers.md)): the nonconforming cantilever (315 and 66 nodes,
1992 solves) 2.45 with `cg`, 2.02 with `jacobi`, 0.97 with `cholesky`, 1.72
with `direct`; the matrix-inclusion split (703 and 262 nodes, 48 solves)
1.19, 1.22, 1.02, and 0.81. For small dynamic subdomains the reused Cholesky
factor as preconditioner is the fastest option measured; for small
quasi-static ones, the direct solver.

On a larger coupled case (`benchmark/coupled.yaml`: a clamped block of 24³
elements and a free block of 16³ elements under Newmark with the matching
Dirichlet-Neumann exchange and Anderson acceleration, 15 Schwarz iterations
per step, 600 subdomain solves over 20 steps) the coupling layer took 5% of
the run on one thread and 12% on 12 threads, the construction of the
projectors of the nonconforming interface included. The two subdomain solves
of each Schwarz iteration are the rest, and the run is 7 times that of the
monolithic bar, in proportion to the number of solves. The coupling layer
runs on one thread. Its largest items per Schwarz iteration were the
d'Alembert reaction of the Dirichlet side, computed on all degrees of
freedom for the Neumann datum and again for the criterion, and the momentum
and the row sums of the mass that the criterion recomputed; since October
24, 2026 they are computed on the interface rows only, and the criterion
interpolates the partner histories on those rows, which took the coupled run
on 12 threads from 21.2 to 17.9 s (`benchmark/README.md`).

### Inexact Newton control of the subdomain solves

The optional block `inexact Newton` sets the relative Newton tolerance of every
implicit subdomain (`Hessian minimizer`) before the subdomain solves of each
Schwarz iteration, replacing the `relative tolerance` of the subdomain
`solver` blocks. At Schwarz iteration k ≥ 1 of a stop the tolerance is
clamp(η r_{k−1}, floor, ceiling), where r_{k−1} is the relative Schwarz
residual of the previous Schwarz iteration of that stop: the larger of the
interface jump and the force residual for matching Dirichlet–Neumann pairs,
the relative change of u + Δt v otherwise. At the first Schwarz iteration of a
stop (k = 0) it is the ceiling. The tolerance and the number of Newton
iterations of every Schwarz iteration are logged with the stopping criterion.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `tolerance ratio` | no | `0.01` | ratio η > 0 of the relative Newton tolerance to the previous relative Schwarz residual; a smaller η solves each subdomain to a smaller Newton tolerance relative to the Schwarz residual. The key `forcing factor`, which set η before October 2026, aborts |
| `minimum tolerance` | no | `1.0e-12` | floor of the relative Newton tolerance |
| `maximum tolerance` | no | `1.0e-2` | ceiling of the relative Newton tolerance, used at the first Schwarz iteration of a stop |
| `tight final solve` | no | `false` | accept Schwarz convergence only at a Schwarz iteration whose Newton tolerance was the floor; a convergence reached above the floor is rejected and the Schwarz iteration continues with the tolerance at the floor for the rest of the stop |

The defaults come from the inexact Schwarz study of October 2026 on the
neohookean cantilever of `studies/dynamic-stability` (two implicit subdomains,
matching exchange, matching relative tolerance 1e-12, 10 000 stops): with
η = 0.01 and no tight final solve the run needed 0.735 of the Newton
iterations of the run with the constant Newton tolerance 1e-12, at a relative
displacement error of 1.4e-12 against that run at the final time, whereas the
constant Newton tolerance 1e-6 needed 0.667 at an error of 6.0e-9 (0.717 and
1.8e-12 for η = 0.01 on the mesh with 8 times the elements and 20 000 stops).
Every Schwarz iteration restarts the subdomain solves from the state at the
start of the stop, not from the previous Schwarz iterate, so a constant Newton
tolerance larger than the Schwarz tolerance leaves an error of the size of the
Newton error, which neither the update criterion nor the matching criterion
detects.

### Relaxation and acceleration

| Key | Required | Default | Meaning |
|---|---|---|---|
| `relaxation` | no | fixed | `aitken recursive` (Irons–Tuck) or `aitken secant` adaptive relaxation, or `anderson` (Anderson acceleration of the transmitted datum, for every relaxed coupling); omit for a fixed factor |
| `anderson depth` | no | `10` | number m of previous iterates combined by Anderson acceleration; it must cover the Schwarz iterations a stop needs (5 is too short with an explicit Dirichlet side) |
| `relaxation parameter` | no | `1.0`; `0.5` under `anderson` | relaxation factor θ: the constant factor under fixed relaxation; under either Aitken method the factor used wherever an adaptive one is not yet available; under `anderson` the mixing parameter β of the update (see below), which is also the factor of the first Schwarz iteration of every stop |
| `aitken N0 parameter` | no | `1` | Schwarz iteration, counting from zero, at which the adaptive factor takes over from `relaxation parameter` (Aitken methods only) |
| `aitken scope` | no | `per interface` | `per interface`: each relaxed interface has its own Aitken factor; `global`: one factor for all relaxed interfaces, from their stacked residuals (Aitken methods only; see below) |
| `interface predictor` | no | `false` | extrapolate the interface state at the start of each step |
| `naive stabilized` | no | `false` | naive interface stabilization |

Relaxation applies to the datum that a nonoverlap or contact coupling
transmits: for `Schwarz DN nonoverlap` and `Schwarz contact` it is the
kinematic datum of the Dirichlet side (the interface displacement, or for a
matching pair the matched quantity); for `Schwarz RR nonoverlap` it is
the Robin datum of the subdomain listed first in `domains`, and only that
datum: other loads on the interface nodes, such as Neumann or pressure
conditions, are applied unrelaxed. The subdomain listed first is solved before
its partner in each Schwarz iteration, so its datum comes from the partner's
previous iteration and carries the state between iterations; Anderson
acceleration converges only when these lagged data are the accelerated ones
(issue #245). `Schwarz overlap matching` relaxes the
constrained datum of each side with its own history. `Schwarz overlap` is not
relaxed: a simulation whose couplings are all `Schwarz overlap` conditions
aborts at setup when `relaxation` or a `relaxation parameter` other than one
is set, since the keys would have no effect. Both Aitken forms work with every relaxed
coupling, and `relaxation parameter` and `aitken N0 parameter` mean the same
thing in each.

`relaxation parameter` is not ignored when `relaxation` names an Aitken method.
It is the factor applied for Schwarz iterations below `aitken N0 parameter`,
and the factor used whenever an adaptive factor cannot be formed, which
includes the first iterations of every step before two iterates exist to
compare.

An adaptive factor cannot be formed either when the difference of two
successive interface iterates, or of two successive interface residuals, is
not resolved: its norm is at most 1e-14 times the norm of the current
interface iterate, which is 45 times the unit roundoff of double precision
arithmetic (2.2e-16). The secant
form then applies `relaxation parameter` and the recursive form keeps its
previous factor. Before October 2026 this floor was absolute (1e-20 on the
squared norm of the difference). On
`examples/nonoverlap/static-same-step/matrix-inclusion` at 5% strain in 10
load steps, with `relaxation parameter: 0.5` and a relative tolerance of
1e-12, the absolute floor was reached at interface differences of about
1e-10 m, the factor 0.5 used in its place made the iteration diverge by a factor of
1.33 per Schwarz iteration, and every load step ended at the limit of 128
Schwarz iterations. With the relative floor the same run converges in 26 to
33 Schwarz iterations per load step.

By default each relaxed interface (one side set of one coupling) has its own
Aitken factor, formed from its own residuals. With `aitken scope: global`
the residuals and iterates of all relaxed interfaces are stacked into one
vector and a single factor is formed from it, by the same secant or recursive
formula, and applied to every interface:
θ = −(Σᵢ dᵢ·δᵢ) / Σᵢ ‖δᵢ‖² for the secant form and
θ = −θ_prev (Σᵢ rᵢ⁽ⁿ⁻¹⁾·δᵢ) / Σᵢ ‖δᵢ‖² for the recursive form, with the sums
over the interfaces. The residual of an interface depends on the data of the
others, which per-interface factors ignore. Each interface contributes its
relaxed datum on the interface: the Robin datum, the projected trace of the
constrained datum of a matching pair, or the projected trace of the partner
displacement of an unconstrained Dirichlet-Neumann or contact pair. For one
relaxed Robin or matching interface the two scopes are identical. The
per-interface factor of an unconstrained pair is formed from the whole
partner displacement field instead of its trace, so all side sets coupled to
the same partner subdomain already share one factor under `per interface`.
The residuals are gathered when the first relaxed
interface of a Schwarz iteration applies its datum, that is, at the end of
the previous Gauss-Seidel pass, and each interface applies the factor to the
partner data it reads at its own application. Interfaces are grouped by the
relaxed quantity (Robin datum, displacement, velocity), with one factor per
group, so that no sum mixes units. The coupled initial acceleration of
matching pairs keeps per-pair factors. Where the global factor does not
apply, the run aborts rather than fall back to other factors: at setup when
`aitken scope: global` is set without an Aitken `relaxation`, with a
`Schwarz overlap matching` coupling, or with a Dirichlet-Neumann or contact
coupling whose partner is a reduced-order model; and at the first stop with
several substeps of a relaxed side (see below), which depends on the time
steps at run time.

Aitken acceleration is applied only to stops with a single substep, that is
when the controller `time step` equals the relaxed subdomain's time step. A
windowed stop couples all of its time slots in one Schwarz iteration, where the adaptive
factors were measured to diverge or to need more iterations than a fixed
factor, so such stops use `relaxation parameter` throughout. This is automatic and needs no input.

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
| `stored_energy`, `kinetic_energy`, `total_energy` | blended energy E: strain energy, kinetic energy (the staggered form ½ vᵀ M_L v − (Δt²/8) aᵀ M_L a for central difference subdomains), and their sum |
| `rate_invariant_total` | Ẽ, the rate invariant: the sum over subdomains of the energy functional evaluated on the velocity and the acceleration, which is the invariant of the time-differentiated equation of motion; in units of energy per time squared |
| `rate_invariant_<subdomain>` | Ẽ of one subdomain, one column per subdomain in the order of `domains` |
| `staggered_kinetic_interface`, `staggered_kinetic_interior` | staggered kinetic energy of the central difference subdomains, summed over the rows of the nodes in any Schwarz side set and over the remaining rows; `NaN` without central difference subdomains |
| `displacement_jump_<D>_<N>`, `velocity_jump_<D>_<N>` | one-sided interface jumps of each Dirichlet–Neumann pair at the end of the stop, D the Dirichlet and N the Neumann subdomain (defined under `Schwarz DN nonoverlap`) |
| `displacement_jump_rms_<D>_<N>`, `velocity_jump_rms_<D>_<N>` | the same jumps as root mean square values over the interface, ‖q_D − Π_D q_N‖_{W_D} / √\|Γ_D\|, in units of displacement and of velocity |
| `impulse_residual_x_<D>_<N>`, `_y_`, `_z_`, `impulse_residual_relative_<D>_<N>` | matching pairs: the trapezoid-in-time sum over the Neumann side's substeps of the force it received, minus Π_Dᵀ times the trapezoid sum of the Dirichlet side's reaction over its own substeps, summed over the interface nodes per component (in units of force times time), and its W_N⁻¹ norm relative to that of the transferred impulse; `NaN` for direct pairs |

For a central difference subdomain, `kinetic_energy` is the product
½ v_{n−½}ᵀ M_L v_{n+½} of the two half-step velocities around the stop n,
with M_L the lumped mass matrix, which equals ½ vᵀ M_L v − (Δt²/8) aᵀ M_L a
for the velocity v and acceleration a at the stop. With this kinetic energy
central difference conserves E exactly on a linear problem without loads. The
product is not a quadratic form of one velocity and can be negative. For a
body released at rest it is −(Δt²/8) aᵀ M_L a at t = 0, so E(0) is below
the initial stored energy. Measured on the linear elastic cantilever of
`studies/dynamic-stability` (case `beam_mono_E_L1`, Δt = 1 µs, 10 ms): the
stored energy at t = 0 is 1506.28 J, the kinetic energy −10.68 J, and
|E(t)/E(0) − 1| at most 1.0e-14. The kinetic energy ½ vᵀ M_L v of the
velocity at the stop gives E(0) equal to the stored energy, but E then varies
by 0.62% over the same run, so it is not written.

Ẽ is the rate invariant, the Newmark discrete invariant of the differentiated equation of motion
M ȧ + K v = ḟ (Prakash and Hjelmstad 2004, Eqs. (56)–(58) and (71)),
Ẽ = ½ aᵀ A a + ½ vᵀ K v with A = M + (Δt²/2)(2β − γ) K, evaluated as
½ aᵀ M a + (Δt²/2)(2β − γ) SE(a) + SE(v), where SE(x) is the strain energy
with the nodal field x in place of the displacement. This is ½ aᵀ M a + SE(v)
for Newmark with β = 1/4, γ = 1/2 (consistent mass) and
½ aᵀ M_L a − (Δt²/4) SE(a) + SE(v) for central difference (lumped mass). The
identity SE(x) = ½ xᵀ K x holds for the small-strain `linear elastic`
material, so Ẽ is written only for subdomains whose materials are all linear
elastic and whose integrator is Newmark without HHT-α or central difference,
and is `NaN` otherwise. On a linear problem without loads, continuity of the
interface velocity conserves Ẽ and continuity of the interface displacement
with β = 1/4, γ = 1/2 conserves E; an undecomposed run conserves both.

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
| `search tolerance` | no | `1.0e-6` | tolerance, in the parametric coordinates of a partner element, for accepting a point as inside it (dimensionless) |

### `Schwarz overlap`

Dirichlet-overlap coupling: the subdomain reads its partner's displacement in
the shared overlap region.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source block` | yes | — | partner element block covering the overlap |
| `weak` | no | `false` | weak (integrated) rather than pointwise transfer |
| `compute overlap L2 relative error` | no | `""` | report the overlap L2 relative error of `disp`, `velo`, or `acce` (drives the overlap-error mesh-swap criterion) |

With `weak: true` the datum on the subdomain's Schwarz boundary Γ is the
L2 projection W⁻¹L of the partner's field, with W the facet mass matrix of
Γ and L the cross-mass matrix of Γ against the partner's volume shape
functions. Since 2026-10-19, L and W are integrated by the segmented rule of
the matching exchange (below): in the plane of each facet the partner's
elements appear as the polygons in which the plane cuts them (the convex
hull of the crossings of the element's corner edges with the plane, and a
face lying on the plane is counted once, from the side of the facet normal),
the facet is clipped by each, and the pieces take the seven-point rule of
degree five with each point located in the element by the inverse
isoparametric map. On affine elements (tetrahedra, parallelepiped
hexahedra) and triangular or parallelogram facets this is exact, the
projection reproduces constants and linear fields, and its transpose
carries a uniform traction exactly: on the 1:0.75 beam overlap the
projector agrees with an independent cut integration to 5e-15, and the
static patch test of the two-sided surface coupling (gist 6478a5bf,
`patch.jl`, C2b) gives a strain error of 1e-12 at both refinement levels,
against 4.4e-2 with the facet's own rule (test
`schwarz-overlap-dynamic-cantilever-weak.jl`; results in
~/Projects/Norma-Overlap-Energy/patch-norma-2026-10-19). A facet the
sections do not cover, or a curved one, takes the facet's Gauss rule with
point location, and the setup log reports the count.

### `Schwarz overlap matching`

Energy-conserving two-sided overlap coupling (since 2026-10-19). Each
subdomain is constrained on its own Schwarz boundary Γ_k to the L2
projection W_k⁻¹L_k of the partner's matched quantity, with L_k the
cross-mass matrix of Γ_k against the partner's volume shape functions by the
segmented rule above, and derives the other two kinematic fields from its
own Newmark relations as the Dirichlet side of a matching pair does. It
receives on all of its nodes the partner's d'Alembert reaction on the
partner's Schwarz boundary through the transpose of the partner's projector,
−Π_jᵀρ_j. Mass, internal force, stiffness, and energy of every quadrature
point inside the partner are weighted by one half in the assembly (constant
weights; a point inside n partners has the weight 1/(1 + n)), and so are the
`Neumann` and `Neumann pressure` loads at their facet integration points, so
that a load on a boundary inside the partner is counted once; the spring term
of a `Robin` condition is not weighted. With these weights the
fixed point of the Schwarz iteration is the blended constrained problem, the
window balance of the rate invariant cancels as for a nonoverlap matching
pair, and the blended energy output and the rate invariant count the overlap
once with the same weights. Both sides are constrained, so both take one step
per window: equal time steps are required, and the coupled initial
acceleration, the matching criterion (`matching relative tolerance`), and
Anderson acceleration apply as for a matching pair, with the Anderson step
of each side taken at its own application on the partner data of the
current iteration. Both sides of the overlap must set the condition.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source block` | yes | — | partner element block that contains this side set |
| `matching` | no | `velocity` | matched quantity, `velocity` or `displacement` |
| `matching at` | no | `every substep` | `window end` on the fine-step side of a pair with different steps (below) |
| `interface solve` | no | `iterative` | `direct` on both sides: two central difference sides at equal steps solve both reactions from one symmetric block system per step, [A₁ −B; −Bᵀ A₂] with A_k = m_Γk⁻¹ + Π_k m_j⁻¹ Π_kᵀ and B the coupling of the rows the two transfers share, factored once; no Schwarz iteration, the velocity constraint on both boundaries to roundoff (test 140) |
| `search tolerance` | no | `1.0e-6` | as for `Schwarz overlap` |

**Different time steps (`matching at: window end`).** Both sides of the
overlap matching coupling are constrained, so with different steps one
constraint sits on the subcycled side. Constrained at every substep to the
partner's interpolated velocity, that side is over-constrained in time, its
substep reactions oscillate, and no single force on the coarse side is
consistent with both the rate invariant (which pairs force increments) and
the impulse (which pairs force averages): this is the arrangement the
nonoverlap role rule excludes (rate invariant conserved, impulse residual of
the order of the impulse, energy varying by its order). With `matching at:
window end` the subcycled side imposes the matched velocity at the end of the
window only; at its interior substeps its interface rows are free and carry
the linear ramp in time between its reaction at the window start (the
converged one of the previous stop) and at the window end (its last pass).
The force increments of its substeps then add up to the change of its
reaction and the window balance cancels, and the impulse of the ramp equals
the impulse the coarse side receives: both conserved, as with the coarse-step
side constrained. The interior substeps respond to the ramp and the end
substep undoes that response, a loop with a gain of about −(m − 1) for m
substeps that the relaxation of the datum does not pass through; under
Anderson acceleration the window-end reaction is a second field of the
interface's step (joint step on datum and reaction, the reaction residual
scaled once per interface to the datum residual), and the iteration converges
at two to four times the count of the coarse-step-constrained arrangement.
Fixed relaxation damps the reaction with the same factor: 0.5 and 0.25 diverge, and 0.1 reaches 4e-12 but not 1e-12 within 256 iterations on the first stop, so Anderson acceleration is the method for this constraint. The matching criterion measures the jump at the window end only and adds to the force residual the difference between the reaction at the window end and the end value of the ramp that the interior substeps applied, in the W⁻¹ norm of the constrained side relative to the larger of the reaction and the momentum rate of its interface (issue #244: in the first windows of the cylinders the reaction is at roundoff, and relative to it alone the residual stayed near 1e-2 after the iteration had converged): the window balance of the rate invariant contains that difference, so a window is converged only when it is at the tolerance. The joint residual is formed from the end value of the ramp applied in the pass, which is kept until the joint step sets the next one, and every window starts from the measured reaction of the previous one. On the beam (October 19, 2026, after this was fixed) the window-end constraint conserves Ẽ to between 2.3e-13 and 7.7e-13 at 1.2 to 1.3 times the iterations of the coarse-step-constrained arrangement; at mesh ratio 1:0.5 with the fine-mesh, fine-step side constrained both pairs converge with mixing 0.5 and 0.25; the explicit pairs stall in no window and the implicit pairs in 6 to 17% of them, near 2e-12. Before the fix the measured reaction overwrote the ramp end value, the joint residual of the reaction was zero, and the same cases needed up to 40 iterations per window, stalled in most windows, and failed at 1:0.5 with mixing 0.5. The
same key applies to the Dirichlet side of a nonoverlap matching pair, which
lifts role rule 1 (coarse-step side as the Dirichlet side) from a requirement
to a preference for the iteration count. Measured on the beam, free part
subcycled 4:1 as the constrained side, 1 ms, Anderson 0.5, 1e-12
(~/Projects/Norma-Overlap-Energy/window-end-2026-10-19): nonoverlap EE
13.5/17 iterations, δẼ 9e-14, impulse residual 2e-16, E within 1.1e-3; II
24.8/41, δẼ 5e-12, impulse 2e-16, E within 4.0e-3; overlap EE 13.5/16
(conforming) and 17.1/23 (1:0.75), δẼ ≤ 1e-11; overlap II 21.4/33 and
23.5/36, δẼ 1.4e-10. Tests 136 and 140.

Measured on the beam overlaps (example `cantilever-matching`, 1:0.75, and
the conforming meshes) over 50 steps of 1 µs with Anderson acceleration at
tolerance 1e-12: E and Ẽ conserved to 1e-11 for the trapezoidal and the
central difference pair (test 140, `schwarz-overlap-matching.jl`); the
transfer is exact (section above) and the two-sided static patch test passes
to roundoff. The force transferred to a side is applied to all of its nodes
and is part of its external forces; the reaction on Γ_k is the residual of
the full equation including it.

### `Schwarz DN nonoverlap`

Non-overlapping Dirichlet–Neumann coupling across a shared interface. The two
sides must take opposite roles. The Dirichlet side receives the partner's
projected interface *displacement* (not the projected current position): on a
curved interface the two sides discretize the geometry as different facet
polyhedra, so a position transfer would add the chordal mismatch between
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

In quasi-statics a subdomain with only Neumann data on its interfaces and no
support of its own has a singular stiffness, so each such subdomain must be
the Dirichlet side of at least one interface. When that forces the Dirichlet
role onto the stiffer side of an interface with a large stiffness ratio, the
iteration can diverge under every relaxation. On the quarter of the three
nested cylinders (`examples/jmp/concentric-cylinders/quasistatic`,
test 154, issue #246) the chains of two and three subdomains, in which the
core stays inside one subdomain, converge in 21 to 25 iterations per step with
Anderson acceleration; the chain of four, in which the core (shear modulus
37.6 GPa) is the Dirichlet side against the filler (72 MPa), inverts elements
or does not converge under Anderson, fixed, and Aitken relaxation, at any load
step, and also with Robin–Robin coupling on that interface with Robin
parameters E/h. Fixed relaxation converges only below about 0.005 (0.004
needs 4017 iterations for the first step, 0.008 diverges), which measures a
largest gain of about −400. The input `cylinders-3.yaml` of that example
records the measurements. Robin–Robin coupling on all three interfaces of the
chain of four converges (`cylinders-4.yaml`): the soft filler side of every
interface has a Robin parameter close to a Dirichlet condition and the stiff
side one near the filler's interface stiffness, so that the Robin term
supports the core, and Anderson acceleration needs 30 to 34 iterations per
step against 21 for the chain of three.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner interface surface |
| `default BC type` | no | `Dirichlet` | this side's role: `Dirichlet` or `Neumann` (the two sides must be opposite) |
| `swap BC types` | no | `false` | swap the Dirichlet/Neumann roles once per controller stop, before its Schwarz iteration |
| `matching` | no | — | selects the matching exchange and its matched quantity, `velocity` or `displacement`: the Dirichlet side imposes the projected matched quantity of the Neumann side and derives the other two kinematic fields from its own Newmark relations, and the Neumann side receives the d'Alembert reaction of the Dirichlet side; must be set on both sides with the same value, and requires the top-level `matching relative tolerance`. Without it the pair uses the standard exchange, which imposes the projected displacement, velocity, and acceleration together. The keys `constrained` and `constraint`, which selected this exchange before October 2026, abort |
| `matching at` | no | `every substep` | `window end`: a subcycled Dirichlet side matches at the window end only and carries the ramp of its own reaction at interior substeps (see `Schwarz overlap matching`) |
| `interface solve` | no | `iterative` | `direct` replaces the Schwarz iteration by one solve of the interface force per stop; only for a pair with `matching: velocity` of two central difference subdomains at equal time steps, set on both sides |

A matching pair sets `matching` on both sides and the controller file sets
`matching relative tolerance`, as in
`examples/nonoverlap/dynamic-same-step/cantilever-dn-nonconforming/`:

```yaml
# cantilever-free.yaml (Dirichlet side)
boundary conditions:
  Schwarz DN nonoverlap:
    - side set: ssx-
      source: cantilever-clamped
      source side set: ssx+
      default BC type: Dirichlet
      matching: velocity

# cantilever-clamped.yaml (Neumann side): the same block with
# side set ssx+, source cantilever-free, source side set ssx-,
# default BC type Neumann, and matching: velocity

# cantilever-multi.yaml (controller): the only coupling is a matching pair,
# so the controller has no `relative tolerance` or `absolute tolerance`
matching relative tolerance: 1.0e-08 # interface jump and force residual
relaxation: anderson
```

**Interface residuals.** At every Schwarz iteration the log reports, for each
Dirichlet–Neumann pair with Dirichlet side D and Neumann side N, the one-sided
jumps ‖q_D − Π_D q_N‖_{W_D} / ‖q_D‖_{W_D} of the velocity and of the
displacement, with Π_D the Dirichlet projector and W_D the boundary mass
matrix of the Dirichlet interface, ‖x‖²_W = Σ_c x_cᵀ W x_c over the three
components, and the force residual
‖Π_Dᵀ r_D − f_N‖_{W_N⁻¹} / max(‖Π_Dᵀ r_D‖_{W_N⁻¹}, ‖C_N M_N v_N‖_{W_N⁻¹}/Δt_N),
with r_D = −(M a + f_int − f_body − f_boundary) on the interface rows of D,
f_N the interface force applied on N in the same iteration, and C_N M_N v_N
the momentum of N on its interface rows. The momentum rate in the
denominator, the force that stops the interface mass in one step, is the scale
of the interface forces when these are at roundoff, as in a body that still
moves rigidly: there the relative change of a lagged reaction (N solved
before D) is a ratio of roundoff values and not a measure of convergence. The
window-end ramp residual (below) and the force residual of
`Schwarz overlap matching`, in the dual norm of the receiving side's lumped
mass on the rows where the transferred force acts, use the same scale. For the
standard exchange these values are reported only.

**Matching exchange.** With `matching` both transfer operators are
built from one cross mass matrix B, Π_D = W_D⁻¹ B, and the Neumann side's
force transfer is Π_Dᵀ. B = ∫ N_D N_Nᵀ dΓ is integrated by
segment-based quadrature: each facet of the integration side is split into
its intersections with the facets of the other side, in the plane of the
facet, and each piece is integrated with a seven-point triangle rule of
degree five. The facets may be TRI3, QUAD4, or straight-sided TRI6 (the
facets of TETRA10) on either side, in any combination: the shape functions
of a triangle are polynomials of degree one or two in the plane for any
triangle shape, those of a parallelogram QUAD4 of degree two, so the product
of two shape functions is of degree at most four on such facets. On
interfaces of planar triangles and parallelograms B is therefore exact for
any mesh ratio, both integration sides give the same matrix, and Π_D
reproduces linear fields (measured on the 1:0.75 beam with the free part
remeshed in TET4 and in TET10 against the HEX8 clamped part, both roles:
constants and linear fields to 5e-15, uniform traction transferred to
7e-16; October 18, 2026, data in
~/Projects/Norma-Schwarz-Coupling/element-topology-2026-10). W_D is integrated facet by facet with the same
quadrature as B, so that the row sums of the two agree exactly on every
facet, also where the physical-coordinate rule is inexact (on the paved cap
facets of the nested cylinders the two quadratures of W_D differ by up to
1e-2; on the beam they agree to roundoff); Π_D therefore reproduces
constants and the force transfer conserves the total force exactly on
every interface
(before 2026-10-18, with the facet's own Gauss rule on the facets of D,
Π_D reproduced linear fields but the entries of B carried a quadrature
error of 9e-3 to 1.1e-2 on the non-nested 1:0.75 beam interface, so the
transposed force transfer distributed a uniform traction over the Neumann
nodes with a strain error of 4e-3 that rose to 7e-3 under one refinement;
with the segmented rule both uniform-strain patch tests are reproduced to
roundoff, measured statically on rigel on October 18, 2026). A facet that
the other side's facets do not cover to a relative 1e-8 of its area (gap,
partial overlap, or a curved interface on which the projection of a
quadrature point fails), or a TRI6 facet whose midside nodes depart from
the edge midpoints by more than 1e-8 of the longest edge (a curved facet),
is integrated with its own Gauss rule and the
containing-facet projection of each point, and the log warns with the count
of such facets. B is then taken on the Neumann side when Π_D reproduces
constants to 1e-12 with it, and otherwise on the Dirichlet side, where its
row sums equal those of W_D exactly for any quadrature, so Π_D reproduces
constants to roundoff and the total force transferred to the Neumann side
equals the total reaction of the Dirichlet side. The log reports the side
used, the partition-of-unity error of Π_D, and, for information, that of the
unused projection onto the Neumann side. The Neumann side receives Π_Dᵀ r_D, the d'Alembert reaction of the
Dirichlet side, so the two interface rows add to the row of the undecomposed
problem. The Dirichlet side computes from its state (u_n, v_n, a_n) at the
start of the step u_pre = u_n + Δt v_n + (½ − β) Δt² a_n and
v_pre = v_n + (1 − γ) Δt a_n (β = 0 for central difference) and imposes, with
`matching: velocity`, v = Π_D v_N, a = (v − v_pre)/(γ Δt),
u = u_pre + β Δt² a, and with `matching: displacement`, u = Π_D u_N,
a = (u − u_pre)/(β Δt²), v = v_pre + γ Δt a. Relaxation acts on the
matched quantity only, and Aitken factors and Anderson coefficients are
formed from its projected interface trace Π_D q_N, not from the partner's
whole field. Anderson acceleration (`relaxation: anderson`, with the defaults
`anderson depth: 10` and `relaxation parameter: 0.5`) is the recommended
setting for matching pairs. With x_k the datum the Dirichlet side used at
Schwarz iteration k, f_k = G(x_k) − x_k the residual (G the partner datum the
iteration returns), t_k the projected interface trace of f_k, and ΔX, ΔF, ΔT
the differences of the last m + 1 iterates, it computes the coefficients
γ = argmin ‖t_k − ΔT γ‖ and the next datum
x_{k+1} = x_k + β f_k − (ΔX + β ΔF) γ, with β the mixing parameter
`relaxation parameter`. The first iteration of a stop has no history, so it is
the fixed relaxation x_1 = x_0 + β f_0. With β = 1 and a depth at least the
iteration count the iterates are those of GMRES on a linear interface map
(test 136 checks this); β < 1 damps the residual term. On the cantilever
(1 ms, tolerance 1e-12) β = 0.5 needed 7.0, 23.1, and 7.9 Schwarz iterations
per stop on the implicit pair, the mixed pair with the explicit side as the
Dirichlet side, and the subcycled 1:0.5 explicit pair with the coarse mesh as
the Dirichlet side, against 8.9, 28.7, and 8.8 at β = 1, with the same
energies; hence the default 0.5. The least-squares problem is solved by QR
with the oldest column dropped while the condition number exceeds 1e10. One
least-squares problem is solved for all relaxed interfaces of the simulation,
with their traces stacked into one vector, and its coefficients update every
interface; interfaces are grouped by the relaxed quantity (interface
displacement, interface velocity, or Robin datum), one problem per group. The
residual of one interface depends on the data of the others, so coefficients
computed for each interface separately fit wrong secant models: on a disk in a
matrix whose interface is coupled as four side sets
(`examples/nonoverlap/static-same-step/matrix-inclusion-split`, test 139) they
inverted an element or did not converge in 500 iterations, against 23 for the
stacked problem and 21 for the interface as one side set (issue #233). Groups
are not coupled: on the chain of four nested cylinders with Robin–Robin
coupling on one interface and Dirichlet–Neumann on the two others, which share
a subdomain with it, Anderson (depth 10) did not converge in 256 iterations for
any Robin parameter while Aitken did (163 to 171 iterations per step); with
Robin–Robin on all three interfaces, one group, it converged in 31 to 32. The
histories hold at most m + 1 columns. A column holds, for each relaxed field
of an interface, the iterate and the residual on the partner degrees of
freedom that the transfer reads: the nodes of the partner side set, or for an
overlap matching interface the nodes in the support of its projector. For a
Robin interface the datum is a field of its own side, nonzero on the nodes of
its side set only, and the step is computed on those nodes. A column holds the
residual trace once per interface. Under a static partner, whose
velocity and acceleration are constant, the displacement is the only relaxed
field of an unconstrained pair. On the inclusion split a column holds 114
footprint entries twice and 258 trace entries, 3.9 KB, against 401.5 KB when
the histories held the whole partner field of all three kinematic fields. The
memory of the histories is logged at every Schwarz iteration, and its peak at
the end of every stop. The
stacked step is computed when the first relaxed interface of an iteration
applies its datum, from the partner data at that moment, so an interface whose
partner is solved after that moment and before its own subdomain receives the
partner data of the previous iteration. A Neumann side solved after the first
relaxed Dirichlet side and before its own Dirichlet partner receives the
reaction that the partner computed in the previous Schwarz iteration; the data
it returns then depend on two Schwarz iterations, and the least-squares problem
fits a relation that does not hold. List every Dirichlet side before its
Neumann partner in `domains`; setup warns otherwise. Two subdomains are not
affected in either order. On the five-subdomain inclusion of issue #233,
cyclic roles on the rectangle–rectangle interfaces or the inclusion listed last
did not converge in 512 Schwarz iterations, against 30 to 38 with every
Dirichlet side first. It keeps its histories per
interface, clears them at every stop and after a stall, and
applies only to stops with a single substep (windowed stops use the fixed
factor `relaxation parameter`, as for Aitken). On the cantilever it needed
5 to 9 iterations per stop where recursive Aitken needed 5.8 to 20.9 and the
fixed factor 0.5 needed 5.9 to 44.5 or diverged, and 14 to 16 with an
explicit Dirichlet side, where Aitken needed 34 to 35 and every fixed factor
diverged. The depth must cover the iterations a stop needs: with an explicit
Dirichlet side a depth of 5 needed 54 to 69 iterations per stop. Recursive
Aitken (`relaxation: aitken recursive`) converged in every case measured as
well and is accepted by the stall rule at fewer stops than the secant form;
where 0.5 is near the optimum, as for identical integrators on conforming
meshes, the fixed factor is equivalent. A pair converges when the jump of its matched
quantity and the force residual are both at or below `matching relative
tolerance`, or when the root mean square jump over the interface is at or
below `matching absolute tolerance` (in the units of the matched quantity;
default 0, which disables the absolute test). When every Schwarz coupling
is a matching pair this criterion alone stops the iteration, including at the
first Schwarz iteration of a stop, and the controller's `relative tolerance`
and `absolute tolerance` are not used; otherwise the update criterion
must hold as well. A residual that stays above the tolerance and has not
fallen 5% below its smallest value for three consecutive iterations is a
stall, treated according to `stalled interface jump action`; when the simulation has other couplings,
only iterations at which the update criterion holds are counted. A
single iteration without decrease, which Aitken factors produce, is not a
stall.

Velocity matching admits different time steps on the two sides, as in
Gravouil and Combescure (2001): the side with the finer step receives the
partner's velocity (Dirichlet side) or reaction (Neumann side) interpolated
linearly in time from the partner's substep history, the relaxation state is
kept per substep time, and the jump and force residual of the stopping rule
are the largest over the substeps of the stop. With different steps the side with the coarser step must be the Dirichlet
side: it imposes the fine side's velocity at the end of the stop, the fine
side receives its reaction interpolated linearly in time, and the interface
terms of the energy balance then cancel (the r = 1 multirate scheme of
Connors, Owen, Kuberry, and Bochev 2024, and that of Prakash and Hjelmstad
2004). On the cantilever with a 4:1 step ratio this kept Ẽ within 4e-12 and
E within 0.4% over 10 ms. With the finer step on the Dirichlet side Ẽ was
conserved only to the accumulated stopping error, E varied by up to 51%, and
the interface impulse residual (energy CSV) was of the order of the impulse.
On nonconforming meshes with identical integrators the finer mesh should also
be the Dirichlet side: on the 1:0.5 beam, giving the fine mesh the coarse step
and the Dirichlet role conserved Ẽ to 4e-12, while the coarse mesh as the
Dirichlet side diverged for explicit sides (gain radius 3.5 at relaxation
0.5) and needed 44.5 iterations per stop for implicit ones. Where the
stability limit of an explicit member does not allow the finer mesh the
coarser step, compute the gain. An explicit coarse side with an implicit
fine side has no good choice: as the Dirichlet side the explicit member
diverges. The
run aborts with
`swap BC types`, with HHT-α, with reduced order models, and with
`matching: displacement` unless both sides are Newmark with the same
β > 0, γ, and time step. The roles are those of the input; at setup the
code warns, without changing them, when the Dirichlet side has the finer
time step, when it is the explicit member of a mixed pair at equal steps,
or when it has fewer interface nodes than the partner with the same
integrator and step.

**Direct interface solve.** For a pair with `matching: velocity` of two central
difference subdomains at equal steps, `interface solve: direct` computes the interface
force without iteration. With m_D and m_N the lumped masses of the interface
rows and λ the force on the Dirichlet side (−Π_Dᵀ λ on the Neumann side), the
matching of the velocities at the end of the step gives, per component,
γ Δt (m_D⁻¹ + Π_D m_N⁻¹ Π_Dᵀ) λ = Π_D v_N,free − v_D,free, where v_free are the
velocities of both sides advanced without interface force. The matrix is
symmetric positive definite, of the size of the Dirichlet interface, and
factored once (dense Cholesky). Each stop advances both sides with free
interface rows, solves for λ, and corrects the interface accelerations and
velocities; the internal force is evaluated at the predictor displacement, so
the solve is exact for nonlinear materials as well. The same matrix gives the
initial interface force at t = 0 from the matching of the accelerations. Every
Schwarz coupling of the simulation must then be a direct pair, and each stop
must be one step of each subdomain (the explicit stable step must not shorten
it). Relaxation and the stopping tolerances do not apply.

The initial acceleration of a matching pair is found by a Schwarz
iteration at t = 0: the Dirichlet side imposes the projected partner
acceleration, the Neumann side receives the d'Alembert reaction, and each side
recomputes its initial acceleration with the interface rows included. The
Dirichlet datum is relaxed with `relaxation parameter`, under Aitken
relaxation with the factor of the selected form (recursive or secant) from the
iteration's own iterates, or under
Anderson acceleration by the same scheme on its own iterates (one
least-squares problem per pair, not stacked), and the
iteration stops when the
acceleration jump and the force residual are at or below `matching relative
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

### `Schwarz RR nonoverlap`

The classical Robin-Robin coupling t + α W u = g, with t the interface
traction, u the interface displacement, W the boundary mass matrix, α the
Robin parameter (force per unit area per unit displacement), and g the datum
assembled from the traction and displacement of the partner. Each side builds
its own transfer operators from the partner trace (per-side transfer), so the
two sides of an interface may use different values of α. The condition does
not absorb interface waves, so in elastodynamics it can increase the energy at the
interface (issue #176; on the cantilever and nested-cylinders benchmarks of
`docs/notes/schwarz-coupling` nearly every dynamic Robin-Robin run ends by
element inversion after exponential energy growth, on overlap and nonoverlap
decompositions alike). It is intended for quasi-statics and for comparison
against the classical Robin-Robin literature, and the run warns when it is used
with a dynamic time integrator. For dynamics use the matching
Dirichlet-Neumann exchange: `Schwarz DN nonoverlap` with `matching: velocity`,
with `interface solve: direct` for two central
difference subdomains at equal steps and `relaxation: anderson` otherwise.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner interface surface |
| `robin parameter` | yes | — | Robin coefficient α (positive); the two sides may use different values |

### `Schwarz contact`

Frictionless or tied contact enforced through Schwarz coupling. Frictionless
contact constrains the Cartesian x component of the interface nodes, so it
supports contact surfaces whose normal lies along x and aborts otherwise; tied
contact has no such restriction. Normals are computed on the current
configuration for finite kinematics and on the reference configuration for
infinitesimal kinematics.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source side set` | yes | — | partner contact surface |
| `friction type` | yes | — | `frictionless` or `tied` |
| `swap BC types` | no | `false` | swap the interface roles once per controller stop, before its Schwarz iteration |

## Mesh swapping (`swaps`)

Both single- and multi-domain runs can replace a mesh mid-run when a criterion
is met. `swaps` is a list; each entry names a `replacement` file and a
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
`initial conditions`, materials with internal state variables (`j2
plasticity`, `crystal plasticity` and `anisotropic hyperelastic damage`, whose
state the snapshot does not hold), `Schwarz contact`, and mesh smoothing.

```yaml
restart:
  index: -1
```

## Canonical examples

- Overlapping Schwarz: `examples/overlap/`
- Overlapping Schwarz with J2 plasticity (circular laser weld, graded pair):
  `examples/overlap/static-same-step/clw/`
- Non-overlapping Schwarz (Dirichlet-Neumann, matching Dirichlet-Neumann,
  and Robin-Robin): `examples/nonoverlap/dynamic-same-step/`,
  `examples/nonoverlap/static-same-step/`; the matching exchange on 2:1
  nonconforming meshes: `examples/nonoverlap/dynamic-same-step/cantilever-dn-nonconforming/`
- Contact: `examples/contact/`
- Adaptive mesh swapping: `examples/adaptive-time-stepping/`,
  `examples/ahead/`

See `docs/notes/schwarz-coupling` for the theory and stability analysis of the
Schwarz couplings.
