# Dynamic stability of Schwarz coupling: the constrained exchange on the cantilever beam

This directory holds the plan and the tools for a parameter study of the
constrained Dirichlet-Neumann exchange, the Schwarz coupling for dynamics
that conserves energy by construction, on the bent cantilever beam. The aim
is to show, with energy histories, the effect of each parameter of the
coupled problem on conservation and on cost: the coupling itself against the
exchange it replaces, the mesh ratio, the Schwarz tolerance, the time steps
of the two sides, the solver of the interface problem, and the refinement.

Everything needed is here. `study.sh` runs the whole study in one command:
`matrix.jl` defines the cases, `generate.jl` writes their inputs and meshes,
`run.jl` runs them, `collect.jl` summarizes them in a table, and `plot.py`
draws one figure per parameter. The runs and their results stay out of the
repository (see `.gitignore`); the findings are compiled in a slide deck for
the team.

## Background

Read these first, in this order:

1. `docs/notes/schwarz-coupling/schwarz-coupling.pdf`, Section 9: what the
   constrained exchange is (its five changes against the Dirichlet-Neumann
   exchange), why it conserves (the pseudo-energy identity and its
   cancellation), and what was measured on the beam, including the Aitken
   and Anderson subsections; then the conclusions. The tables there are the
   results this study extends; reproducing one of its numbers is the best
   check that the setup is right.
2. The input reference for multidomain runs, `docs/src/reference/multidomain.md`,
   for the meaning of every key in the generated inputs, and the example the
   inputs are modeled on, `examples/nonoverlap/dynamic-same-step/cantilever-dn`.

## The problem

A linear elastic aluminum bar, 254 × 25.4 × 12.7 mm (E = 6.895 GPa,
ν = 0.25, ρ = 2768 kg/m³), clamped at x = 0 and released at rest from the
bent shape u_y = 0.393701 x². It then vibrates in bending with a period of
about 10 ms; the total energy is 1506 J and must stay constant. The bar is
split at x = 127 mm into a clamped part and a free part, each meshed with
HEX8 elements of size 6.35 mm / L, where L is the refinement level; the
clamped part may be coarsened by the mesh ratio r (its size is
6.35 mm / (L r)). Horizon 10 ms, stops of 1 µs / L. At every level the
element size and the time step are divided together, which keeps the
Courant number of every run the same; the explicit runs use a Courant
number of 0.5.

Two energies are written at every stop (`run-energy.csv`): the energy E,
kinetic plus strain, which the plots show, and the pseudo-energy, the energy
functional applied to the velocity and the acceleration (`pseudo_energy_total`,
in J/s²), which is the quantity that the velocity constraint conserves by
theorem and which the summary table reports.

## The factors

The study varies one factor at a time around a default case: the
constrained exchange, level 1, conforming meshes, Schwarz tolerance 1e-12,
equal time steps, Anderson acceleration. Each factor is a tier, and each
tier is one figure.

| Tier | Factor | Levels | Pairs | Cases |
|---|---|---|---|---|
| A | coupling | monolithic; Dirichlet-Neumann exchange (the baseline, with the tolerances of the examples); constrained exchange | II, IE, EI, EE | 10 |
| B | mesh ratio | 1, 0.75, 0.5 | II, EE | 6 |
| C | Schwarz tolerance | 1e-6, 1e-8, 1e-10, 1e-12 | II, EE | 8 |
| D | time steps | equal; 4:1 with the clamped part (coarse step) as the Dirichlet side; 4:1 with the free part (fine step) as the Dirichlet side | II, EE, IE | 9 |
| E | solver | Anderson (depth 10); Aitken, recursive; fixed relaxation 0.5; direct interface solve (EE only) | II, EE | 7 |
| F | refinement | levels 1, 2, 4, with the monolithic references of each level | II, EE | 12 |

The integrator pair names the integrator of each part, the clamped part
first: I is implicit Newmark (trapezoidal rule), E is explicit central
difference. The default cases are shared between the tiers, so the study
has 39 distinct cases. Tiers A to E take minutes per case (2 to 20 min on 4
threads; the subcycled and the fixed-relaxation cases up to an hour); tier
F at level 4 takes days per case and is optional.

### Fixed settings

The generated inputs fix everything else; each setting is one line in
`matrix.jl` or `generate.jl`, and changing one is an experiment of its own,
to be done in a separate runs directory (`--runs`), never in the middle of
a tier.

- The Dirichlet side (the side that receives the velocity) is the implicit
  member of a mixed pair, otherwise the free part, which is the finer mesh
  when the clamped part is coarsened; with different time steps it is the
  factor under study. These are the role rules of the coupling note; with
  the explicit member as the Dirichlet side only Anderson acceleration
  converges.
- The constrained exchange stops on the interface velocity jump and the
  interface force residual at the relative tolerance of the case, with no
  absolute floor; the baseline keeps the displacement criterion of the
  examples (relative 1e-8, absolute 2.54e-8). At most 256 Schwarz
  iterations per stop; a stalled iteration is accepted with a warning.
- Newton on implicit subdomains: relative 1e-10, absolute 2.54e-8, direct
  linear solver.
- The free part is listed first in `run.yaml`, so the relaxation acts on the
  clamped side, as in the repository examples.

Names: `beam_<coupling>_<pair>_L<level>_r<ratio>` with the suffixes `_s4c`
or `_s4f` (4:1 steps, clamped or free part as the Dirichlet side), `_tol<n>`
(tolerance 1e-n), and `_aitken`, `_fixed`, or `_direct` for a solver other
than Anderson; for example `beam_no-cd_EE_L1_r100_s4f`. To list the cases of
a tier without writing anything:

    julia --project=../.. generate.jl --list D

## Running the study on rigel

Setup, once (`JULIA_PKG_SERVER=""` for every Pkg operation on rigel):

    cd ~/Repos/Norma.jl && git fetch && git checkout main && git merge --ff-only origin/main
    JULIA_PKG_SERVER="" julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
    cd test && julia --project=.. runtests.jl 135 136      # the energy output and the constrained exchange

Then, from this directory:

1. **Smoke check** (about ten minutes): every case of tiers A to E for 20
   stops, in a scratch directory, to confirm that the setup, the summary,
   and the figures work before the long runs start.

       ./study.sh --runs smoke --final-time 2.0e-5 A B C D E

2. **The study**, detached, with the log to follow:

       setsid nohup ./study.sh > study.log 2>&1 &
       tail -f study.log

   The default is 8 cases at a time with 4 threads each. The script runs the
   tiers in order, and after each tier rewrites `summary.csv` and the
   figures, so the results of a tier can be looked at while the next one
   runs. It is resumable: if rigel reboots, the same command continues from
   the cases that have not completed.

3. **One case by hand**, to look at it in detail:

       julia --project=../.. -t 4 run_case.jl runs/beam_no-cd_II_L1_r100

   Exodus output is off to save disk. To look at a run in ParaView,
   regenerate that case with `--exodus-interval 1.0e-5` and run it again.

Keep a log of every run set: the date, the Norma commit (`git log -1
--oneline`; `study.log` records it), the machine, the command, and anything
unexpected. When a run fails in a way that the outcome classes do not
explain, keep its directory and its log and ask.

## What to look at

`summary.csv` has one row per case (the header of `collect.jl` lists the
columns): the energy ratio at 1, 2, 5, and 10 ms, its extremes, the largest
change of the energy and of the pseudo-energy, the interface jumps at the
end, the impulse residual of subcycled cases, the Schwarz iterations per
stop, the wall time, and an outcome:

- **conserving to roundoff**: completed, and the energy ratio stayed within
  1e-8 of 1 over the whole run;
- **conserving**: completed, and the ratio stayed within [0.95, 1.05] at the
  end and never exceeded 1.05;
- **dissipating**: completed, but ended below 0.95;
- **growing**: completed, but exceeded 1.05 at some time;
- **failed**: stopped before the horizon, usually by an inverted element
  after growth, or by the Schwarz iteration reaching its limit.

`figures/` has one figure per factor, a panel per integrator pair, a curve
per level of the factor, the monolithic reference dashed; `*_log.png` are
the same on a logarithmic axis, where exponential growth is a straight line.

The monolithic references must conserve to roundoff; if they do not, stop
and find out why before looking at anything else. The expected outcomes,
from the coupling note (reproduce them first):

1. **Coupling**: the baseline inverts an element within 1 ms with implicit
   subdomains and oscillates within a few percent with explicit ones; the
   constrained exchange conserves to roundoff for every pair.
2. **Mesh ratio**: the constrained exchange conserves at every ratio.
3. **Tolerance**: the energy drift grows with the tolerance; at 1e-12 the
   pseudo-energy drifts by about 1e-10 over the 10⁴ stops. The question is
   how the drift scales with the tolerance and with the number of stops.
4. **Time steps**: with the clamped (coarse-step) part as the Dirichlet side
   the pseudo-energy is conserved and the energy stays within 0.4%; with the
   free (fine-step) part as the Dirichlet side the pseudo-energy is still
   conserved but the energy swings by 12% (II) and 31% (EE) over the bending
   period.
5. **Solver**: all converged solvers give the same energies; they differ in
   iterations per stop (Anderson 5 to 16, Aitken up to 35, fixed relaxation
   diverges with an explicit Dirichlet side) and in wall time; the direct
   solve conserves to roundoff with no iteration.
6. **Refinement**: the question is whether the drift per stop and the
   iteration count change with the mesh at fixed Courant number.

## Presenting the results

Compile the results in a slide deck for the team, one slide per factor with
its figure and the rows of `summary.csv` behind it, and add to it as each
tier is completed rather than at the end. Every number on a slide should
trace back to a row of `summary.csv`; note on the slide the Norma commit the
runs used. Close with the trends in one sentence each, the agreement or
disagreement with the tables of the coupling note, the cost of each solver,
and the open questions.

## Beyond the matrix

Natural next experiments, each a change of one fixed setting in a separate
runs directory: the displacement constraint (`constraint: displacement`,
trapezoidal Newmark pairs only), other step ratios (2:1, 8:1), the mixed
pairs at ratio 0.5 with subcycling, the Anderson depth and mixing parameter,
and the nested cylinders (`examples/jmp/concentric-cylinders` and its linear
version) with the same factors.

## Files

| File | Purpose |
|---|---|
| `study.sh` | the whole study in one command, resumable |
| `matrix.jl` | the cases, the problem parameters, the factors, and the tiers |
| `meshes.jl` | beam meshes (written directly, no Cubit) |
| `generate.jl` | writes the case directories |
| `run.jl`, `run_case.jl` | run the cases, several at a time |
| `collect.jl` | summary table |
| `plot.py` | one energy figure per factor (needs matplotlib) |
