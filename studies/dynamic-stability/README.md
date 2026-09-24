# Dynamic stability of Schwarz coupling: a parameter study

This directory holds the plan and the tools for a systematic study of the
energy behavior of Schwarz coupling in dynamics, on two problems, the bent
cantilever beam and the nested cylinders. The aim is to map trends: which
combinations of decomposition, transmission condition, time integrators,
and refinement conserve energy, which dissipate it, and which grow it until
the run fails.

Everything needed is here: `matrix.jl` defines the cases, `generate.jl`
writes their inputs and meshes, `run.jl` runs them, `collect.jl` summarizes
them in a table, and `plot.py` draws their energy histories. The runs and
their results stay out of the repository (see `.gitignore`); the findings
are compiled in a slide deck for the team (see Presenting the results).

## Background

Read these first, in this order:

1. `docs/notes/schwarz-coupling/schwarz-coupling.pdf`, Sections 1 to 3
   (the couplings and the energy metrics) and the measurement sections on
   the beam and the cylinders. The tables there are the results this study
   extends: overlap Dirichlet coupling grows energy until the mesh inverts,
   nonoverlap Dirichlet-Neumann with explicit subdomains conserves it, the
   paired impedance coupling is sign-definite but dissipates on the
   cylinders.
2. The input reference for multidomain runs, `docs/src/reference/multidomain.md`,
   for the meaning of every key in the generated inputs.
3. The examples the generated inputs are modeled on:
   `examples/overlap/dynamic-same-step/cantilever*`,
   `examples/nonoverlap/dynamic-same-step/cantilever-*`, and
   `examples/jmp/concentric-cylinders`.

## The two problems

**Cantilever beam.** A linear elastic aluminum bar, 254 × 25.4 × 12.7 mm
(E = 6.895 GPa, ν = 0.25, ρ = 2768 kg/m³), clamped at x = 0 and released at
rest from the bent shape u_y = 0.393701 x². It then vibrates in bending; the
total energy is 1506 J and must stay constant. The bar is split at
x = 127 mm; the overlap decomposition extends each part by 12.7 mm past the
split (overlap width 25.4 mm). The free part is meshed at h = 6.35 mm / L
(HEX8), the clamped part at (6.35 mm / L) / r, where L is the refinement
level and r the mesh ratio: r = 1 gives conforming meshes, r = 0.75 a
coarser clamped part. Horizon 10 ms, time step 1 µs / L.

**Nested cylinders.** A steel shell, a soft filler, and an offset core,
launched at −35 m/s against a support that fixes the bottom face in z
(`examples/jmp/concentric-cylinders`, where the geometry and materials are
described). The Schwarz interface is a closed surface inside the filler
around the core, with a lateral surface and two caps; the inner and outer
meshes are nonconforming (1:0.8). The meshes come from one Cubit journal
whose element sizes are divided by the level L. Horizon 2 ms, time step
0.5 µs / L. The filler has one fiftieth of the shell's wave impedance, so
this case is much harder on the coupling than the beam.

At every level the element size and the time step are divided together,
which keeps the Courant number of every run the same; the explicit runs use
a Courant number of 0.5.

## The matrix

Each run is one combination of:

| Factor | Levels |
|---|---|
| problem | beam, cylinders |
| coupling | `ov-dir` overlap Dirichlet, `ov-imp` overlap impedance, `no-dn` nonoverlap Dirichlet-Neumann, `no-imp` nonoverlap paired impedance |
| integrator pair | `II`, `IE`, `EI`, `EE`: implicit Newmark (I) or explicit central difference (E) on each subdomain, the subdomain with the support first (clamped part of the beam, outer part of the cylinders) |
| refinement level | 1, 2, 4 (beam); 1, 2 (cylinders) |
| mesh ratio | 1.0, 0.75 (beam only) |

Every set of coupled runs is compared with a monolithic reference, one mesh
and no coupling, at the same level and with each integrator (`mono`, `I` or
`E`). The runs are grouped in tiers, to be done in order:

| Tier | Cases | What it adds | Estimated cost |
|---|---|---|---|
| A | 34 | beam, level 1: 4 couplings × 4 pairs × 2 ratios, 2 references | 2 to 20 min per run |
| B | 34 | beam, level 2: the same, refined | about 16 times tier A per run |
| C | 5 | cylinders, level 1, explicit on both sides | 20 min to 3 h per run on 4 threads |
| D | 13 | cylinders, level 1, with implicit subdomains (`II`, `IE`, `EI`) | 10 to 55 h per run on 2 threads |
| E | 5 | cylinders, level 2, explicit on both sides | about 16 times tier C per run |
| F | 34 | beam, level 4 | days per run; optional |

The costs are estimates from short runs, measured with four runs at a time
on a 16-core workstation; record the actual wall time (it is in the summary)
and revise them after tier A. The implicit cylinder runs are expensive
because every Newton iteration factors the stiffness of a subdomain of 20 to
60 thousand elements: 10 to 20 s per stop for Dirichlet-Neumann and overlap
Dirichlet, 30 to 50 s for the impedance couplings, and about 2 min for the
implicit monolithic reference. Run tier D with more threads per case and
fewer cases at a time, start with `II` and the reference, and consider a
horizon of 1 ms (`--final-time 1.0e-3`), which on the explicit runs already
separates growth from conservation. A level 2 cylinder run spends about
nine minutes in setup before its first stop, then about 3 s per stop on
4 threads for the explicit Dirichlet-Neumann coupling. Refining once multiplies the
elements by 8 and the steps by 2, so every level costs about 16 times the
one before.

Names: `<problem>_<coupling>_<pair>_L<level>[_r<ratio>]`, for example
`beam_no-imp_IE_L1_r075` or `cyl_ov-dir_EE_L2`. To list the cases of a tier
without writing anything:

    julia --project=../.. generate.jl --list A

### Fixed settings

The generated inputs fix everything else at the values the earlier studies
used, so the results can be compared with the coupling note:

- Schwarz tolerances: relative 1e-8 and absolute 2.54e-8 on the beam,
  relative 1e-6 and absolute 1e-4 on the cylinders; at most 128 and 256
  Schwarz iterations per stop.
- Relaxation of the Schwarz iteration: none for the overlap couplings;
  recursive Aitken for Dirichlet-Neumann; for the paired impedance a fixed
  factor of 0.5 when either subdomain is explicit (Aitken diverges on the
  explicit cylinders) and recursive Aitken on the implicit beam (a tenth of
  the iterations of the fixed factor).
- Robin parameter: 0 for the overlap impedance; 2.0e9 (beam) and 2.8e9
  (cylinders) for the paired impedance, the same on both sides. It changes
  the convergence rate of the Schwarz iteration, not the converged answer.
- The subdomain without the support is listed first in `run.yaml`, so the
  relaxation acts on the support side, as in the repository examples.
- Newton on implicit subdomains: relative 1e-10 and absolute 2.54e-8 on the
  beam, 1e-8 and 1e-5 on the cylinders (whose residual stalls near 1e-6 N
  in round-off), with a direct linear solver on the cylinders.

Each of these is a single line in `matrix.jl` or `generate.jl`. Changing one
is an experiment of its own: do it in a separate runs directory (`--runs`),
never in the middle of a tier.

## Workflow

All commands run from this directory. `../..` is the Norma project.

1. **Check the setup.** Generate a tier with a very short horizon in a scratch
   directory and run it; every case should complete in about a minute, most of
   it compilation.

       julia --project=../.. generate.jl --runs smoke --final-time 2.0e-5 A C
       julia --project=../.. run.jl --runs smoke --jobs 4 --threads 2 A C
       julia --project=../.. collect.jl --runs smoke --output smoke.csv

2. **Generate a tier.**

       julia --project=../.. generate.jl A

   Each case directory under `runs/` holds `run.yaml` (the multidomain
   input), one input per subdomain, the meshes, and `case.txt`. Open a few and
   read them against the input reference: this is exactly what Norma runs.
   The beam meshes are written by `meshes.jl`; the cylinder meshes at level 1
   are the committed ones, and at level 2 `generate.jl` runs the Cubit journal
   (put `cubit` on PATH or set `CUBIT` to its path; it takes a while, once
   per level).

3. **Run it.**

       julia --project=../.. run.jl --jobs 4 --threads 4 A

   Choose `--jobs` × `--threads` no larger than the number of cores. Each run
   writes `run.log`, `run-energy.csv` (the blended energy at every stop), and
   `status.txt`. Rerunning the same command skips the cases that completed, so
   an interrupted tier resumes where it stopped; `--force` reruns them. A case
   can also be run by hand, which is the way to look at one in detail:

       julia --project=../.. -t 4 run_case.jl runs/beam_ov-dir_II_L1_r075

   Exodus output is off to save disk (a cylinder output file reaches a
   gigabyte). To look at a run in ParaView, regenerate that case with
   `--exodus-interval 1.0e-5` and run it again.

4. **Summarize.**

       julia --project=../.. collect.jl
       python3 plot.py            # add --log to see exponential growth as a line

   `summary.csv` has one row per case (see the header of `collect.jl` for the
   columns) and `figures/` one figure per problem, level, and ratio.

## What to look at

For each run the summary gives the energy ratio E(t)/E(0) at four
checkpoints, its maximum and minimum, the first times it exceeds 1.1 and 2,
the Schwarz iterations per stop, and an outcome:

- **conserving**: completed, and the ratio stayed within [0.95, 1.05] at the
  end and never exceeded 1.05;
- **dissipating**: completed, but ended below 0.95;
- **growing**: completed, but exceeded 1.05 at some time;
- **failed**: stopped before the horizon, usually by an inverted element after
  growth, or by the Schwarz iteration reaching its limit.

The monolithic references must conserve (to a fraction of a percent); if
they do not, stop and find out why before looking at anything else.

Questions to answer, tier by tier:

1. Which couplings conserve, dissipate, or grow, and does the answer depend
   on the integrator pair? On the mesh ratio?
2. For the growing runs, is the growth exponential (straight on the
   logarithmic plot)? What is its rate (the e-folding time, from the slope),
   and how does it change with refinement? Growth that speeds up with
   refinement at a fixed Courant number points at the spatial coupling;
   growth that does not points elsewhere.
3. For the dissipating runs, does the loss shrink with refinement (a
   discretization error that converges) or stay the same (a property of the
   coupling)?
4. Do the two problems agree? Where they do not, what is different about
   them (impedance contrast, interface shape, loading)?
5. How many Schwarz iterations does each coupling need, and does that change
   over the run? A count that rises before a failure is a warning sign.

Compare with the tables of the coupling note wherever the cases coincide
(level 1): reproducing a published number is the best check that the setup
is right.

## Suggested order and milestones

1. Read the background and run the smoke check. Open one generated case of
   each coupling and explain every key to yourself.
2. Run the two beam references of tier A and confirm they conserve. Run one
   coupled case by hand and follow it in `run.log`.
3. Run the rest of tier A, collect, plot, and draft the first slides:
   what conserves, what grows, what dissipates, and how fast.
4. Tier B, then compare the two levels (question 2 and 3 above).
5. Tier C, then compare the cylinders with the beam (question 4).
6. Tiers D and E as time allows; tier F only if there is a clear reason.

Keep a log of every run set: the date, the Norma commit (`git log -1
--oneline`), the machine, the command, and anything unexpected. When a run
fails in a way that the outcome classes do not explain, keep its directory
and its log and ask.

## Presenting the results

Compile the results in a slide deck to present to the team, and add to it
as each tier is completed rather than at the end. A suggested outline:

1. The question and the two problems, with a picture of each mesh and
   interface.
2. The matrix: the factors, the tiers that were run, and the fixed settings.
3. For each problem and level, the figure from `plot.py` (linear and, where
   there is growth, logarithmic) and a summary table of the outcomes: one
   row per coupling, one column per integrator pair, each cell the outcome
   and the energy ratio at the horizon or the time of failure.
4. The trends: the effect of the coupling, of the integrators, of the mesh
   ratio, and of refinement, each stated in one sentence with the figure
   that shows it.
5. Agreement and disagreement with the tables of the coupling note, and
   between the two problems.
6. The cost of each coupling (Schwarz iterations per stop and wall time).
7. Open questions and proposed next experiments.

Every number on a slide should trace back to a row of `summary.csv`; note on
the slide the Norma commit the runs used.

## Beyond the matrix

Once the matrix is done, natural next experiments, each a change of one
fixed setting in a separate runs directory: the relaxation (fixed factors,
the two Aitken forms), the Robin parameter of the paired impedance, the
overlap width of the beam, more mesh ratios, the classical Robin-Robin
coupling (`Schwarz RR nonoverlap`), and different time steps on the two
subdomains (subcycling).

## Files

| File | Purpose |
|---|---|
| `matrix.jl` | the cases, the problem parameters, and the tiers |
| `meshes.jl` | beam meshes (written directly) and cylinder meshes (Cubit) |
| `generate.jl` | writes the case directories |
| `run.jl`, `run_case.jl` | run the cases, several at a time |
| `collect.jl` | summary table |
| `plot.py` | energy histories (needs matplotlib) |
