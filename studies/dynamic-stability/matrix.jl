# The case matrix of the dynamic stability study: one definition of the
# problem, the factors, and the tiers, shared by generate.jl, run.jl,
# collect.jl, and plot.py.  Change the study here, not in the scripts.
#
# The subject is the constrained Dirichlet-Neumann exchange on the
# cantilever beam (coupling "no-cd": d'Alembert reaction, adjoint transfers
# from one cross-mass matrix, velocity constraint, stopping on the interface
# jump and force residual, coupled initial acceleration;
# docs/notes/schwarz-coupling, Section 9).  The plain Dirichlet-Neumann
# exchange ("no-dn") is the baseline that fails, and the monolithic run
# ("mono") the reference.
#
# The study varies one factor at a time around a default case, so that each
# figure of plot.py shows the effect of one factor with the others fixed.

# A case is one run.  The two letters of `pair` name the time integrator of
# each subdomain, the clamped part first: I is implicit Newmark (trapezoidal
# rule), E is explicit central difference.  A monolithic reference has a
# single letter.
#
#   level      refinement level: element size and time step divided by it
#   ratio      clamped element size = free size / ratio (1 is conforming)
#   tolerance  relative Schwarz tolerance: on the interface velocity jump and
#              the force residual (constrained), on the displacement update
#              (baseline)
#   steps      substeps of the free part per stop: 1 (equal steps) or 4
#   role       with steps > 1, the Dirichlet side: "coarse" (the clamped
#              part, which takes the coarse step) or "fine" (the free part,
#              which takes the fine step); "" at equal steps
#   solver     "anderson" (depth 10, mixing 0.5), "aitken" (recursive, from
#              0.5), "fixed" (relaxation 0.5), or "direct" (the direct
#              interface solve, two explicit sides at equal steps only)
struct Case
    coupling::String   # "mono", "no-dn", "no-cd"
    pair::String       # "II", "IE", "EI", "EE"; "I" or "E" for "mono"
    level::Int
    ratio::Float64
    tolerance::Float64
    steps::Int
    role::String
    solver::String
end

const DEFAULT_TOLERANCE = 1.0e-12
const DEFAULT_SOLVER = "anderson"

Case(coupling, pair, level, ratio; tolerance=DEFAULT_TOLERANCE, steps=1, role="", solver=DEFAULT_SOLVER) =
    Case(coupling, pair, level, ratio, tolerance, steps, role, solver)

const COUPLINGS = Dict(
    "mono" => "monolithic reference",
    "no-dn" => "Dirichlet-Neumann exchange (Schwarz DN nonoverlap), the baseline",
    "no-cd" => "constrained Dirichlet-Neumann exchange (Schwarz DN nonoverlap, constrained: true)",
)
const PAIRS = ["II", "IE", "EI", "EE"]
# The two pairs with one integrator on both sides anchor every factor; the
# mixed pairs enter the factors where the Dirichlet role matters.
const ANCHOR_PAIRS = ["II", "EE"]

# Physical and numerical parameters of the beam at refinement level 1.  At
# level L the element sizes and the time step are divided by L, which keeps
# the Courant number of every run the same.
const BEAM = (
    length = 0.254,           # m
    height = 0.0254,
    width = 0.0127,
    h = 6.35e-3,              # element size of the free part (and of the reference)
    split = 0.127,            # interface at x = split
    elastic_modulus = 6.895e9,
    poissons_ratio = 0.25,
    density = 2768.0,
    release = "0.393701 * x * x",   # initial bent shape u_y(x)
    time_step = 1.0e-6,
    final_time = 1.0e-2,
    maximum_iterations = 256,
    # The baseline exchange keeps the tolerances of the repository examples,
    # under which its failures were first measured.
    baseline_relative_tolerance = 1.0e-8,
    baseline_absolute_tolerance = 2.54e-8,
    newton_relative_tolerance = 1.0e-10,   # implicit subdomains
    newton_absolute_tolerance = 2.54e-8,
    checkpoints = [1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2],
)

function case_name(c::Case)
    name = "beam_$(c.coupling)_$(c.pair)_L$(c.level)"
    c.coupling == "mono" && return name
    name *= "_r" * lpad(round(Int, 100 * c.ratio), 3, '0')
    c.steps > 1 && (name *= "_s$(c.steps)" * (c.role == "coarse" ? "c" : "f"))
    c.tolerance != DEFAULT_TOLERANCE && (name *= "_tol" * string(round(Int, -log10(c.tolerance))))
    c.solver != DEFAULT_SOLVER && (name *= "_" * c.solver)
    return name
end

reference_cases(level) = [Case("mono", i, level, NaN) for i in ("I", "E")]

# The factors, each a tier.  Every tier includes the default case of its
# anchor pairs, so a tier can be plotted on its own.
const TIERS = [
    "A" => "coupling: monolithic, Dirichlet-Neumann, constrained; the four integrator pairs",
    "B" => "mesh ratio: 1, 0.75, 0.5",
    "C" => "Schwarz tolerance: 1e-6, 1e-8, 1e-10, 1e-12",
    "D" => "time steps: equal, 4:1 with the coarse-step side as Dirichlet, 4:1 with the fine-step side",
    "E" => "solver: Anderson, Aitken, fixed relaxation, direct solve",
    "F" => "refinement: levels 1, 2, 4 at fixed Courant number (long runs)",
]
const TIER_FACTOR = Dict("A" => "coupling", "B" => "ratio", "C" => "tolerance", "D" => "steps", "E" => "solver",
                         "F" => "level")

function tier_cases(tier::AbstractString)
    if tier == "A"
        return [[Case(cpl, pair, 1, 1.0) for cpl in ("no-dn", "no-cd") for pair in PAIRS]; reference_cases(1)]
    elseif tier == "B"
        return [Case("no-cd", pair, 1, r) for r in (1.0, 0.75, 0.5) for pair in ANCHOR_PAIRS]
    elseif tier == "C"
        return [Case("no-cd", pair, 1, 1.0; tolerance=tol) for tol in (1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12)
                for pair in ANCHOR_PAIRS]
    elseif tier == "D"
        cases = [Case("no-cd", pair, 1, 1.0) for pair in ("II", "EE", "IE")]
        append!(cases, [Case("no-cd", pair, 1, 1.0; steps=4, role=role) for role in ("coarse", "fine")
                        for pair in ("II", "EE", "IE")])
        return cases
    elseif tier == "E"
        cases = [Case("no-cd", pair, 1, 1.0; solver=s) for s in ("anderson", "aitken", "fixed") for pair in ANCHOR_PAIRS]
        push!(cases, Case("no-cd", "EE", 1, 1.0; solver="direct"))
        return cases
    elseif tier == "F"
        return [[Case("no-cd", pair, L, 1.0) for L in (1, 2, 4) for pair in ANCHOR_PAIRS];
                reference_cases(1); reference_cases(2); reference_cases(4)]
    end
    error("Unknown tier $tier; the tiers are $(join(first.(TIERS), ", "))")
end

all_cases() = unique(reduce(vcat, [tier_cases(first(t)) for t in TIERS]))

# The case with a given name, across all tiers, or nothing.
function find_case(name::AbstractString)
    for c in all_cases()
        case_name(c) == name && return c
    end
    return nothing
end
