# The case matrix of the dynamic stability study: one definition of the
# problems, the factors, and the tiers, shared by generate.jl, run.jl, and
# collect.jl.  Change the study here, not in the scripts that use it.

# A case is one run.  The two letters of `pair` name the time integrator of
# each subdomain, the subdomain that holds the support first: the clamped
# part of the beam and the outer part of the cylinders.  I is implicit
# Newmark, E is explicit central difference.  A monolithic reference has a
# single letter.
struct Case
    problem::String   # "beam" or "cyl"
    coupling::String  # "mono", "ov-dir", "ov-imp", "no-dn", "no-imp"
    pair::String      # "II", "IE", "EI", "EE"; "I" or "E" for "mono"
    level::Int        # refinement level: element size and time step divided by it
    ratio::Float64    # beam only: clamped element size = free size / ratio
end

# The four couplings: the decomposition (overlap or nonoverlap) crossed with
# the transmission condition (standard or impedance).
const COUPLINGS = Dict(
    "ov-dir" => "overlap, Dirichlet (Schwarz overlap)",
    "ov-imp" => "overlap, impedance (Schwarz impedance overlap)",
    "no-dn" => "nonoverlap, Dirichlet-Neumann (Schwarz DN nonoverlap)",
    "no-imp" => "nonoverlap, paired impedance (Schwarz impedance nonoverlap)",
)
const COUPLING_ORDER = ["ov-dir", "ov-imp", "no-dn", "no-imp"]
const PAIRS = ["II", "IE", "EI", "EE"]

# Physical and numerical parameters of the two problems at refinement level 1.
# At level L the element sizes and the time step are divided by L, which
# keeps the Courant number of every run the same.
const BEAM = (
    length = 0.254,           # m
    height = 0.0254,
    width = 0.0127,
    h = 6.35e-3,              # element size of the free part (and of the reference)
    split = 0.127,            # nonoverlap interface at x = split
    overlap = 0.0254,         # overlap width, centered on the split
    elastic_modulus = 6.895e9,
    poissons_ratio = 0.25,
    density = 2768.0,
    release = "0.393701 * x * x",   # initial bent shape u_y(x)
    time_step = 1.0e-6,
    final_time = 1.0e-2,
    robin_parameter = 2.0e9,  # paired impedance, same on both sides
    relative_tolerance = 1.0e-8,
    absolute_tolerance = 2.54e-8,
    maximum_iterations = 128,
    newton_relative_tolerance = 1.0e-10,   # implicit subdomains
    newton_absolute_tolerance = 2.54e-8,
    checkpoints = [1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2],
)
const CYLINDERS = (
    time_step = 5.0e-7,
    final_time = 2.0e-3,
    launch_velocity = -35.0,  # m/s in z
    robin_parameter = 2.8e9,  # paired impedance: filler E / 76 mm, as in the coupling note
    relative_tolerance = 1.0e-6,
    absolute_tolerance = 1.0e-4,
    maximum_iterations = 256,
    # The residual of a cylinder subdomain stalls near 1e-6 N in round-off,
    # so the Newton tolerances are looser than on the beam.
    newton_relative_tolerance = 1.0e-8,
    newton_absolute_tolerance = 1.0e-5,
    checkpoints = [2.5e-4, 5.0e-4, 1.0e-3, 2.0e-3],
)

problem_parameters(problem::AbstractString) = problem == "beam" ? BEAM : CYLINDERS

function case_name(c::Case)
    name = "$(c.problem)_$(c.coupling)_$(c.pair)_L$(c.level)"
    if c.problem == "beam" && c.coupling != "mono"
        name *= "_r" * lpad(round(Int, 100 * c.ratio), 3, '0')
    end
    return name
end

coupled_cases(problem, level, pairs; ratios=[NaN]) =
    [Case(problem, coupling, pair, level, ratio) for ratio in ratios for coupling in COUPLING_ORDER for pair in pairs]

reference_cases(problem, level, integrators) = [Case(problem, "mono", i, level, NaN) for i in integrators]

# The tiers, in the order to run them.  Each is a complete set of cases with
# the references it is compared against.
const BEAM_RATIOS = [1.0, 0.75]
const TIERS = [
    "A" => "beam, level 1: every coupling and integrator pair at two mesh ratios",
    "B" => "beam, level 2: the same matrix refined",
    "C" => "cylinders, level 1: every coupling, explicit on both sides",
    "D" => "cylinders, level 1: every coupling with implicit subdomains",
    "E" => "cylinders, level 2: every coupling, explicit on both sides",
    "F" => "beam, level 4: the beam matrix refined twice (long runs)",
]

function tier_cases(tier::AbstractString)
    tier == "A" && return [coupled_cases("beam", 1, PAIRS; ratios=BEAM_RATIOS); reference_cases("beam", 1, ["I", "E"])]
    tier == "B" && return [coupled_cases("beam", 2, PAIRS; ratios=BEAM_RATIOS); reference_cases("beam", 2, ["I", "E"])]
    tier == "C" && return [coupled_cases("cyl", 1, ["EE"]); reference_cases("cyl", 1, ["E"])]
    tier == "D" && return [coupled_cases("cyl", 1, ["II", "IE", "EI"]); reference_cases("cyl", 1, ["I"])]
    tier == "E" && return [coupled_cases("cyl", 2, ["EE"]); reference_cases("cyl", 2, ["E"])]
    tier == "F" && return [coupled_cases("beam", 4, PAIRS; ratios=BEAM_RATIOS); reference_cases("beam", 4, ["I", "E"])]
    error("Unknown tier $tier; the tiers are $(join(first.(TIERS), ", "))")
end

all_cases() = reduce(vcat, [tier_cases(first(t)) for t in TIERS])

# Every case with a given name, across all tiers.
function find_case(name::AbstractString)
    for c in all_cases()
        case_name(c) == name && return c
    end
    return nothing
end
