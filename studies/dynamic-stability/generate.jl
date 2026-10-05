# Write the input files and meshes of the dynamic stability study.
#
#   julia --project=../.. generate.jl [options] TIER|CASE ...
#
# Every argument is a tier letter (see matrix.jl) or a case name; each case
# gets a directory under the runs directory with run.yaml (the top-level
# multidomain input), one input per subdomain, its meshes, and case.txt
# (the parameters of the case, which collect.jl and plot.py read).  An
# existing case directory is rewritten unless it has completed.
#
# Options:
#   --runs DIR           directory of the case directories (default: runs)
#   --final-time T       end every case at T instead of the beam's final
#                        time, for a quick check of a setup (for example 1.0e-5)
#   --exodus-interval T  write Exodus output every T seconds (default 0: none;
#                        the energy history is written in any case)
#   --list               print the names of the selected cases and write nothing

include(joinpath(@__DIR__, "matrix.jl"))
include(joinpath(@__DIR__, "meshes.jl"))

# A number as YAML reads it as a float: YAML 1.1 needs a sign in the
# exponent, so 2.0e9 would be read as a string, while 2.0e+9 is a number.
function num(x::Real)
    s = string(Float64(x))
    return occursin(r"e\d", s) ? replace(s, "e" => "e+") : s
end

function parse_arguments(args)
    options = Dict{String,Any}(
        "runs" => joinpath(@__DIR__, "runs"), "final time" => nothing, "exodus interval" => 0.0, "list" => false
    )
    selections = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--runs"
            options["runs"] = abspath(args[i + 1]); i += 2
        elseif arg == "--final-time"
            options["final time"] = parse(Float64, args[i + 1]); i += 2
        elseif arg == "--exodus-interval"
            options["exodus interval"] = parse(Float64, args[i + 1]); i += 2
        elseif arg == "--list"
            options["list"] = true; i += 1
        else
            push!(selections, arg); i += 1
        end
    end
    isempty(selections) && error("Give one or more tiers ($(join(first.(TIERS), ", "))) or case names")
    return selections, options
end

function selected_cases(selections)
    cases = Case[]
    for s in selections
        c = find_case(s)
        append!(cases, c === nothing ? tier_cases(s) : [c])
    end
    return unique(cases)
end

function time_integrator(kind::Char, time_step)
    kind == 'I' && return """
        time integrator:
          type: Newmark
          β: 0.25
          γ: 0.5
          time step: $(num(time_step))
        """
    return """
        time integrator:
          type: central difference
          time step: $(num(time_step))
          CFL: 0.5
          γ: 0.5
        """
end

function solver(kind::Char)
    kind == 'E' && return """
        solver:
          type: explicit solver
          step: explicit
        """
    return """
        solver:
          type: Hessian minimizer
          step: full Newton
          linear solver: direct
          minimum iterations: 1
          maximum iterations: 16
          relative tolerance: $(num(BEAM.newton_relative_tolerance))
          absolute tolerance: $(num(BEAM.newton_absolute_tolerance))
        """
end

materials(domain) = """
    model:
      type: solid mechanics
      material:
        blocks:
          $domain: elastic
        elastic:
          model: linear elastic
          elastic modulus: $(num(BEAM.elastic_modulus))
          Poisson's ratio: $(num(BEAM.poissons_ratio))
          density: $(num(BEAM.density))
    """

const INITIAL_CONDITIONS = """
    initial conditions:
      displacement:
        - node set: nsall
          component: y
          function: "$(BEAM.release)"
    """

const SUPPORT = """
      Dirichlet:
        - node set: nsx-
          component: x
          function: "0.0"
        - node set: nsx-
          component: y
          function: "0.0"
        - node set: nsx-
          component: z
          function: "0.0"
    """

# The Dirichlet side of a pair (the side that receives the velocity).  At
# equal steps: the implicit member of a mixed pair, otherwise the free part,
# which is the finer mesh when the clamped part is coarsened.  With
# different steps the role is the factor under study: "coarse" is the
# clamped part (coarse step), "fine" the free part (fine step).  With the
# explicit member as the Dirichlet side the gain of the iteration has
# eigenvalues down to -13 and only Anderson acceleration converges
# (docs/notes/schwarz-coupling).
function dirichlet_domain(c::Case)
    c.steps > 1 && return c.role == "coarse" ? "clamped" : "free"
    clamped_kind, free_kind = c.pair[1], c.pair[2]
    clamped_kind == 'I' && free_kind == 'E' && return "clamped"
    return "free"
end

function coupling_condition(c::Case, domain)
    partner = domain == "clamped" ? "free" : "clamped"
    side_set = domain == "clamped" ? "ssx+" : "ssx-"
    partner_side_set = partner == "clamped" ? "ssx+" : "ssx-"
    bc_type = domain == dirichlet_domain(c) ? "Dirichlet" : "Neumann"
    text = """
          Schwarz DN nonoverlap:
            - side set: $side_set
              source: $partner
              source side set: $partner_side_set
              default BC type: $bc_type
        """
    c.coupling == "no-dn" && return text
    text *= """
              constrained: true
              constraint: velocity
        """
    c.solver == "direct" && (text *= "      interface solve: direct\n")
    return text
end

# The time step of a subdomain: the free part takes the fine step when the
# case subcycles.
subdomain_step(c::Case, domain) = BEAM.time_step / c.level / (domain == "free" ? c.steps : 1)

function subdomain_input(c::Case, domain, kind::Char)
    text = "type: single\ninput mesh file: $domain.g\noutput mesh file: $domain.e\n"
    text *= materials(domain)
    text *= time_integrator(kind, subdomain_step(c, domain))
    text *= INITIAL_CONDITIONS
    conditions = ""
    domain in ("beam", "clamped") && (conditions *= SUPPORT)
    c.coupling != "mono" && (conditions *= coupling_condition(c, domain))
    isempty(conditions) || (text *= "boundary conditions:\n" * conditions)
    return text * solver(kind)
end

# Relaxation of the Schwarz iteration.  The baseline keeps recursive Aitken,
# the setting of the repository examples.  For the constrained exchange the
# solver is a factor: Anderson acceleration (depth 10, mixing 0.5) is the
# default and converges in every case measured; recursive Aitken converges
# everywhere but needs more iterations where the best factor is far from
# 0.5; the fixed factor 0.5 diverges with an explicit Dirichlet side; the
# direct solve needs no iteration.
function relaxation(c::Case)
    c.coupling == "no-dn" && return "relaxation: aitken recursive\n"
    c.solver == "anderson" && return "relaxation: anderson\nanderson depth: 10\nrelaxation parameter: 0.5\n"
    c.solver == "aitken" && return "relaxation: aitken recursive\nrelaxation parameter: 0.5\n"
    return "relaxation parameter: 0.5\n"
end

# Schwarz tolerances: the constrained criterion is relative, on the
# interface velocity jump and the force residual, with no absolute floor;
# the baseline keeps the displacement criterion of the examples.
function tolerances(c::Case)
    c.coupling == "no-dn" && return (BEAM.baseline_relative_tolerance, BEAM.baseline_absolute_tolerance)
    return (c.tolerance, 1.0e-15)
end

function run_input(c::Case, domains, final_time, exodus_interval)
    relative, absolute = tolerances(c)
    list = join(("\"$d.yaml\"" for d in domains), ", ")
    text = """
        type: multi
        domains: [$list]
        blended energy output: true
        Exodus output interval: $(num(exodus_interval))
        CSV output interval: 0
        initial time: 0.0
        final time: $(num(final_time))
        time step: $(num(BEAM.time_step / c.level))
        minimum iterations: 1
        maximum iterations: $(BEAM.maximum_iterations)
        relative tolerance: $(num(relative))
        absolute tolerance: $(num(absolute))
        """
    c.coupling == "mono" && return text
    return text * relaxation(c) * "stalled interface jump action: warn\n"
end

function write_case(c::Case, options)
    dir = joinpath(options["runs"], case_name(c))
    if isfile(joinpath(dir, "status.txt")) && occursin("status: completed", read(joinpath(dir, "status.txt"), String))
        return dir, false
    end
    isdir(dir) && rm(dir; recursive=true)
    mkpath(dir)
    final_time = something(options["final time"], BEAM.final_time)
    if c.coupling == "mono"
        domains, kinds = ["beam"], [c.pair[1]]
    else
        # The free part is listed first, so the Schwarz relaxation acts on
        # the clamped side, as in the repository examples.
        domains, kinds = ["free", "clamped"], [c.pair[2], c.pair[1]]
    end
    for (domain, kind) in zip(domains, kinds)
        write(joinpath(dir, "$domain.yaml"), subdomain_input(c, domain, kind))
    end
    write(joinpath(dir, "run.yaml"), run_input(c, domains, final_time, options["exodus interval"]))
    write_beam_meshes(c, dir)
    write(joinpath(dir, "case.txt"), """
        name: $(case_name(c))
        coupling: $(c.coupling)
        pair: $(c.pair)
        level: $(c.level)
        ratio: $(c.ratio)
        tolerance: $(num(c.tolerance))
        steps: $(c.steps)
        role: $(c.role)
        solver: $(c.solver)
        dirichlet side: $(c.coupling == "mono" ? "" : dirichlet_domain(c))
        time step: $(num(BEAM.time_step / c.level))
        final time: $(num(final_time))
        """)
    return dir, true
end

function main(args)
    selections, options = parse_arguments(args)
    cases = selected_cases(selections)
    if options["list"]
        foreach(c -> println(case_name(c)), cases)
        println("$(length(cases)) cases")
        return
    end
    written = 0
    for c in cases
        dir, fresh = write_case(c, options)
        println("  ", basename(dir), fresh ? "" : "  (completed, kept)")
        written += fresh
    end
    println("$written cases written under $(options["runs"]) ($(length(cases) - written) completed and kept)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
