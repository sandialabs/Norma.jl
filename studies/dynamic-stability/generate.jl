# Write the input files and meshes of the dynamic stability study.
#
#   julia --project=../.. generate.jl [options] TIER|CASE ...
#
# Every argument is a tier letter (see matrix.jl) or a case name; each case
# gets a directory under the runs directory with run.yaml (the top-level
# multidomain input), one input per subdomain, its meshes, and case.txt
# (the parameters of the case, which collect.jl reads).  An existing case
# directory is rewritten.
#
# Options:
#   --runs DIR           directory of the case directories (default: runs)
#   --final-time T       end every case at T instead of its problem's final
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

# Names of the two subdomains of a problem: the one that holds the support
# first, as in the integrator pair.
support_domain(problem) = problem == "beam" ? "clamped" : "outer"
other_domain(problem) = problem == "beam" ? "free" : "inner"

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

# Newton with a direct linear solver on the cylinders, whose subdomains are
# large enough for conjugate gradients to dominate the cost.
function solver(kind::Char, problem)
    p = problem_parameters(problem)
    kind == 'E' && return """
        solver:
          type: explicit solver
          step: explicit
        """
    linear = problem == "cyl" ? "  linear solver: direct\n" : ""
    return """
        solver:
          type: Hessian minimizer
          step: full Newton
          minimum iterations: 1
          maximum iterations: 16
          relative tolerance: $(num(p.newton_relative_tolerance))
          absolute tolerance: $(num(p.newton_absolute_tolerance))
        """ * linear
end

function materials(problem, domain)
    if problem == "beam"
        return """
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
    end
    blocks = Dict(
        "cylinders" => ["core", "filler", "shell"], "inner" => ["core", "filler"], "outer" => ["filler", "shell"]
    )
    definitions = Dict(
        "core" => """
                core:
                  model: neohookean
                  elastic modulus: 100.0e+09
                  Poisson's ratio: 0.33
                  density: 4915.0
            """,
        "filler" => """
                filler:
                  model: neohookean
                  bulk modulus: 5.99e+09
                  shear modulus: 7.1658e+07
                  density: 1843.0
            """,
        "shell" => """
                shell:
                  model: neohookean
                  elastic modulus: 200.0e+09
                  Poisson's ratio: 0.28963
                  density: 7822.8
            """,
    )
    text = "model:\n  type: solid mechanics\n  material:\n    blocks:\n"
    for b in blocks[domain]
        text *= "      $b: $b\n"
    end
    return text * join(definitions[b] for b in blocks[domain])
end

function initial_conditions(problem)
    problem == "beam" && return """
        initial conditions:
          displacement:
            - node set: nsall
              component: y
              function: "$(BEAM.release)"
        """
    return """
        initial conditions:
          velocity:
            - node set: nsall
              component: z
              function: "$(CYLINDERS.launch_velocity)"
        """
end

function support(problem)
    problem == "beam" && return """
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
    return """
          Dirichlet:
            - node set: nsz-
              component: z
              function: "0.0"
        """
end

# The coupling condition of one subdomain toward its partner.
function coupling_condition(c::Case, domain)
    beam = c.problem == "beam"
    partner = domain == support_domain(c.problem) ? other_domain(c.problem) : support_domain(c.problem)
    side_set = beam ? (domain == "clamped" ? "ssx+" : "ssx-") : "ssint"
    partner_side_set = beam ? (partner == "clamped" ? "ssx+" : "ssx-") : "ssint"
    partner_block = beam ? partner : "filler"
    c.coupling == "ov-dir" && return """
          Schwarz overlap:
            - side set: $side_set
              source: $partner
              source block: $partner_block
        """
    c.coupling == "ov-imp" && return """
          Schwarz impedance overlap:
            - side set: $side_set
              source: $partner
              source block: $partner_block
              robin parameter: 0.0
        """
    c.coupling == "no-dn" && return """
          Schwarz DN nonoverlap:
            - side set: $side_set
              source: $partner
              source side set: $partner_side_set
              default BC type: $(domain == support_domain(c.problem) ? "Neumann" : "Dirichlet")
        """
    return """
          Schwarz impedance nonoverlap:
            - side set: $side_set
              source: $partner
              source side set: $partner_side_set
              robin parameter: $(num(problem_parameters(c.problem).robin_parameter))
        """
end

function subdomain_input(c::Case, domain, kind::Char, time_step)
    text = "type: single\ninput mesh file: $domain.g\noutput mesh file: $domain.e\n"
    text *= materials(c.problem, domain)
    text *= time_integrator(kind, time_step)
    text *= initial_conditions(c.problem)
    conditions = ""
    domain in ("beam", "cylinders", support_domain(c.problem)) && (conditions *= support(c.problem))
    c.coupling != "mono" && (conditions *= coupling_condition(c, domain))
    isempty(conditions) || (text *= "boundary conditions:\n" * conditions)
    return text * solver(kind, c.problem)
end

# Relaxation of the Schwarz iteration: none for the overlap couplings,
# recursive Aitken for Dirichlet-Neumann, and for the paired impedance a
# fixed factor of 0.5 whenever a subdomain is explicit (Aitken diverges on
# the explicit cylinders) and recursive Aitken on the implicit beam, where
# it needs a tenth of the iterations (see docs/notes/schwarz-coupling).
function relaxation(c::Case)
    c.coupling == "no-dn" && return "relaxation: aitken recursive\n"
    c.coupling == "no-imp" && return (c.pair == "II" && c.problem == "beam") ? "relaxation: aitken recursive\n" :
                                     "relaxation parameter: 0.5\n"
    return ""
end

function run_input(c::Case, domains, time_step, final_time, exodus_interval)
    p = problem_parameters(c.problem)
    list = join(("\"$d.yaml\"" for d in domains), ", ")
    return """
        type: multi
        domains: [$list]
        blended energy output: true
        Exodus output interval: $(num(exodus_interval))
        CSV output interval: 0
        initial time: 0.0
        final time: $(num(final_time))
        time step: $(num(time_step))
        minimum iterations: 1
        maximum iterations: $(p.maximum_iterations)
        relative tolerance: $(num(p.relative_tolerance))
        absolute tolerance: $(num(p.absolute_tolerance))
        """ * relaxation(c)
end

function write_case(c::Case, options)
    dir = joinpath(options["runs"], case_name(c))
    isdir(dir) && rm(dir; recursive=true)
    mkpath(dir)
    p = problem_parameters(c.problem)
    time_step = p.time_step / c.level
    final_time = something(options["final time"], p.final_time)
    if c.coupling == "mono"
        domains = [c.problem == "beam" ? "beam" : "cylinders"]
        kinds = [c.pair[1]]
    else
        # The subdomain without the support is listed first, so the Schwarz
        # relaxation acts on the support side, as in the repository examples.
        domains = [other_domain(c.problem), support_domain(c.problem)]
        kinds = [c.pair[2], c.pair[1]]
    end
    for (domain, kind) in zip(domains, kinds)
        write(joinpath(dir, "$domain.yaml"), subdomain_input(c, domain, kind, time_step))
    end
    write(joinpath(dir, "run.yaml"), run_input(c, domains, time_step, final_time, options["exodus interval"]))
    if c.problem == "beam"
        write_beam_meshes(c, dir)
    else
        link_cylinder_meshes(c, dir, joinpath(dirname(options["runs"]), "meshes"))
    end
    write(joinpath(dir, "case.txt"), """
        name: $(case_name(c))
        problem: $(c.problem)
        coupling: $(c.coupling)
        pair: $(c.pair)
        level: $(c.level)
        ratio: $(c.ratio)
        time step: $(num(time_step))
        final time: $(num(final_time))
        """)
    return dir
end

function main(args)
    selections, options = parse_arguments(args)
    cases = selected_cases(selections)
    if options["list"]
        foreach(c -> println(case_name(c)), cases)
        println("$(length(cases)) cases")
        return
    end
    for c in cases
        println("  ", basename(write_case(c, options)))
    end
    println("$(length(cases)) cases written under $(options["runs"])")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
