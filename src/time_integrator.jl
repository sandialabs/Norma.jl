# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

function adaptive_stepping_parameters(integrator_params::Parameters)
    has_minimum = haskey(integrator_params, "minimum time step")
    has_decrease = haskey(integrator_params, "decrease factor")
    has_maximum = haskey(integrator_params, "maximum time step")
    has_increase = haskey(integrator_params, "increase factor")
    has_any = has_minimum || has_decrease || has_maximum || has_increase
    has_all = has_minimum && has_decrease && has_maximum && has_increase
    if has_any == true && has_all == false
        norma_abort(
            "Adaptive time stepping requires 4 parameters: minimum and maximum time steps and decrease and increase factors",
        )
    elseif has_any == true && has_all == true
        minimum_time_step = integrator_params["minimum time step"]
        decrease_factor = integrator_params["decrease factor"]
        maximum_time_step = integrator_params["maximum time step"]
        increase_factor = integrator_params["increase factor"]
    else
        minimum_time_step = maximum_time_step = integrator_params["time step"]
        decrease_factor = increase_factor = 1.0
    end
    if minimum_time_step > maximum_time_step
        norma_abortf(
            "Minimum time step %.4e must be less than or equal to maximum %.4e.", minimum_time_step, maximum_time_step
        )
    end
    if decrease_factor > 1.0
        norma_abortf("Decrease factor %.4e must be less than or equal to one.", decrease_factor)
    end
    if increase_factor < 1.0
        norma_abortf("Increase factor %.4e must be greater than or equal to one.", increase_factor)
    end
    return minimum_time_step, decrease_factor, maximum_time_step, increase_factor
end

function QuasiStatic(params::Parameters, model::Model)
    integrator_params = params["time integrator"]
    time_step = integrator_params["time step"]
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(
        integrator_params
    )
    time = prev_time = -Inf
    num_dof = length(model.free_dofs)
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    prev_disp = Float64[]
    prev_velo = Float64[]
    prev_acce = Float64[]
    prev_∂Ω_f = Float64[]
    stored_energy = 0.0
    initial_equilibrium = get(integrator_params, "initial equilibrium", false)
    return QuasiStatic(
        prev_time,
        time,
        time_step,
        minimum_time_step,
        maximum_time_step,
        decrease_factor,
        increase_factor,
        displacement,
        velocity,
        acceleration,
        prev_disp,
        prev_velo,
        prev_acce,
        prev_∂Ω_f,
        stored_energy,
        initial_equilibrium,
    )
end

function Newmark(params::Parameters, model::Model)
    integrator_params = params["time integrator"]
    time_step = integrator_params["time step"]
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(
        integrator_params
    )
    time = prev_time = -Inf
    β = integrator_params["β"]
    γ = integrator_params["γ"]
    hht_alpha = Float64(get(integrator_params, "HHT α", get(integrator_params, "HHT alpha", 0.0)))
    if hht_alpha < 0.0 || hht_alpha > 1.0 / 3.0
        norma_abort("`HHT alpha` must be in [0, 1/3]; got $hht_alpha.")
    end
    if hht_alpha > 0.0
        γ = 0.5 + hht_alpha
        β = 0.25 * (1.0 + hht_alpha)^2
        norma_logf(
            0,
            :setup,
            "HHT-α dissipation: α = %.4g, overriding γ = %.4g, β = %.4g (ρ∞ = %.3f)",
            hht_alpha, γ, β, (1.0 - hht_alpha) / (1.0 + hht_alpha),
        )
    end
    num_dof = length(model.free_dofs)
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    prev_disp = Float64[]
    prev_velo = Float64[]
    prev_acce = Float64[]
    prev_∂Ω_f = Float64[]
    disp_pre = zeros(num_dof)
    velo_pre = zeros(num_dof)
    stored_energy = 0.0
    kinetic_energy = 0.0
    return Newmark(
        prev_time,
        time,
        time_step,
        minimum_time_step,
        maximum_time_step,
        decrease_factor,
        increase_factor,
        β,
        γ,
        hht_alpha,
        zeros(num_dof),
        displacement,
        velocity,
        acceleration,
        disp_pre,
        velo_pre,
        prev_disp,
        prev_velo,
        prev_acce,
        prev_∂Ω_f,
        stored_energy,
        kinetic_energy,
    )
end

function CentralDifference(params::Parameters, model::Model)
    integrator_params = params["time integrator"]
    time_step = integrator_params["time step"]
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(
        integrator_params
    )
    time = prev_time = -Inf
    CFL = integrator_params["CFL"]
    γ = integrator_params["γ"]
    num_dof = length(model.free_dofs)
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    prev_disp = Float64[]
    prev_velo = Float64[]
    prev_acce = Float64[]
    prev_∂Ω_f = Float64[]
    stored_energy = 0.0
    kinetic_energy = 0.0
    return CentralDifference(
        prev_time,
        time,
        time_step,
        minimum_time_step,
        maximum_time_step,
        decrease_factor,
        increase_factor,
        CFL,
        γ,
        displacement,
        velocity,
        acceleration,
        prev_disp,
        prev_velo,
        prev_acce,
        prev_∂Ω_f,
        stored_energy,
        kinetic_energy,
    )
end

function create_time_integrator(params::Parameters, model::SolidMechanics)
    integrator_params = params["time integrator"]
    integrator_name = integrator_params["type"]
    if integrator_name == "quasi static"
        return QuasiStatic(params, model)
    elseif integrator_name == "Newmark"
        return Newmark(params, model)
    elseif integrator_name == "central difference"
        return CentralDifference(params, model)
    else
        norma_abort("Unknown type of time integrator : $integrator_name")
    end
end

function is_static(integrator::TimeIntegrator)
    return integrator isa QuasiStatic
end

function is_dynamic(integrator::TimeIntegrator)
    return is_static(integrator) == false
end

function initialize(integrator::QuasiStatic, solver::Solver, model::SolidMechanics; trust_schwarz::Bool=false)
    if integrator.initial_equilibrium == true
        norma_log(0, :equilibrium, "Establishing Initial Equilibrium...")
        solve(integrator, solver, model)
        if model.failed == true
            norma_abort("Finite element model failed to establish initial equlibrium")
        end
    end
    return nothing
end

function predict(integrator::QuasiStatic, solver::Solver, model::SolidMechanics)
    return nothing
end

function correct(integrator::QuasiStatic, solver::Solver, model::SolidMechanics)
    return nothing
end

function initialize(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics; trust_schwarz::Bool=false)
    norma_log(0, :acceleration, "Computing Initial Acceleration...")
    free = model.free_dofs
    # HHT-α: make the blend a no-op until the first predictor runs.
    integrator.hht_alpha > 0.0 && (integrator.hht_disp_prev .= integrator.displacement)
    evaluate(model, integrator, solver)
    if model.failed == true
        norma_abort("Finite element model failed to initialize")
    end
    # The Robin and impedance-Schwarz self terms belong to the internal force
    # (they are added there by the Newton assembly in solver.jl), while the
    # partner terms already sit in boundary_force from apply_bcs. Omitting the
    # self terms here leaves an unbalanced interface force ∝ α in the initial
    # acceleration: a one-step impulsive kick, visible from any start whose
    # interface displacement or velocity is nonzero (e.g. a stress-free
    # finite-rotation initial condition).
    apply_robin_bcs_internal_force!(model)
    internal_force = model.internal_force + build_impedance_schwarz_force(model)
    external_force = model.body_force + model.boundary_force
    inertial_force = external_force - internal_force
    kinetic_energy = 0.5 * dot(integrator.velocity, model.mass, integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    integrator.stored_energy = model.strain_energy
    atol = solver.linear_solver_absolute_tolerance
    rtol = solver.linear_solver_relative_tolerance
    # Mff * a_free = (F_ext - F_int)_free - Mfc * a_fixed. The prescribed
    # (Dirichlet) DOFs generally have a nonzero acceleration of their own
    # (e.g. a moving boundary condition), already set on model/integrator
    # .acceleration by apply_bcs. With a consistent (non-lumped) mass matrix,
    # free and fixed DOFs adjacent to the same elements are inertially
    # coupled through the Mfc off-diagonal block, so that contribution has to
    # be subtracted from the right-hand side. This correction is applied
    # unconditionally, including on a fresh (non-restarted) run: when a_fixed
    # is exactly zero (e.g. a run starting at rest), the Mfc*a_fixed term
    # vanishes and this is a no-op, but a fresh run with a nonzero prescribed
    # boundary acceleration at t=0 needs it just as much as a restart does.
    fixed = .!free
    # Schwarz-coupled fixed DOFs are excluded from this correction by
    # default (trust_schwarz=false): this function runs once, at the very
    # start of the simulation (or of a restart), before any subdomain has
    # ever taken a predictor step or had its own initial acceleration
    # computed. The Schwarz coupling machinery (coupling_weak_dbc /
    # contact_weak_dbc / overlap variants) derives its acceleration either
    # by reading the coupled side's integrator.acceleration directly, or —
    # under Aitken-secant relaxation — via integrator.disp_pre, both of
    # which are only physically meaningful once the coupled side has a real
    # acceleration of its own. The first time this function runs for any
    # subdomain, apply_bcs() has only ever seen every subdomain's untouched
    # (zero) initial acceleration, so whatever value it left on
    # Schwarz-coupled DOFs is not trustworthy and must not be allowed into
    # this one-shot linear solve. A plain (non-Schwarz) Dirichlet BC's
    # prescribed acceleration is exact at any time, including t=0, and is
    # always trusted.
    #
    # trust_schwarz=true is used for a second, refinement pass over all
    # subdomains of a restarted multi-domain (Schwarz) simulation — see
    # initialize(sim::MultiDomainSimulation) in simulation.jl. By the time
    # that second pass runs, every subdomain already has a real initial
    # acceleration from the first pass, apply_bcs() has been re-run so
    # Schwarz-coupled DOFs now carry a real (not stale-zero) coupled
    # acceleration, and it is safe — indeed necessary for an accurate
    # restart, since the true acceleration at a Schwarz interface generally
    # is not zero mid-simulation — to include them in the correction below.
    # This is a single Jacobi-style correction, not a fully converged fixed
    # point; any remaining small inconsistency is absorbed by the real
    # Newton/Schwarz iterations that run on the first actual control step.
    schwarz_fixed = falses(length(free))
    for bc in model.boundary_conditions
        if bc isa SolidMechanicsSchwarzBoundaryCondition
            for i_global in bc.global_from_local_map
                schwarz_fixed[(3 * i_global - 2):(3 * i_global)] .= true
            end
        end
    end
    trusted_fixed = trust_schwarz ? fixed : (fixed .& .!schwarz_fixed)
    # Read the prescribed (Dirichlet / Schwarz-coupled) accelerations from
    # model.acceleration rather than integrator.acceleration. apply_bcs()
    # always writes prescribed values onto model.acceleration. For a
    # top-level SolidMechanics simulation, integrator.acceleration is the
    # same aliased memory (initialize_storage, simulation.jl), so this is
    # unchanged there. But this function also runs, via the RomModel
    # dispatch in opinf_time_integrator.jl, with model set to a ROM
    # subdomain's fom_model and integrator set to its fom_integrator; that
    # fom_integrator.acceleration is a separate array that is never
    # aliased to fom_model.acceleration, so it stays zero on fixed DOFs
    # even after the Schwarz BCs have written real prescribed values into
    # fom_model.acceleration. Reading model.acceleration (flattened to
    # match the node/component ordering used everywhere else) is
    # therefore correct in both cases.
    prescribed_acceleration = vec(model.acceleration)
    rhs_free = inertial_force[free] - model.mass[free, trusted_fixed] * prescribed_acceleration[trusted_fixed]
    integrator.acceleration[free] = solve_linear(model.mass[free, free], rhs_free, atol, rtol)
    return nothing
end

function predict(integrator::Newmark, solver::Solver, model::SolidMechanics)
    free = model.free_dofs
    fixed = .!free
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    u = integrator.displacement
    # HHT-α: the internal force is evaluated at (1-ᾱ) u_{n+1} + ᾱ u_n;
    # capture u_n before the predictor overwrites it.
    integrator.hht_alpha > 0.0 && (integrator.hht_disp_prev .= u)
    v = integrator.velocity
    a = integrator.acceleration
    u_pre = integrator.disp_pre
    v_pre = integrator.velo_pre
    u_pre[free] = u[free] += Δt * v[free] + (0.5 - β) * Δt * Δt * a[free]
    v_pre[free] = v[free] += (1.0 - γ) * Δt * a[free]
    u_pre[fixed] = u[fixed]
    v_pre[fixed] = v[fixed]
    # Set acceleration consistent with predicted displacement (u = u_pre → a = 0)
    a[free] .= 0.0
    return nothing
end

function correct(integrator::Newmark, solver::Solver, model::SolidMechanics)
    free = model.free_dofs
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    integrator.displacement = solver.solution
    u = integrator.displacement
    u_pre = integrator.disp_pre
    v_pre = integrator.velo_pre
    integrator.acceleration[free] = (u[free] - u_pre[free]) / β / Δt / Δt
    integrator.velocity[free] = v_pre[free] + γ * Δt * integrator.acceleration[free]
    return nothing
end

function initialize(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics; trust_schwarz::Bool=false)
    norma_log(0, :acceleration, "Computing Initial Acceleration...")
    free = model.free_dofs
    set_time_step(integrator, model)
    evaluate(model, integrator, solver)
    if model.failed == true
        norma_abort("The finite element model has failed to initialize")
    end
    # Same interface-force balance as the Newmark initializer above: the Robin
    # and impedance-Schwarz self terms must accompany the partner terms that
    # apply_bcs left in boundary_force.
    apply_robin_bcs_internal_force!(model)
    internal_force = model.internal_force + build_impedance_schwarz_force(model)
    external_force = model.body_force + model.boundary_force
    kinetic_energy = 0.5 * model.lumped_mass ⋅ (integrator.velocity .* integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    integrator.stored_energy = model.strain_energy
    inertial_force = external_force - internal_force
    integrator.acceleration[free] = inertial_force[free] ./ model.lumped_mass[free]
    return nothing
end

function predict(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    free = model.free_dofs
    set_time_step(integrator, model)
    Δt = integrator.time_step
    γ = integrator.γ
    u = integrator.displacement
    v = integrator.velocity
    a = integrator.acceleration
    u[free] += Δt * v[free] + 0.5 * Δt * Δt * a[free]
    v[free] += (1.0 - γ) * Δt * a[free]
    return nothing
end

function correct(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    Δt = integrator.time_step
    γ = integrator.γ
    integrator.acceleration = solver.solution
    a = integrator.acceleration
    free = model.free_dofs
    integrator.velocity[free] += γ * Δt * a[free]
    return nothing
end

include("opinf/opinf_time_integrator.jl")
