# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
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
    predictor_name = String(get(integrator_params, "predictor", "constant"))
    if predictor_name ∉ ("constant", "linear", "tangent")
        norma_abortf("Unknown `predictor: %s`. Valid values are `constant`, `linear` and `tangent`.", predictor_name)
    end
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
        Symbol(predictor_name),
        Float64[],
        false,
        Float64[],
        0.0,
        Float64[],
        0.0,
        spzeros(0, 0),
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
            hht_alpha,
            γ,
            β,
            (1.0 - hht_alpha) / (1.0 + hht_alpha),
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
    stable_name = get(integrator_params, "stable time step", "length")
    stable_method = if stable_name == "length"
        :length
    elseif stable_name == "eigenvalue"
        :eigenvalue
    else
        norma_abort("Unknown `stable time step: $(stable_name)`. Valid values are `length` (default) and `eigenvalue`.")
    end
    stable_interval = Int(get(integrator_params, "stable time step interval", stable_method === :length ? 1 : 100))
    stable_interval ≥ 1 || norma_abort("`stable time step interval` must be at least 1, got $(stable_interval).")
    if model isa SolidMechanics && any(block -> block.element_type == TETRA10, model.blocks)
        norma_log(
            0,
            :warning,
            "Central difference on TETRA10 elements: their row-sum lumped masses are negative at the corners, so the " *
            "estimate of the critical step uses the scaled diagonal of the consistent mass and is not a bound. " *
            "The composite TETRA10 (`composite tet 10`) has positive lumped masses.",
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
        stable_method,
        stable_interval,
        NaN,
        0,
        Float64[],
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

function initialize(integrator::QuasiStatic, solver::Solver, model::SolidMechanics; include_schwarz_rows::Bool=false)
    if integrator.initial_equilibrium == true
        norma_log(0, :equilibrium, "Establishing Initial Equilibrium...")
        solve(integrator, solver, model)
        if model.failed == true
            norma_abort("Finite element model failed to establish initial equlibrium")
        end
    else
        # A periodic cell starts from the affine field of F̄ at the initial
        # time, which the initial equilibrium would predict. From the corners
        # alone, `tangent` would leave F̄ - I on the elements next to them at
        # the first step, and `linear` would extrapolate it into the second.
        affine = periodic_affine_increment!(model)
        if affine !== nothing
            integrator.displacement .+= affine
            update_periodic_dependents!(model)
        end
    end
    return nothing
end

# The quasi-static predictor. A step starts from the previous solution with the
# new values of the prescribed DOFs, so without a predictor (`constant`) the
# first Newton iterate puts the whole increment of a loaded face on the layer
# of elements next to it. A path-dependent material evaluated at that strain
# can fail where the converged step would not: the flow rule of crystal
# plasticity has a bounded range. `linear` adds to the free DOFs the increment
# of the last accepted step, scaled by the ratio of the time steps; the first
# step has none to extrapolate. `tangent` adds the solution of
#
#     K_ff Δu_f = -K_fc Δu_c,
#
# the linear response to the increment Δu_c of the prescribed DOFs, with the
# stiffness K of the last accepted step, or before the first one with the
# stiffness at the start of the step. A predictor changes the iterate Newton
# starts from, not the solution it converges to.
#
# The increment of the prescribed DOFs needs the displacement the step starts
# from, which start_step! records before apply_bcs: the last accepted
# solution, the start of the step after a failed attempt, or the start of the
# stop in a Schwarz iteration.
start_step!(::TimeIntegrator) = nothing
accept_step!(::TimeIntegrator, ::Solver) = nothing

function start_step!(integrator::QuasiStatic)
    integrator.predictor == :constant && return nothing
    integrator.step_start_disp = copy(integrator.displacement)
    integrator.step_started = true
    return nothing
end

function accept_step!(integrator::QuasiStatic, solver::Solver)
    integrator.predictor == :constant && return nothing
    length(integrator.step_start_disp) == length(integrator.displacement) || return nothing
    integrator.last_increment = integrator.displacement - integrator.step_start_disp
    integrator.last_time_step = integrator.time - integrator.prev_time
    if integrator.predictor == :tangent
        # The Hessian of a quasi-static step is the stiffness of the model,
        # whose values the next evaluation rewrites; the predictor of the next
        # step runs before that evaluation, so the reference holds the
        # converged stiffness of this step when it is read.
        integrator.last_stiffness = solver.hessian
    end
    return nothing
end

# The history starts again when a simulation takes over the state of another
# (transfer_state!), as it starts with a run: a variant that an in-memory
# restart switches back in holds the history of its last step before it was
# switched out, not of the step before the switch.
reset_predictor!(::TimeIntegrator) = nothing

function reset_predictor!(integrator::QuasiStatic)
    integrator.step_started = false
    integrator.last_increment = Float64[]
    integrator.last_time_step = 0.0
    integrator.last_stiffness = spzeros(0, 0)
    return nothing
end

# The tangent predictor solves with the stiffness, which only the Hessian
# minimizer assembles.
function validate_predictor(integrator::QuasiStatic, solver::Solver)
    if integrator.predictor == :tangent && !(solver isa HessianMinimizer)
        norma_abort("`predictor: tangent` requires the `Hessian minimizer` solver.")
    end
    return nothing
end

validate_predictor(::TimeIntegrator, ::Solver) = nothing

function predict(integrator::QuasiStatic, solver::Solver, model::SolidMechanics)
    empty!(integrator.predicted_increment)
    # On a periodic cell, the increment of F̄ applied to the whole cell. It is
    # the prediction of `constant`, and of a step the predictor does not
    # predict, such as the first step of `linear`. A step the predictor
    # predicts gets the whole field from it, so the increment is not added
    # twice.
    affine = periodic_affine_increment!(model)
    # Only the solve of a step uses the predictor; the initial equilibrium
    # solve has no step to start from.
    started = integrator.step_started
    integrator.step_started = false
    δ = if started && integrator.predictor == :linear
        linear_prediction(integrator, model)
    elseif started && integrator.predictor == :tangent
        tangent_prediction(integrator, solver, model)
    else
        nothing
    end
    δ === nothing && (δ = affine)
    δ === nothing && return nothing
    integrator.unpredicted_residual = evaluate_unpredicted_residual(integrator, solver, model)
    integrator.displacement .+= δ
    integrator.predicted_increment = δ
    return nothing
end

# The predicted increment of the free DOFs, zero on the prescribed ones, or
# nothing when there is no prediction.
function linear_prediction(integrator::QuasiStatic, model::SolidMechanics)
    increment = integrator.last_increment
    integrator.last_time_step > 0.0 && length(increment) == length(model.free_dofs) || return nothing
    ratio = (integrator.time - integrator.prev_time) / integrator.last_time_step
    return ratio .* increment .* model.free_dofs
end

function tangent_prediction(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    u = integrator.displacement
    u_start = integrator.step_start_disp
    free = model.free_dofs
    fixed = .!free
    Δu_fixed = u[fixed] - u_start[fixed]
    any(!iszero, Δu_fixed) || return nothing
    K = integrator.last_stiffness
    if size(K, 1) != length(u)
        u_fixed = u[fixed]
        u[fixed] = u_start[fixed]
        evaluate(integrator, solver, model)
        u[fixed] = u_fixed
        if model.failed == true
            norma_log(4, :solve, "Tangent predictor: the start of the step failed to evaluate; no prediction")
            model.failed = false
            return nothing
        end
        K = solver.hessian
    end
    δ = zeros(length(u))
    K_free = reduced_matrix!(solver, K, free)
    δ[free] = solve_linear(solver, K_free, -(K[free, fixed] * Δu_fixed), symmetric_tangent(model))
    return δ
end

# The relative tolerance of a solve is measured from its initial residual,
# which an accurate predictor makes small: for a linear problem the predicted
# iterate is the solution, and measured from its residual the relative
# tolerance asks for less than the roundoff. The solve measures it instead from
# the residual the iterate would have without the predictor, if larger. The
# Hessian minimizer estimates that residual as R - K δ, to first order, with δ
# the predicted increment. A matrix-free solver has no K, so it evaluates the
# residual of the unpredicted iterate before the prediction is added: one more
# evaluation per predicted step. If that evaluation fails, the tolerance is
# measured from the predicted iterate.
evaluate_unpredicted_residual(::QuasiStatic, ::Solver, ::SolidMechanics) = 0.0

function evaluate_unpredicted_residual(integrator::QuasiStatic, solver::MatrixFree, model::SolidMechanics)
    # The evaluation brings the nodes of an exact Surface condition onto their
    # surfaces; the iterate stays as it was.
    u = copy(integrator.displacement)
    evaluate(integrator, solver, model)
    integrator.displacement .= u
    residual = model.failed ? 0.0 : norm(solver.gradient[model.free_dofs])
    model.failed = false
    return residual
end

unpredicted_residual_norm(::TimeIntegrator, ::Solver, ::Model) = 0.0

function unpredicted_residual_norm(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    δ = integrator.predicted_increment
    length(δ) == length(model.free_dofs) || return 0.0
    return norm((solver.gradient - solver.hessian * δ)[model.free_dofs])
end

function unpredicted_residual_norm(integrator::QuasiStatic, ::MatrixFree, model::SolidMechanics)
    length(integrator.predicted_increment) == length(model.free_dofs) || return 0.0
    return integrator.unpredicted_residual
end

function correct(integrator::QuasiStatic, solver::Solver, model::SolidMechanics)
    return nothing
end

function initialize(
    integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics; include_schwarz_rows::Bool=false
)
    norma_log(0, :acceleration, "Computing Initial Acceleration...")
    free = model.free_dofs
    # HHT-α: make the blend a no-op until the first predictor runs.
    integrator.hht_alpha > 0.0 && (integrator.hht_disp_prev .= integrator.displacement)
    evaluate(model, integrator, solver)
    if model.failed == true
        norma_abort("Finite element model failed to initialize")
    end
    # The Robin and Robin-Robin Schwarz self terms belong to the internal force
    # (they are added there by the Newton assembly in solver.jl), while the
    # partner terms already sit in boundary_force from apply_bcs. Omitting the
    # self terms here leaves an unbalanced interface force ∝ α in the initial
    # acceleration: a one-step impulsive kick, visible from any start whose
    # interface displacement or velocity is nonzero (e.g. a stress-free
    # finite-rotation initial condition).
    apply_robin_bcs_internal_force!(model)
    internal_force = model.internal_force + build_robin_schwarz_force(model)
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
    # The rows held by a Schwarz coupling are prescribed rows whose value is
    # the partner's projected acceleration. They are excluded from the
    # correction by default (include_schwarz_rows=false): this function runs
    # once at the start of the simulation (or of a restart), before any
    # subdomain has computed its own initial acceleration, so apply_bcs()
    # has written the partner's zero acceleration onto those rows, and
    # including them would subtract a zero datum that is not the partner's
    # acceleration. The acceleration prescribed by a plain Dirichlet
    # boundary condition is exact at any time, including t = 0, and is
    # always included.
    #
    # include_schwarz_rows=true is used for the refinement passes over all
    # subdomains of a multidomain simulation (initialize(sim::MultiDomainSimulation)
    # and coupled_initial_acceleration! in simulation.jl). By then every
    # subdomain has an initial acceleration from the first pass and
    # apply_bcs() has written the partner's acceleration onto the coupling
    # rows, so the rows carry the datum the correction needs; on a restart
    # the acceleration at a Schwarz interface is in general not zero, and
    # the pass is required for a correct restart. One pass is a single
    # Jacobi correction, not a converged fixed point; the remaining
    # inconsistency is removed by the Newton and Schwarz iterations of the
    # first stop.
    schwarz_fixed = falses(length(free))
    for bc in model.boundary_conditions
        if bc isa SolidMechanicsSchwarzBoundaryCondition
            for i_global in bc.global_from_local_map
                schwarz_fixed[(3 * i_global - 2):(3 * i_global)] .= true
            end
        end
    end
    correction_rows = include_schwarz_rows ? fixed : (fixed .& .!schwarz_fixed)
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
    rhs_free = inertial_force[free] - model.mass[free, correction_rows] * prescribed_acceleration[correction_rows]
    # With `linear solver: direct` the mass matrix is factored, as the Newton
    # systems are. Conjugate gradients at the linear solver tolerances leaves a
    # relative error of about 1e-10 in the acceleration, which the coupled
    # initial acceleration of a matching pair repeats at every iteration and
    # cannot reduce below (a floor of 5e-10 on the beam against a tolerance of
    # 1e-12). The consistent mass is not always positive definite: the
    # default 5-point rule of the TETRA10 has a negative weight. Cholesky then
    # fails, and LU factors it, as solve_linear_direct does.
    mass_free = reduced_matrix!(solver, model.mass, free)
    integrator.acceleration[free] = if solver.linear_solver == "direct"
        factor = cholesky(Symmetric(mass_free); perm=fill_reducing_ordering(mass_free), check=false)
        issuccess(factor) ? factor \ rhs_free : lu(mass_free) \ rhs_free
    else
        solve_linear(mass_free, rhs_free, atol, rtol)
    end
    return nothing
end

function predict(integrator::Newmark, solver::Solver, model::SolidMechanics)
    free = model.free_dofs
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
    # One pass over the degrees of freedom: the masked forms u[free] copy the
    # vector at every use. The acceleration of the free degrees of freedom is
    # set consistent with the predicted displacement (u = u_pre → a = 0).
    @inbounds for i in eachindex(u)
        if free[i]
            u[i] += Δt * v[i] + (0.5 - β) * Δt * Δt * a[i]
            v[i] += (1.0 - γ) * Δt * a[i]
            a[i] = 0.0
        end
        u_pre[i] = u[i]
        v_pre[i] = v[i]
    end
    return nothing
end

function correct(integrator::Newmark, solver::Solver, model::SolidMechanics)
    free = model.free_dofs
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    integrator.displacement = solver.solution
    u = integrator.displacement
    v = integrator.velocity
    a = integrator.acceleration
    u_pre = integrator.disp_pre
    v_pre = integrator.velo_pre
    @inbounds for i in eachindex(u)
        free[i] || continue
        a[i] = (u[i] - u_pre[i]) / β / Δt / Δt
        v[i] = v_pre[i] + γ * Δt * a[i]
    end
    return nothing
end

function initialize(
    integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics; include_schwarz_rows::Bool=false
)
    norma_log(0, :acceleration, "Computing Initial Acceleration...")
    free = model.free_dofs
    set_time_step(integrator, model)
    evaluate(model, integrator, solver)
    if model.failed == true
        norma_abort("The finite element model has failed to initialize")
    end
    # Same interface-force balance as the Newmark initializer above: the Robin
    # and Robin-Robin Schwarz self terms must accompany the partner terms that
    # apply_bcs left in boundary_force.
    apply_robin_bcs_internal_force!(model)
    internal_force = model.internal_force + build_robin_schwarz_force(model)
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
    # advance_time already capped the step at the stable step and advanced
    # the time by that same step.
    Δt = integrator.time_step
    γ = integrator.γ
    u = integrator.displacement
    v = integrator.velocity
    a = integrator.acceleration
    @inbounds for i in eachindex(u)
        free[i] || continue
        u[i] += Δt * v[i] + 0.5 * Δt * Δt * a[i]
        v[i] += (1.0 - γ) * Δt * a[i]
    end
    return nothing
end

function correct(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    Δt = integrator.time_step
    γ = integrator.γ
    integrator.acceleration = solver.solution
    a = integrator.acceleration
    v = integrator.velocity
    free = model.free_dofs
    @inbounds for i in eachindex(v)
        free[i] && (v[i] += γ * Δt * a[i])
    end
    return nothing
end

include("opinf/opinf_time_integrator.jl")
