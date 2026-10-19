# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Matching Dirichlet-Neumann exchange (`matching: velocity`) on the
# conforming cantilever of examples/nonoverlap/dynamic-same-step/cantilever-dn.
# The free subdomain is the Dirichlet side. Checked: the force transfer of the
# Neumann side is the transpose of the Dirichlet projector; at convergence the
# force applied on the Neumann side is the transposed projection of the
# Dirichlet side's d'Alembert reaction; the Dirichlet side's interface
# acceleration and displacement satisfy its own Newmark relations over the
# last step; and the rate invariant Ẽ is conserved over 50 steps.
# The initial acceleration of a constrained pair is found by a Schwarz
# iteration at t = 0 (coupled_initial_acceleration!), so conservation holds
# from step 0; the last two testsets check the first step and the absence of
# the alternating interface acceleration difference that the per-subdomain
# initial acceleration left in the explicit pair.

using YAML
using LinearAlgebra
using Exodus

const constrained_dn_example = "../examples/nonoverlap/dynamic-same-step/cantilever-dn"

function constrained_dn_subdomain(
    name::String,
    explicit::Bool,
    constraint::String,
    dt::Float64,
    own_dt::Float64=dt;
    swap_roles::Bool=false,
    interface_solve::String="iterative",
    edit_bc=nothing,
    edit_sub=nothing,
)
    sub = YAML.load_file("$constrained_dn_example/$name.yaml"; dicttype=Norma.Parameters)
    edit_sub === nothing || edit_sub(name, sub)
    if explicit
        sub["time integrator"] = Norma.Parameters(
            "type" => "central difference", "time step" => own_dt, "CFL" => 1.0, "γ" => 0.5
        )
        sub["solver"] = Norma.Parameters("type" => "explicit solver", "step" => "explicit")
    else
        sub["time integrator"]["time step"] = own_dt
        sub["solver"]["linear solver"] = "direct"
        sub["solver"]["linear solver relative tolerance"] = 1.0e-14
        sub["solver"]["maximum iterations"] = 4
        sub["solver"]["absolute tolerance"] = 1.0e-6
    end
    bc = sub["boundary conditions"]["Schwarz DN nonoverlap"][1]
    if swap_roles
        bc["default BC type"] = bc["default BC type"] == "Dirichlet" ? "Neumann" : "Dirichlet"
    end
    bc["matching"] = constraint
    bc["interface solve"] = interface_solve
    edit_bc === nothing || edit_bc(bc)
    YAML.write_file("$name.yaml", sub)
    return nothing
end

function run_constrained_dn(
    explicit::Bool;
    constraint="velocity",
    num_steps=50,
    dt=1.0e-6,
    theta=0.5,
    free_substeps=1,
    clamped_dirichlet=false,
    interface_solve="iterative",
    mesh_dir=constrained_dn_example,
    initialize_only=false,
    relaxation=nothing,
    edit_bc=nothing,
    edit_params=nothing,
    clamped_mesh="cantilever-clamped.g",
    free_mesh="cantilever-free.g",
    edit_sub=nothing,
)
    cp("$mesh_dir/$clamped_mesh", "cantilever-clamped.g"; force=true)
    cp("$mesh_dir/$free_mesh", "cantilever-free.g"; force=true)
    free_dt = dt / free_substeps
    constrained_dn_subdomain(
        "cantilever-free",
        explicit,
        constraint,
        dt,
        free_dt;
        swap_roles=clamped_dirichlet,
        interface_solve,
        edit_bc,
        edit_sub,
    )
    constrained_dn_subdomain(
        "cantilever-clamped", explicit, constraint, dt; swap_roles=clamped_dirichlet, interface_solve, edit_bc, edit_sub
    )
    params = Norma.Parameters(
        "type" => "multi",
        "name" => "constrained-dn",
        "domains" => ["cantilever-free.yaml", "cantilever-clamped.yaml"],
        "Exodus output interval" => 1.0,
        "initial time" => 0.0,
        "final time" => num_steps * dt,
        "time step" => dt,
        "minimum iterations" => 1,
        "maximum iterations" => 64,
        "matching relative tolerance" => 1.0e-12,
        "relaxation parameter" => theta,
        "blended energy output" => true,
    )
    if relaxation !== nothing
        params["relaxation"] = relaxation
        params["aitken N0 parameter"] = 1
    end
    edit_params === nothing || edit_params(params)
    files = [
        "cantilever-clamped.g",
        "cantilever-free.g",
        "cantilever-clamped.e",
        "cantilever-free.e",
        "cantilever-clamped.yaml",
        "cantilever-free.yaml",
        "constrained-dn-energy.csv",
    ]
    if initialize_only
        # The state at t = 0 after the coupled initial acceleration.
        sim = Norma.create_simulation(params)
        Norma.sync_control_time(sim)
        Norma.initialize(sim)
        for subsim in sim.subsims
            Exodus.close(subsim.params["input_mesh"])
            Exodus.close(subsim.params["output_mesh"])
        end
        foreach(f -> rm(f; force=true), files)
        return sim, nothing, nothing
    end
    sim = Norma.run(params)
    rows = readlines("constrained-dn-energy.csv")
    header = split(rows[1], ",")
    data = [parse.(Float64, split(row, ",")) for row in rows[2:end]]
    for f in [
        "cantilever-clamped.g",
        "cantilever-free.g",
        "cantilever-clamped.e",
        "cantilever-free.e",
        "cantilever-clamped.yaml",
        "cantilever-free.yaml",
        "constrained-dn-energy.csv",
    ]
        rm(f; force=true)
    end
    return sim, header, data
end

function dn_bc_of(subsim)
    for bc in subsim.model.boundary_conditions
        bc isa Norma.SolidMechanicsNonOverlapSchwarzBoundaryCondition && return bc
    end
    return nothing
end

function energy_drift(header, data, column)
    i = findfirst(==(column), header)
    e = [row[i] for row in data]
    return maximum(abs.(e ./ e[1] .- 1.0))
end

invariant_drift(header, data) = energy_drift(header, data, "rate_invariant_total")

@testset "Matching DN: Operators, Reaction, Newmark Relations (II)" begin
    sim, header, data = run_constrained_dn(false; num_steps=5)
    @test sim.failed == false
    free, clamped = sim.subsims[1], sim.subsims[2]
    bc_D = dn_bc_of(free)
    bc_N = dn_bc_of(clamped)
    @test bc_D.is_dirichlet && !bc_N.is_dirichlet
    @test bc_D.constraint == :velocity && bc_N.constraint == :velocity
    # One cross matrix: the Neumann force transfer is the transpose of the
    # Dirichlet projector, and W_D Π_D is the cross matrix seen from either side.
    @test bc_N.neumann_projector == transpose(bc_D.dirichlet_projector)
    @test bc_D.square_projector * bc_D.dirichlet_projector ≈ transpose(bc_N.square_projector * bc_N.dirichlet_projector) rtol =
        1.0e-12
    # Force applied on the Neumann side against Π_Dᵀ r_D, r_D the Dirichlet
    # side's reaction -(M a + f_int - f_body - f_boundary) on its interface rows.
    model_D = free.model
    r_global = -(model_D.internal_force + Norma.dalembert_inertia_minus_loads(model_D))
    r_D = reshape(Norma.extract_local_vector(bc_D, r_global, 3), 3, :)
    expected = vec(transpose(transpose(bc_D.dirichlet_projector) * transpose(r_D)))
    @test norm(bc_N.transferred_force - expected) ≤ 1.0e-10 * norm(expected)
    residuals = Norma.dn_interface_residuals(sim)
    @test length(residuals) == 1
    r = residuals[1][2]
    @test r.force_residual ≤ 1.0e-10
    @test r.velocity_jump ≤ 1.0e-10
    # Newmark relations of the Dirichlet side over the last step, from its
    # state at the start of the last stop.
    integrator = free.integrator
    Δt = integrator.time_step
    β, γ = integrator.β, integrator.γ
    u_n = sim.controller.stop_disp[1]
    v_n = sim.controller.stop_velo[1]
    a_n = sim.controller.stop_acce[1]
    for i_global in bc_D.global_from_local_map, comp in 1:3
        k = 3 * (i_global - 1) + comp
        u_pre = u_n[k] + Δt * v_n[k] + (0.5 - β) * Δt^2 * a_n[k]
        v_pre = v_n[k] + (1.0 - γ) * Δt * a_n[k]
        a = integrator.acceleration[k]
        @test integrator.velocity[k] ≈ v_pre + γ * Δt * a rtol = 1.0e-12 atol = 1.0e-12
        @test integrator.displacement[k] ≈ u_pre + β * Δt^2 * a rtol = 1.0e-12 atol = 1.0e-15
    end
end

@testset "Matching DN: Rate Invariant Conservation (II)" begin
    sim, header, data = run_constrained_dn(false)
    @test sim.failed == false
    @test length(data) == 51
    @test invariant_drift(header, data) < 1.0e-10
end

@testset "Matching DN: Rate Invariant Conservation (EE)" begin
    sim, header, data = run_constrained_dn(true)
    @test sim.failed == false
    @test length(data) == 51
    @test invariant_drift(header, data) < 1.0e-10
    # Central difference on the Dirichlet side: β = 0, so the imposed interface
    # displacement is the predictor value.
    integrator = sim.subsims[1].integrator
    bc_D = dn_bc_of(sim.subsims[1])
    Δt = integrator.time_step
    u_n = sim.controller.stop_disp[1]
    v_n = sim.controller.stop_velo[1]
    a_n = sim.controller.stop_acce[1]
    for i_global in bc_D.global_from_local_map, comp in 1:3
        k = 3 * (i_global - 1) + comp
        @test integrator.displacement[k] ≈ u_n[k] + Δt * v_n[k] + 0.5 * Δt^2 * a_n[k] rtol = 1.0e-12 atol = 1.0e-15
        v_pre = v_n[k] + 0.5 * Δt * a_n[k]
        @test integrator.velocity[k] ≈ v_pre + 0.5 * Δt * integrator.acceleration[k] rtol = 1.0e-12 atol = 1.0e-12
    end
end

@testset "Matching DN: Displacement Matching Requires Newmark" begin
    @test_throws Exception run_constrained_dn(true; constraint="displacement", num_steps=1)
    for f in [
        "cantilever-clamped.g",
        "cantilever-free.g",
        "cantilever-clamped.e",
        "cantilever-free.e",
        "cantilever-clamped.yaml",
        "cantilever-free.yaml",
        "constrained-dn-energy.csv",
    ]
        rm(f; force=true)
    end
end

# The retired keys `constrained` and `constraint` abort instead of being
# ignored, which would run the standard exchange; the top-level key
# `matching relative tolerance` is required when a pair sets `matching`; the
# retired top-level key `constraint absolute tolerance` aborts. The controller's
# `relative tolerance` and `absolute tolerance` are required when a coupling
# other than a matching pair exists (here the pair without `matching`, which
# is the standard exchange), and ignored when every coupling is a matching pair.
@testset "Matching DN: Retired Keys and Required Tolerances" begin
    cases = [
        (edit_bc=bc -> (delete!(bc, "matching"); bc["constrained"]=true), edit_params=nothing),
        (edit_bc=bc -> (bc["constraint"] = "velocity"), edit_params=nothing),
        (edit_bc=nothing, edit_params=p -> delete!(p, "matching relative tolerance")),
        (edit_bc=nothing, edit_params=p -> (p["constraint absolute tolerance"] = 1.0e-9)),
        (edit_bc=nothing, edit_params=p -> (p["matching relative tolerance"] = -1.0)),
        (edit_bc=bc -> delete!(bc, "matching"), edit_params=nothing),
    ]
    for case in cases
        @test_throws Exception run_constrained_dn(false; num_steps=1, case.edit_bc, case.edit_params)
    end
    controller_tolerances = p -> (p["relative tolerance"]=1.0; p["absolute tolerance"]=1.0)
    sim, _, _ = run_constrained_dn(false; num_steps=2, edit_params=controller_tolerances)
    @test sim.controller.interface_residual ≤ 1.0e-12
    for f in [
        "cantilever-clamped.g",
        "cantilever-free.g",
        "cantilever-clamped.e",
        "cantilever-free.e",
        "cantilever-clamped.yaml",
        "cantilever-free.yaml",
        "constrained-dn-energy.csv",
    ]
        rm(f; force=true)
    end
end

# Displacement matching (`matching: displacement`) on the trapezoidal pair.
# The Dirichlet side imposes Π_D u_N and derives its interface velocity and
# acceleration from its own Newmark relations. With the same Newmark parameters
# and step on both sides and coupled initial data, these relations are linear in
# (u, v, a), so matching the displacement at every step is equivalent to
# matching the velocity: the two runs agree to the Schwarz tolerance, and Ẽ is
# conserved.
@testset "Matching DN: Displacement Matching (II)" begin
    sim, header, data = run_constrained_dn(false; constraint="displacement")
    @test sim.failed == false
    @test length(data) == 51
    bc_D, bc_N = dn_bc_of(sim.subsims[1]), dn_bc_of(sim.subsims[2])
    @test bc_D.constraint == :displacement && bc_N.constraint == :displacement
    r = Norma.dn_interface_residuals(sim)[1][2]
    @test r.displacement_jump ≤ 1.0e-10
    @test r.velocity_jump ≤ 1.0e-10
    @test r.force_residual ≤ 1.0e-10
    @test invariant_drift(header, data) < 1.0e-10
    # Newmark relations of the Dirichlet side over the last step.
    integrator = sim.subsims[1].integrator
    Δt = integrator.time_step
    β, γ = integrator.β, integrator.γ
    u_n = sim.controller.stop_disp[1]
    v_n = sim.controller.stop_velo[1]
    a_n = sim.controller.stop_acce[1]
    for i_global in bc_D.global_from_local_map, comp in 1:3
        k = 3 * (i_global - 1) + comp
        a = integrator.acceleration[k]
        @test integrator.displacement[k] ≈ u_n[k] + Δt * v_n[k] + (0.5 - β) * Δt^2 * a_n[k] + β * Δt^2 * a rtol =
            1.0e-12 atol = 1.0e-15
        @test integrator.velocity[k] ≈ v_n[k] + (1.0 - γ) * Δt * a_n[k] + γ * Δt * a rtol = 1.0e-12 atol = 1.0e-12
    end
    fields_u = [vcat(vec(s.model.displacement), vec(s.model.velocity)) for s in sim.subsims]
    sim_v, _, _ = run_constrained_dn(false)
    @test sim_v.failed == false
    fields_v = [vcat(vec(s.model.displacement), vec(s.model.velocity)) for s in sim_v.subsims]
    for (a, b) in zip(fields_u, fields_v)
        @test norm(a - b) ≤ 1.0e-9 * norm(b)
    end
end

@testset "Matching DN: Coupled Initial Acceleration, First Step (II)" begin
    sim, header, data = run_constrained_dn(false; num_steps=2)
    @test sim.failed == false
    i_e1 = findfirst(==("total_energy"), header)
    i_invariant = findfirst(==("rate_invariant_total"), header)
    @test abs(data[2][i_invariant] / data[1][i_invariant] - 1.0) < 1.0e-9
    @test abs(data[2][i_e1] / data[1][i_e1] - 1.0) < 1.0e-9
    r = Norma.dn_interface_residuals(sim)[1][2]
    @test r.acceleration_jump ≤ 1.0e-10
end

@testset "Matching DN: Coupled Initial Acceleration, First Steps (EE)" begin
    sim, header, data = run_constrained_dn(true; num_steps=10)
    @test sim.failed == false
    i_invariant = findfirst(==("rate_invariant_total"), header)
    @test abs(data[2][i_invariant] / data[1][i_invariant] - 1.0) < 1.0e-9
    # Without the coupled initial acceleration the staggered energy E of the
    # explicit pair alternated by about 1e-3 from step to step.
    @test energy_drift(header, data, "total_energy") < 1.0e-9
    i_u = findfirst(h -> startswith(h, "displacement_jump"), header)
    @test maximum(row[i_u] for row in data) < 1.0e-12
end

# Different time steps: the free side takes four substeps per stop and the
# clamped side, with the coarse step, is the Dirichlet side. This is the r = 1
# multirate scheme of Connors, Owen, Kuberry, and Bochev (2024) and the scheme
# of Prakash and Hjelmstad (2004): the fine side receives the coarse side's
# reaction interpolated linearly between the two end values of the window, the
# coarse side imposes the fine side's velocity at the window end, and the
# interface terms of the rate invariant balance cancel (Connors et al. Eq.
# (103)), so Ẽ is constant to the Schwarz tolerance accumulated over the
# stops.
self_name(bc) = Norma.self_subsim_of(bc).name

function dn_bc_pair(sim)
    bcs = [dn_bc_of(s) for s in sim.subsims]
    bc_D = first(filter(b -> b.is_dirichlet, bcs))
    bc_N = first(filter(b -> !b.is_dirichlet, bcs))
    return bc_D, bc_N
end

for explicit in (false, true)
    @testset "Matching DN: Subcycled 4:1 ($(explicit ? "EE" : "II"))" begin
        sim, header, data = run_constrained_dn(explicit; num_steps=20, free_substeps=4, clamped_dirichlet=true)
        @test sim.failed == false
        @test length(data) == 21
        i_invariant = findfirst(==("rate_invariant_total"), header)
        invariant = [row[i_invariant] for row in data]
        @test maximum(abs.(invariant ./ invariant[1] .- 1.0)) < 1.0e-10
        bc_D, bc_N = dn_bc_pair(sim)
        @test self_name(bc_D) == "cantilever-clamped"
        # The force on the fine (Neumann) side at each substep of the last stop
        # is the linear interpolation between the window's end values.
        t_start = sim.controller.prev_time
        t_end = sim.controller.time
        entries = filter(e -> e[1] >= t_start - 1.0e-15, bc_N.transferred_force_history)
        @test length(entries) == 5
        f_start, f_end = entries[1][2], entries[end][2]
        for (t, f) in entries
            s = (t - t_start) / (t_end - t_start)
            @test norm(f - ((1.0 - s) * f_start + s * f_end)) ≤ 1.0e-10 * norm(f_end)
        end
        # The coarse (Dirichlet) side imposes the fine side's velocity at the
        # window end, and the force residual is below the tolerance at every
        # fine substep.
        r_end = Norma.dn_interface_residuals(sim)[1][2]
        @test r_end.velocity_jump ≤ 1.0e-10
        r_sub = Norma.dn_substep_residuals(sim, bc_D)
        @test r_sub.force_residual ≤ 1.0e-10
        # The interface impulse over the last window matches.
        net, relative = Norma.dn_impulse_residual(sim, bc_D)
        @test relative ≤ 1.0e-10
    end
    @testset "Matching DN: Subcycled 4:1 ($(explicit ? "EE" : "II")), t = 0" begin
        sim, _, _ = run_constrained_dn(explicit; free_substeps=4, clamped_dirichlet=true, initialize_only=true)
        r = Norma.dn_interface_residuals(sim)[1][2]
        @test r.force_residual ≤ 1.0e-11
        @test r.acceleration_jump ≤ 1.0e-11
    end
end

# Direct interface solve of the explicit pair (`interface solve: direct`): the
# interface force is solved from the velocity constraint once per stop. It must
# reproduce the iteration converged to 1e-12 and conserve Ẽ. Conforming meshes
# and the 2:1 nonconforming meshes of the cantilever-dn-nonconforming example.
const nonconforming_beam = "../examples/nonoverlap/dynamic-same-step/cantilever-dn-nonconforming"
for (label, mesh_dir) in (("conforming", constrained_dn_example), ("2:1 nonconforming", nonconforming_beam))
    @testset "Matching DN: Direct Interface Solve, EE $label" begin
        runs = Dict{String,Any}()
        for solve in ("iterative", "direct")
            sim, header, data = run_constrained_dn(true; interface_solve=solve, mesh_dir)
            @test sim.failed == false
            fields = [vcat(vec(s.model.displacement), vec(s.model.velocity)) for s in sim.subsims]
            runs[solve] = (header, data, fields)
        end
        header, it_data, it_fields = runs["iterative"]
        _, dir_data, dir_fields = runs["direct"]
        for (a, b) in zip(it_fields, dir_fields)
            @test norm(a - b) ≤ 1.0e-12 * norm(a)
        end
        for column in ("total_energy", "rate_invariant_total")
            i = findfirst(==(column), header)
            @test maximum(abs(r[i] / q[i] - 1.0) for (r, q) in zip(dir_data, it_data)) ≤ 1.0e-12
        end
        @test invariant_drift(header, dir_data) < 1.0e-12
    end
end

# Aitken relaxation of the constrained exchange forms its factor from the
# projected interface trace Π_D v_N that the Dirichlet side receives. Formed
# from the partner's whole velocity field, whose interior values follow a map
# of zero gain, the factor reached 0.9996 and 0.99999 in the third and fourth
# iterations of every stop of the implicit pair, where the interface gain is
# near -1, and those iterations did not contract (at a step of 5e-7 s, 8.2
# iterations per stop and stalls at the solver floor, against 5.8 and none from
# the trace; at 1e-6 s the recursive form stopped on the stall rule at the
# first stop). At 1e-6 s the trace gives 8 iterations per stop, against 7 for
# the fixed factor 0.5, which is near the optimum for this pair.
for relaxation in ("aitken secant", "aitken recursive")
    @testset "Matching DN: $relaxation on the interface trace (II)" begin
        sim, header, data = run_constrained_dn(false; num_steps=10, relaxation)
        @test sim.failed == false
        iterations = sim.controller.schwarz_iters[1:10]
        @test maximum(iterations) ≤ 8
        @test invariant_drift(header, data) < 1.0e-10
    end
end

# Anderson acceleration (`relaxation: anderson`) of the constrained datum.
@testset "Matching DN: Anderson acceleration (II)" begin
    sim, header, data = run_constrained_dn(false; num_steps=10, relaxation="anderson")
    @test sim.failed == false
    # The Aitken forms take at most 8 iterations per stop on this case.
    @test maximum(sim.controller.schwarz_iters[1:10]) ≤ 8
    @test invariant_drift(header, data) < 1.0e-10
end

# Interface gain of the conforming implicit pair at a stop, G = -Π H_N Πᵀ H_D⁻¹
# (H = C M̃⁻¹ Cᵀ with M̃ = M + β Δt² K), built by probing with unit interface
# data, and the equivalence of Anderson acceleration with mixing 1 and full
# depth to GMRES on the affine interface map x -> A x + b (Walker and Ni 2011,
# Theorem 2.2): the Anderson iterate x_{k+1} equals A x_k + b with x_k the k-th
# GMRES iterate for (I - A) x = b from the same start.
function interface_flexibility(subsim, bc, Δt, β)
    model = subsim.model
    Norma.evaluate(model, subsim.integrator, subsim.solver)
    Mt = model.mass + β * Δt^2 * model.stiffness
    fixed = Norma.prescribed_dofs(model)
    free = findall(.!fixed)
    position = Dict(d => i for (i, d) in enumerate(free))
    dofs = [3 * (n - 1) + c for n in bc.global_from_local_map for c in 1:3]
    F = cholesky(Symmetric(Matrix(Mt[free, free])))
    H = zeros(length(dofs), length(dofs))
    for (j, d) in enumerate(dofs)
        e = zeros(length(free))
        e[position[d]] = 1.0
        x = F \ e
        H[:, j] = [x[position[q]] for q in dofs]
    end
    return H
end

function gmres_iterates(A::Matrix{Float64}, b::Vector{Float64}, x0::Vector{Float64}, k_max::Int)
    n = length(b)
    L = Matrix(1.0I, n, n) - A
    r0 = b - L * x0
    V = zeros(n, k_max + 1)
    Hh = zeros(k_max + 1, k_max)
    V[:, 1] = r0 / norm(r0)
    iterates = Vector{Float64}[]
    for k in 1:k_max
        w = L * V[:, k]
        for i in 1:k
            Hh[i, k] = dot(V[:, i], w)
            w -= Hh[i, k] * V[:, i]
        end
        Hh[k + 1, k] = norm(w)
        V[:, k + 1] = w / Hh[k + 1, k]
        e1 = zeros(k + 1)
        e1[1] = norm(r0)
        y = Hh[1:(k + 1), 1:k] \ e1
        push!(iterates, x0 + V[:, 1:k] * y)
    end
    return iterates
end

@testset "Matching DN: Anderson acceleration is GMRES on the interface map" begin
    sim, _, _ = run_constrained_dn(false; initialize_only=true)
    free, clamped = sim.subsims[1], sim.subsims[2]
    bc_D, bc_N = dn_bc_of(free), dn_bc_of(clamped)
    Δt = free.integrator.time_step
    H_D = interface_flexibility(free, bc_D, Δt, free.integrator.β)
    H_N = interface_flexibility(clamped, bc_N, Δt, clamped.integrator.β)
    P = kron(bc_D.dirichlet_projector, Matrix(1.0I, 3, 3))
    A = -(P * H_N * transpose(P)) / H_D
    n = size(A, 1)
    b = [sin(0.37 * i) for i in 1:n]
    x0 = [cos(0.11 * i) for i in 1:n]
    G(x) = A * x + b
    k_max = 6
    history = Norma.AndersonHistory()
    x = x0
    anderson = Vector{Float64}[]
    for k in 0:k_max
        g = G(x)
        x = Norma.anderson_step!(history, x, g, g - x, 1.0, k_max + 1)
        push!(anderson, x)
    end
    gmres = gmres_iterates(A, b, x0, k_max)
    for k in 1:k_max
        @test norm(anderson[k + 1] - G(gmres[k])) ≤ 1.0e-8 * norm(G(gmres[k]))
    end
end

# The cross-mass matrix B of a matching pair is integrated by segment-based
# quadrature: each facet of the integration side is split into its
# intersections with the facets of the other side, so that B is exact on
# planar interfaces whatever the mesh ratio. On the 1:0.75 meshes (clamped
# part 8.467 mm, free part 6.35 mm) the facets are not nested; with the
# facet's own Gauss rule (before 2026-10-18) the matrix integrated on one side
# gave the projection onto the other side a partition-of-unity error of about
# 1e-2, and the Dirichlet projector did not reproduce linear fields. Checked:
# both projections reproduce constants and the interface coordinates to
# roundoff, and the force transfer carries the consistent nodal loads of a
# uniform traction exactly.
@testset "Matching DN: Cross Mass Matrix by Segmented Quadrature" begin
    for clamped_dirichlet in (false, true)
        sim, _, _ = run_constrained_dn(
            false;
            mesh_dir=nonconforming_beam,
            clamped_mesh="cantilever-clamped-ratio-0.75.g",
            clamped_dirichlet,
            initialize_only=true,
        )
        free, clamped = sim.subsims[1], sim.subsims[2]
        bc_D = clamped_dirichlet ? dn_bc_of(clamped) : dn_bc_of(free)
        bc_N = clamped_dirichlet ? dn_bc_of(free) : dn_bc_of(clamped)
        model_D = clamped_dirichlet ? clamped.model : free.model
        model_N = clamped_dirichlet ? free.model : clamped.model
        @test bc_D.is_dirichlet && !bc_N.is_dirichlet
        P_D = bc_D.dirichlet_projector
        P_N = bc_N.dirichlet_projector
        @test size(P_D, 1) != size(P_D, 2)
        @test maximum(abs.(P_D * ones(size(P_D, 2)) .- 1.0)) ≤ 1.0e-12
        @test maximum(abs.(P_N * ones(size(P_N, 2)) .- 1.0)) ≤ 1.0e-12
        X_D = model_D.reference[:, bc_D.global_from_local_map]
        X_N = model_N.reference[:, bc_N.global_from_local_map]
        scale = maximum(abs.(X_D))
        @test maximum(abs.(X_D' - P_D * X_N')) ≤ 1.0e-12 * scale
        @test maximum(abs.(X_N' - P_N * X_D')) ≤ 1.0e-12 * scale
        W_D = Norma.get_square_projection_matrix(model_D, bc_D)
        W_N = Norma.get_square_projection_matrix(model_N, bc_N)
        # Consistent nodal loads of a uniform traction on D transferred to N.
        loads_D = W_D * ones(size(W_D, 1))
        loads_N = W_N * ones(size(W_N, 1))
        @test maximum(abs.(bc_N.neumann_projector * loads_D - loads_N)) ≤ 1.0e-12 * maximum(abs.(loads_N))
    end
end

# Quadratic facets: the TRI6 facets of TETRA10 elements are integrated by the
# segmented rule when they are straight-sided, against the QUAD4 facets of a
# HEX8 side (the 1:0.75 clamped part against a TET10 free part) and against
# TRI6 facets (the TET10 cuboids of the AHeaD example). A curved TRI6 facet
# is refused and takes the point rule.
@testset "Matching DN: Segmented Quadrature on Quadratic Facets" begin
    function checks(sim, clamped_dirichlet)
        a, b = sim.subsims[1], sim.subsims[2]
        D, N = clamped_dirichlet ? (b, a) : (a, b)
        bc_D, bc_N = dn_bc_of(D), dn_bc_of(N)
        @test bc_D.is_dirichlet && !bc_N.is_dirichlet
        P_D, P_N = bc_D.dirichlet_projector, bc_N.dirichlet_projector
        @test maximum(abs.(P_D * ones(size(P_D, 2)) .- 1.0)) ≤ 1.0e-12
        @test maximum(abs.(P_N * ones(size(P_N, 2)) .- 1.0)) ≤ 1.0e-12
        X_D = D.model.reference[:, bc_D.global_from_local_map]
        X_N = N.model.reference[:, bc_N.global_from_local_map]
        scale = maximum(abs.(X_D))
        @test maximum(abs.(X_D' - P_D * X_N')) ≤ 1.0e-12 * scale
        @test maximum(abs.(X_N' - P_N * X_D')) ≤ 1.0e-12 * scale
        W_D = Norma.get_square_projection_matrix(D.model, bc_D)
        W_N = Norma.get_square_projection_matrix(N.model, bc_N)
        loads_D = W_D * ones(size(W_D, 1))
        loads_N = W_N * ones(size(W_N, 1))
        @test maximum(abs.(bc_N.neumann_projector * loads_D - loads_N)) ≤ 1.0e-12 * maximum(abs.(loads_N))
        return nothing
    end
    function captured(f)
        mktemp() do path, io
            sim = redirect_stdout(io) do
                f()
            end
            flush(io)
            seekstart(io)
            sim, replace(read(io, String), r"\s+" => " ")
        end
    end
    # HEX8 against TET10 on the beam, either side as the Dirichlet side. The
    # default five-point rule of TETRA10 gives an indefinite consistent mass,
    # which the direct factorization of the initial acceleration rejects; the
    # fourteen-point rule of degree four is positive definite.
    function tet10_quadrature(name, sub)
        name == "cantilever-free" && (sub["model"]["num integration points"] = Dict("cantilever_free" => 14))
        return nothing
    end
    for clamped_dirichlet in (false, true)
        sim, out = captured() do
            sim, _, _ = run_constrained_dn(
                false;
                mesh_dir=nonconforming_beam,
                clamped_mesh="cantilever-clamped-ratio-0.75.g",
                free_mesh="cantilever-free-tet10.g",
                clamped_dirichlet,
                initialize_only=true,
                edit_sub=tet10_quadrature,
            )
            sim
        end
        @test !occursin("Interface transfer", out)
        checks(sim, clamped_dirichlet)
    end
    # TET10 against TET10 on the cuboids of the AHeaD example.
    tet10_dir = "../examples/ahead/nonoverlap/cuboid"
    for name in ("cuboid-1", "cuboid-2")
        sub = YAML.load_file("$tet10_dir/dynamic-tet10/$name.yaml"; dicttype=Norma.Parameters)
        sub["input mesh file"] = "$tet10_dir/cuboid-tet10-$(name[end]).g"
        sub["boundary conditions"]["Schwarz DN nonoverlap"][1]["matching"] = "velocity"
        YAML.write_file("$name.yaml", sub)
    end
    params = YAML.load_file("$tet10_dir/dynamic-tet10/cuboids.yaml"; dicttype=Norma.Parameters)
    params["name"] = "cuboids-tet10"
    delete!(params, "relative tolerance")
    delete!(params, "absolute tolerance")
    params["matching relative tolerance"] = 1.0e-12
    sim, out = captured() do
        sim = Norma.create_simulation(params)
        Norma.sync_control_time(sim)
        Norma.initialize(sim)
        sim
    end
    @test !occursin("Interface transfer", out)
    # domains: cuboid-2 (Neumann) first, cuboid-1 (Dirichlet) second.
    checks(sim, true)
    for subsim in sim.subsims
        Exodus.close(subsim.params["input_mesh"])
        Exodus.close(subsim.params["output_mesh"])
    end
    foreach(f -> rm(f; force=true), ["cuboid-1.yaml", "cuboid-2.yaml", "cuboid-1.e", "cuboid-2.e"])
    # A TRI6 facet with a midside node off the edge midpoint is curved.
    straight = [0.0 1.0 0.0 0.5 0.5 0.0; 0.0 0.0 1.0 0.0 0.5 0.5; 0.0 0.0 0.0 0.0 0.0 0.0]
    @test Norma.is_segmentable_facet(straight)
    @test Norma.facet_midside_departure(straight) == 0.0
    curved = copy(straight)
    curved[3, 5] = 1.0e-6
    @test Norma.facet_midside_departure(curved) ≈ 1.0e-6 / sqrt(2.0)
    @test !Norma.is_segmentable_facet(curved)
end

# A subcycled Dirichlet side with `matching at: window end` imposes the matched
# velocity at the end of the window only and carries at its interior substeps
# the ramp of its own window reaction: the window balance of the rate invariant
# cancels and the impulse of the ramp equals the impulse the coarse side
# receives, so the fine-step side can be the Dirichlet side, which the substep
# constraint excludes (rate invariant conserved, impulse residual of the order
# of the impulse). The beam with the free part subcycled 4:1 as the Dirichlet
# side, 50 windows, Anderson acceleration on the datum and the reaction.
@testset "Matching DN: Window-End Constraint on a Subcycled Dirichlet Side" begin
    for explicit in (false, true)
        sim, header, data = run_constrained_dn(
            explicit;
            free_substeps=4,
            clamped_dirichlet=false,
            relaxation="anderson",
            edit_bc=bc -> (bc["matching at"] = "window end"),
            edit_params=p -> (p["anderson depth"]=10; p["relaxation parameter"]=0.5),
        )
        @test !sim.failed
        i_E = findfirst(==("total_energy"), header)
        i_R = findfirst(==("rate_invariant_total"), header)
        i_I = findfirst(startswith("impulse_residual_relative"), header)
        E0, R0 = data[1][i_E], data[1][i_R]
        @test maximum(abs(row[i_E] / E0 - 1.0) for row in data) ≤ 1.0e-2
        @test maximum(abs(row[i_R] / R0 - 1.0) for row in data) ≤ 1.0e-10
        @test maximum(row[i_I] for row in data[2:end]) ≤ 1.0e-12
        # The criterion tests the end value of the ramp that the interior
        # substeps applied against the reaction at the end of the window; in
        # the last, converged window the two agree to the tolerance.
        bc_D, _ = dn_bc_pair(sim)
        @test bc_D.window_end_constraint
        @test length(bc_D.window_ramp_end) == length(bc_D.window_end_reaction_raw)
        r = bc_D.window_end_reaction_raw
        @test norm(r - bc_D.window_ramp_end) ≤ 1.0e-10 * norm(r)
    end
end

# The Dirichlet role is the input's; at setup the code warns when it is
# assigned against the role rules: the fine-step side as the Dirichlet side,
# or, at equal steps and integrators, the side with fewer interface nodes.
@testset "Matching DN: Warning on the Dirichlet Role" begin
    function setup_output(; kwargs...)
        mktemp() do path, io
            redirect_stdout(io) do
                run_constrained_dn(false; initialize_only=true, kwargs...)
            end
            flush(io)
            seekstart(io)
            # The log wraps long lines; collapse the whitespace before matching.
            replace(read(io, String), r"\s+" => " ")
        end
    end
    # Coarser mesh (clamped part, 2:1) as the Dirichlet side at equal steps.
    out = setup_output(; mesh_dir=nonconforming_beam, clamped_dirichlet=true)
    @test occursin("fewer interface nodes", out)
    out = setup_output(; mesh_dir=nonconforming_beam, clamped_dirichlet=false)
    @test !occursin("[WARN]", out)
    # Fine-step side (free part, 4:1) as the Dirichlet side.
    out = setup_output(; free_substeps=4, clamped_dirichlet=false)
    @test occursin("finer time step", out)
    out = setup_output(; free_substeps=4, clamped_dirichlet=true)
    @test !occursin("finer time step", out)
end
