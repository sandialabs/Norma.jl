# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Constrained Dirichlet-Neumann exchange (`constrained: true`) on the
# conforming cantilever of examples/nonoverlap/dynamic-same-step/cantilever-dn.
# The free subdomain is the Dirichlet side. Checked: the force transfer of the
# Neumann side is the transpose of the Dirichlet projector; at convergence the
# force applied on the Neumann side is the transposed projection of the
# Dirichlet side's d'Alembert reaction; the Dirichlet side's interface
# acceleration and displacement satisfy its own Newmark relations over the
# last step; and E2 of the differentiated system is conserved over 50 steps.
# The initial acceleration of a constrained pair is found by a Schwarz
# iteration at t = 0 (coupled_initial_acceleration!), so conservation holds
# from step 0; the last two testsets check the first step and the absence of
# the alternating interface acceleration difference that the per-subdomain
# initial acceleration left in the explicit pair.

using YAML
using LinearAlgebra

const constrained_dn_example = "../examples/nonoverlap/dynamic-same-step/cantilever-dn"

function constrained_dn_subdomain(name::String, explicit::Bool, constraint::String, dt::Float64)
    sub = YAML.load_file("$constrained_dn_example/$name.yaml"; dicttype=Norma.Parameters)
    if explicit
        sub["time integrator"] = Norma.Parameters(
            "type" => "central difference", "time step" => dt, "CFL" => 1.0, "γ" => 0.5
        )
        sub["solver"] = Norma.Parameters("type" => "explicit solver", "step" => "explicit")
    else
        sub["time integrator"]["time step"] = dt
        sub["solver"]["linear solver"] = "direct"
        sub["solver"]["linear solver relative tolerance"] = 1.0e-14
        sub["solver"]["maximum iterations"] = 4
        sub["solver"]["absolute tolerance"] = 1.0e-6
    end
    bc = sub["boundary conditions"]["Schwarz DN nonoverlap"][1]
    bc["constrained"] = true
    bc["constraint"] = constraint
    YAML.write_file("$name.yaml", sub)
    return nothing
end

function run_constrained_dn(explicit::Bool; constraint="velocity", num_steps=50, dt=1.0e-6, theta=0.5)
    for f in ["cantilever-clamped.g", "cantilever-free.g"]
        cp("$constrained_dn_example/$f", f; force=true)
    end
    constrained_dn_subdomain("cantilever-free", explicit, constraint, dt)
    constrained_dn_subdomain("cantilever-clamped", explicit, constraint, dt)
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
        "relative tolerance" => 1.0e-12,
        "absolute tolerance" => 1.0e-15,
        "relaxation parameter" => theta,
        "blended energy output" => true,
    )
    sim = Norma.run(params)
    rows = readlines("constrained-dn-energy.csv")
    header = split(rows[1], ",")
    data = [parse.(Float64, split(row, ",")) for row in rows[2:end]]
    for f in ["cantilever-clamped.g", "cantilever-free.g", "cantilever-clamped.e", "cantilever-free.e",
              "cantilever-clamped.yaml", "cantilever-free.yaml", "constrained-dn-energy.csv"]
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

e2_drift(header, data) = energy_drift(header, data, "e2_total")

@testset "Constrained DN: Operators, Reaction, Newmark Relations (II)" begin
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
    @test bc_D.square_projector * bc_D.dirichlet_projector ≈
        transpose(bc_N.square_projector * bc_N.dirichlet_projector) rtol = 1.0e-12
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

@testset "Constrained DN: E2 Conservation (II)" begin
    sim, header, data = run_constrained_dn(false)
    @test sim.failed == false
    @test length(data) == 51
    @test e2_drift(header, data) < 1.0e-10
end

@testset "Constrained DN: E2 Conservation (EE)" begin
    sim, header, data = run_constrained_dn(true)
    @test sim.failed == false
    @test length(data) == 51
    @test e2_drift(header, data) < 1.0e-10
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

@testset "Constrained DN: Displacement Constraint Requires Newmark" begin
    @test_throws Exception run_constrained_dn(true; constraint="displacement", num_steps=1)
    for f in ["cantilever-clamped.g", "cantilever-free.g", "cantilever-clamped.e", "cantilever-free.e",
              "cantilever-clamped.yaml", "cantilever-free.yaml", "constrained-dn-energy.csv"]
        rm(f; force=true)
    end
end

@testset "Constrained DN: Coupled Initial Acceleration, First Step (II)" begin
    sim, header, data = run_constrained_dn(false; num_steps=2)
    @test sim.failed == false
    i_e1 = findfirst(==("total_energy"), header)
    i_e2 = findfirst(==("e2_total"), header)
    @test abs(data[2][i_e2] / data[1][i_e2] - 1.0) < 1.0e-9
    @test abs(data[2][i_e1] / data[1][i_e1] - 1.0) < 1.0e-9
    r = Norma.dn_interface_residuals(sim)[1][2]
    @test r.acceleration_jump ≤ 1.0e-10
end

@testset "Constrained DN: Coupled Initial Acceleration, First Steps (EE)" begin
    sim, header, data = run_constrained_dn(true; num_steps=10)
    @test sim.failed == false
    i_e2 = findfirst(==("e2_total"), header)
    @test abs(data[2][i_e2] / data[1][i_e2] - 1.0) < 1.0e-9
    # Without the coupled initial acceleration the staggered energy E1 of the
    # explicit pair alternated by about 1e-3 from step to step.
    @test energy_drift(header, data, "total_energy") < 1.0e-9
    i_u = findfirst(h -> startswith(h, "displacement_jump"), header)
    @test maximum(row[i_u] for row in data) < 1.0e-12
end
