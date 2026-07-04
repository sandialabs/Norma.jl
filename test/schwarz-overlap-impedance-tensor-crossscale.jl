# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Cross-scaling test for the P/S-split tensor impedance of the
# impedance-overlap Schwarz BC on a BIMATERIAL pair: each subdomain's
# impedances must be the NEIGHBOR's characteristic values
# Z_p = √(ρ(λ + 2μ)), Z_s = √(ρμ) (optimized-Schwarz cross-scaling: the
# optimal transmission operator approximates the neighbor's
# Dirichlet-to-Neumann map). With identical materials this is
# indistinguishable from using one's own impedances, so the two subdomains
# here get materials whose impedances differ by a factor of 4
# (E and ρ both ×4 ⇒ same wave speeds, 4× impedance).

@testset "Schwarz Overlap Impedance Tensor Cross-Scaling" begin
    example = "../examples/overlap/dynamic-same-step/cantilever-conforming"
    for f in ("cantilever-clamped.g", "cantilever-free.g")
        cp(joinpath(example, f), f; force=true)
    end

    E_soft = 6.895e9
    ρ_soft = 2768.0
    E_stiff = 4 * E_soft
    ρ_stiff = 4 * ρ_soft
    ν = 0.25

    function subdomain_yaml(name, other, block, other_block, ss, other_ss, E_mod, ρ, clamped)
        dirichlet = clamped ? """
        boundary conditions:
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
          Schwarz impedance overlap:
        """ : """
        boundary conditions:
          Schwarz impedance overlap:
        """
        return """
        type: single
        input mesh file: cantilever-$name.g
        output mesh file: cantilever-$name.e
        model:
          type: solid mechanics
          material:
            blocks:
              cantilever_$block: hyperelastic
            hyperelastic:
              model: linear elastic
              elastic modulus: $E_mod
              Poisson's ratio: $ν
              density: $ρ
        time integrator:
          type: Newmark
          time step: 1.0e-06
          β: 0.25
          γ: 0.5
        $(dirichlet)    - side set: $ss
              source: cantilever-$other-impedance
              source block: cantilever_$other_block
              source side set: $other_ss
              robin parameter: 0.0
        solver:
          type: Hessian minimizer
          step: full Newton
          minimum iterations: 1
          maximum iterations: 16
          relative tolerance: 1.0e-10
          absolute tolerance: 2.54e-08
        """
    end

    write(
        "cantilever-clamped-impedance.yaml",
        subdomain_yaml("clamped", "free", "clamped", "free", "ssx+", "ssx-", E_soft, ρ_soft, true),
    )
    write(
        "cantilever-free-impedance.yaml",
        subdomain_yaml("free", "clamped", "free", "clamped", "ssx-", "ssx+", E_stiff, ρ_stiff, false),
    )
    write(
        "cantilever-impedance.yaml",
        """
        type: multi
        domains: ["cantilever-free-impedance.yaml", "cantilever-clamped-impedance.yaml"]
        Exodus output interval: 0
        CSV output interval: 0
        initial time: 0.0
        final time: 2.0e-6
        time step: 1.0e-6
        minimum iterations: 1
        maximum iterations: 128
        relative tolerance: 1.0e-08
        absolute tolerance: 2.54e-08
        """,
    )

    sim = Norma.run("cantilever-impedance.yaml")
    model_free = sim.subsims[1].model      # stiff material
    model_clamped = sim.subsims[2].model   # soft material

    for f in (
        "cantilever-impedance.yaml",
        "cantilever-clamped-impedance.yaml",
        "cantilever-free-impedance.yaml",
        "cantilever-clamped.g",
        "cantilever-free.g",
        "cantilever-clamped.e",
        "cantilever-free.e",
    )
        rm(f; force=true)
    end

    @test sim.failed == false

    function ps_impedances(E_mod, ρ)
        λ_lame = E_mod * ν / ((1 + ν) * (1 - 2ν))
        μ_lame = E_mod / (2 * (1 + ν))
        return sqrt(ρ * (λ_lame + 2μ_lame)), sqrt(ρ * μ_lame)
    end
    Zp_soft, Zs_soft = ps_impedances(E_soft, ρ_soft)
    Zp_stiff, Zs_stiff = ps_impedances(E_stiff, ρ_stiff)
    @test Zp_stiff ≈ 4 * Zp_soft rtol = 1.0e-12  # sanity: contrast is 4x

    imp_bc(model) = first(
        bc for bc in model.boundary_conditions if
        bc isa Norma.SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition
    )

    # The soft (clamped) subdomain must carry the STIFF neighbor's impedances
    # and vice versa.
    @test imp_bc(model_clamped).impedance ≈ Zp_stiff rtol = 1.0e-12
    @test imp_bc(model_clamped).impedance_shear ≈ Zs_stiff rtol = 1.0e-12
    @test imp_bc(model_free).impedance ≈ Zp_soft rtol = 1.0e-12
    @test imp_bc(model_free).impedance_shear ≈ Zs_soft rtol = 1.0e-12

    # Consistent stress recovery is force-enabled on both sides.
    @test model_clamped.recovery_data isa Norma.ConsistentRecovery
    @test model_free.recovery_data isa Norma.ConsistentRecovery

    # Node-aligned interface: the consistent-traction partner traction is active.
    @test imp_bc(model_clamped).traction_patch isa Norma.ConsistentTractionPatch
    @test imp_bc(model_free).traction_patch isa Norma.ConsistentTractionPatch
end
