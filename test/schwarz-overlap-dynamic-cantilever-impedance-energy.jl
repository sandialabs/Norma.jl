# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Exodus

# Energy-conservation regression test for the impedance-overlap Schwarz
# transmission condition
#
#   t - t_p + Z (u̇ - u̇_p) + α (u - u_p) = 0,
#
# with the partner traction t_p interpolated from consistently recovered
# nodal stress and the P/S-split tensor impedance Z built from the partner's
# material. The benchmark is a bent 3D cantilever released from a parabolic
# initial displacement, split into two overlapping subdomains with CONFORMING
# (node-aligned) meshes, so the interface transfer is exact and the energy
# behavior isolates the transmission condition itself.
#
# The energy is measured in a partition-of-unity (PU) metric that counts the
# overlap region once: elements whose centroid lies in the region covered by
# both meshes carry weight 1/2. Summing the two subdomains' energies instead
# double-counts the overlap and swings by ~15% as the pulse transits it, even
# for provably conservative coupling, which masks real defects.
#
# Discrimination at t = 3.0e-4 (300 steps), PU energy relative to the initial
# strain energy of the release:
#   - current formulation (consistent-flux t_p):  E/E0 ≈ 1.00 (0.08% at 1 ms)
#   - recovered-stress t_p (consistent recovery): E/E0 ≈ 0.99
#   - lumped instead of consistent recovery:      E/E0 ≈ 0.95
#   - missing partner-traction term (pre-fix):    E/E0 ≈ 0.67
# The 3% assertion below fails both historical defects. A converged classical
# DBC overlap run on these meshes conserves the PU energy to 0.1% (exact
# reference); with the consistent-flux partner traction the impedance
# condition matches that reference.

# PU-weighted total mechanical energy of the two-subdomain state. Strain
# energy comes from the per-element stored energy; kinetic energy is
# recomputed per element with hex8 shape functions and 2x2x2 Gauss
# quadrature (a consistent-mass kinetic energy — nodal lumping over-counts
# by tens of percent when the velocity field carries mesh-scale content).
# The PU weight is 1/2 when the element centroid lies inside the
# intersection of the two meshes' x-extents.
const _hex8_signs = [
    -1.0 -1.0 -1.0; 1.0 -1.0 -1.0; 1.0 1.0 -1.0; -1.0 1.0 -1.0;
    -1.0 -1.0 1.0; 1.0 -1.0 1.0; 1.0 1.0 1.0; -1.0 1.0 1.0
]
_hex8_N(ξ) = [
    0.125 * (1 + _hex8_signs[n, 1] * ξ[1]) * (1 + _hex8_signs[n, 2] * ξ[2]) * (1 + _hex8_signs[n, 3] * ξ[3])
    for n in 1:8
]
function _hex8_dN(ξ)
    dN = zeros(8, 3)
    for n in 1:8
        a, b, c = _hex8_signs[n, 1], _hex8_signs[n, 2], _hex8_signs[n, 3]
        dN[n, 1] = 0.125 * a * (1 + b * ξ[2]) * (1 + c * ξ[3])
        dN[n, 2] = 0.125 * b * (1 + a * ξ[1]) * (1 + c * ξ[3])
        dN[n, 3] = 0.125 * c * (1 + a * ξ[1]) * (1 + b * ξ[2])
    end
    return dN
end
const _gauss_pts = [[i / sqrt(3.0), j / sqrt(3.0), k / sqrt(3.0)] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)]

function pu_energy(models::Vector, mesh_files::Vector{String}, ρ::Float64)
    lo = maximum(minimum(m.reference[1, :]) for m in models)
    hi = minimum(maximum(m.reference[1, :]) for m in models)
    total = 0.0
    for (model, mesh_file) in zip(models, mesh_files)
        mesh = ExodusDatabase(mesh_file, "r")
        blocks = Exodus.read_sets(mesh, Block)
        try
            for (b, block) in enumerate(blocks)
                conn = Norma.get_block_connectivity(mesh, block.id)
                num_elems, num_nodes_per_elem = size(conn)
                for e in 1:num_elems
                    # get_block_connectivity's matrix is meant for LINEAR
                    # indexing (element-major raw order), matching its use in
                    # model.jl; conn[e, n] would scramble the node lists.
                    nodes = [conn[(e - 1) * num_nodes_per_elem + n] for n in 1:num_nodes_per_elem]
                    X = model.reference[:, nodes]
                    V = model.velocity[:, nodes]
                    centroid_x = sum(X[1, :]) / num_nodes_per_elem
                    w = (centroid_x > lo + 1.0e-12 && centroid_x < hi - 1.0e-12) ? 0.5 : 1.0
                    kinetic = 0.0
                    for ξ in _gauss_pts
                        detJ = det(X * _hex8_dN(ξ))
                        vg = V * _hex8_N(ξ)
                        kinetic += 0.5 * ρ * (vg[1]^2 + vg[2]^2 + vg[3]^2) * detJ
                    end
                    total += w * (model.stored_energy[b][e] + kinetic)
                end
            end
        finally
            Exodus.close(mesh)
        end
    end
    return total
end

@testset "Schwarz Overlap Dynamic Cantilever Impedance Energy" begin
    example = "../examples/overlap/dynamic-same-step/cantilever-conforming"
    files = [
        "cantilever-impedance.yaml",
        "cantilever-clamped-impedance.yaml",
        "cantilever-free-impedance.yaml",
        "cantilever-clamped.g",
        "cantilever-free.g",
    ]
    for f in files
        cp(joinpath(example, f), f; force=true)
    end
    sim = Norma.run("cantilever-impedance.yaml")
    model_free = sim.subsims[1].model
    model_clamped = sim.subsims[2].model

    @test sim.failed == false

    # The BC forces consistent nodal stress recovery: the partner-traction
    # term is built from recovered stress, and lumped recovery's O(h) error
    # was measured to dominate the interface energy error.
    for model in (model_free, model_clamped)
        @test model.recovery_data isa Norma.ConsistentRecovery
        @test size(model.recovered_stress, 1) == 6
    end

    # P/S-split tensor impedance from the (here identical) partner material:
    # Z_p = √(ρ(λ + 2μ)), Z_s = √(ρμ).
    E_mod = 6.895e9
    ν = 0.25
    ρ = 2768.0
    λ_lame = E_mod * ν / ((1 + ν) * (1 - 2ν))
    μ_lame = E_mod / (2 * (1 + ν))
    for model in (model_free, model_clamped)
        bcs = filter(
            bc -> bc isa Norma.SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition,
            model.boundary_conditions,
        )
        @test length(bcs) == 1
        @test bcs[1].impedance ≈ sqrt(ρ * (λ_lame + 2μ_lame)) rtol = 1.0e-12
        @test bcs[1].impedance_shear ≈ sqrt(ρ * μ_lame) rtol = 1.0e-12
        # Node-aligned interface: the consistent-flux partner traction must
        # be active (the O(h) recovered-stress dissipation is not).
        @test bcs[1].flux_patch isa Norma.FluxTractionPatch
    end

    # PU energy at the final time vs. the initial strain energy of the
    # parabolic release. E0 is the (deterministic) FE strain energy of the
    # initial condition assembled on the committed meshes, with the overlap
    # counted once; it changes only if the meshes, material, or initial
    # condition change.
    E0 = 1506.28
    E_final = pu_energy([model_clamped, model_free], ["cantilever-clamped.g", "cantilever-free.g"], ρ)
    @test abs(E_final / E0 - 1.0) < 0.03

    for f in files
        rm(f; force=true)
    end
    rm("cantilever-clamped.e"; force=true)
    rm("cantilever-free.e"; force=true)
end
