# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Test: quasi-static uniaxial tension of a unit cube with J2 plasticity
# (finite deformation, multiplicative split, Hencky elasticity, linear hardening).
#
# Parameters: E = 200 GPa, ν = 0.25, σy = 1 GPa, H = 20 GPa
# Loading:    ε_z ≈ 0.01 applied over 10 quasi-static steps
#
# Physical expectations:
#   • avg_disp[3] ≈ half the applied BC displacement (linear distribution, 0 → 0.01)
#   • yielding occurs at ε_z^yield = σy/E = 0.005; thereafter stress hardens
#   • σzz must satisfy  σy < σzz < σy + H·ε_z  (between perfectly-plastic and elastic)
#   • lateral (σxx, σyy) and shear stresses are near zero (uniaxial-like BCs)
#
# Regression values obtained from the reference implementation:
#   avg_disp[3]   ≈ 4.939e-3
#   avg_stress[3] ≈ 1.052e9 Pa  (below small-strain uniaxial 1.091 GPa due to
#                                finite-deformation and von Mises constraint effects)

@testset "Single Static Solid J2 Plasticity Cube" begin
    cp("../examples/materials/j2/cube.yaml", "cube.yaml"; force=true)
    cp("../examples/materials/j2/cube.g", "cube.g"; force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml")
    rm("cube.g")
    rm("cube.e")

    avg_disp = average_components(integrator.displacement)
    avg_stress = average_components(model.stress)

    # Average z-displacement ≈ half the applied BC value (linear field from 0 to 0.01)
    @test avg_disp[3] ≈ 5.000e-3 rtol = 1.0e-3

    # σzz must be above yield stress (confirming plastic loading occurred)
    @test avg_stress[3] > 1.0e9

    # σzz regression check (Simo-Hughes J2 with correct state handling)
    @test avg_stress[3] ≈ 1.088e9 rtol = 1.0e-3

    # Lateral and shear stresses should be small relative to σzz
    @test avg_stress[1] ≈ 0.0 atol = 2.0e7
    @test avg_stress[2] ≈ 0.0 atol = 2.0e7
    @test avg_stress[4] ≈ 0.0 atol = 2.0e7
    @test avg_stress[5] ≈ 0.0 atol = 2.0e7
    @test avg_stress[6] ≈ 0.0 atol = 2.0e7
end
