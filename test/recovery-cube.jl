# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Exodus

function _read_recovered_stress(exo_path::String)
    exo = ExodusDatabase(exo_path, "r")
    last_step = Exodus.read_number_of_time_steps(exo)
    components = Dict{String,Vector{Float64}}()
    for name in ("sigma_xx_n", "sigma_yy_n", "sigma_zz_n", "sigma_yz_n", "sigma_xz_n", "sigma_xy_n")
        components[name] = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, name))
    end
    Exodus.close(exo)
    return components
end

@testset "Recovery — lumped, linear elastic, uniaxial cube" begin
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    open("cube.yaml", "w") do io
        write(io, """
type: single
input mesh file: cube.g
output mesh file: cube.e
CSV output interval: 0
model:
  type: solid mechanics
  stress recovery: lumped
  material:
    blocks:
      cube: elastic
    elastic:
      model: linear elastic
      elastic modulus: 1.0e+09
      Poisson's ratio: 0.25
      density: 1000.0
time integrator:
  type: quasi static
  initial time: 0.0
  final time: 1.0
  time step: 0.1
boundary conditions:
  Dirichlet:
    - node set: nsx-
      component: x
      function: "0.0"
    - node set: nsy-
      component: y
      function: "0.0"
    - node set: nsz-
      component: z
      function: "0.0"
    - node set: nsz+
      component: z
      function: "1.0 * t"
solver:
  type: Hessian minimizer
  step: full Newton
  minimum iterations: 1
  maximum iterations: 16
  relative tolerance: 1.0e-14
  absolute tolerance: 1.0e-08
  linear solver relative tolerance: 1.0e-14
  linear solver absolute tolerance: 1.0e-08
""")
    end
    Norma.run("cube.yaml")
    σ = _read_recovered_stress("cube.e")
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    # Uniaxial linear-elastic: σ_zz ≈ E·ε ≈ 1e9, all others ≈ 0.
    @test all(isapprox.(σ["sigma_zz_n"], 1.0e9; rtol=1.0e-6))
    for c in ("sigma_xx_n", "sigma_yy_n", "sigma_xy_n", "sigma_yz_n", "sigma_xz_n")
        @test maximum(abs.(σ[c])) < 1.0e3   # ~1e-6 relative to σ_zz ≈ 1e9
    end
end

@testset "Recovery — consistent, linear elastic, uniaxial cube" begin
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    open("cube.yaml", "w") do io
        write(io, """
type: single
input mesh file: cube.g
output mesh file: cube.e
CSV output interval: 0
model:
  type: solid mechanics
  stress recovery: consistent
  material:
    blocks:
      cube: elastic
    elastic:
      model: linear elastic
      elastic modulus: 1.0e+09
      Poisson's ratio: 0.25
      density: 1000.0
time integrator:
  type: quasi static
  initial time: 0.0
  final time: 1.0
  time step: 0.1
boundary conditions:
  Dirichlet:
    - node set: nsx-
      component: x
      function: "0.0"
    - node set: nsy-
      component: y
      function: "0.0"
    - node set: nsz-
      component: z
      function: "0.0"
    - node set: nsz+
      component: z
      function: "1.0 * t"
solver:
  type: Hessian minimizer
  step: full Newton
  minimum iterations: 1
  maximum iterations: 16
  relative tolerance: 1.0e-14
  absolute tolerance: 1.0e-08
  linear solver relative tolerance: 1.0e-14
  linear solver absolute tolerance: 1.0e-08
""")
    end
    Norma.run("cube.yaml")
    σ = _read_recovered_stress("cube.e")
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    # Consistent recovery is exact for fields lying in the basis: σ is
    # piecewise-constant per element, but the L2-projection onto the
    # nodal basis still reproduces a uniform σ_zz at every node.
    @test all(isapprox.(σ["sigma_zz_n"], 1.0e9; rtol=1.0e-9))
    for c in ("sigma_xx_n", "sigma_yy_n", "sigma_xy_n", "sigma_yz_n", "sigma_xz_n")
        @test maximum(abs.(σ[c])) < 1.0e1   # tighter — consistent is exact
    end
end

@testset "Recovery — internal variables (J2, eqps)" begin
    cp("../examples/materials/j2/cube.g", "cube.g"; force=true)
    open("cube.yaml", "w") do io
        write(io, """
type: single
input mesh file: cube.g
output mesh file: cube.e
CSV output interval: 0
model:
  type: solid mechanics
  stress recovery: lumped
  recover internal variables: true
  material:
    blocks:
      cube: plastic
    plastic:
      model: j2 plasticity
      elastic modulus: 200.0e9
      Poisson's ratio: 0.25
      density: 1000.0
      yield stress: 1.0e9
      hardening modulus: 20.0e9
time integrator:
  type: quasi static
  initial time: 0.0
  final time: 1.0
  time step: 0.1
boundary conditions:
  Dirichlet:
    - node set: nsx-
      component: x
      function: "0.0"
    - node set: nsy-
      component: y
      function: "0.0"
    - node set: nsz-
      component: z
      function: "0.0"
    - node set: nsz+
      component: z
      function: "0.01 * t"
solver:
  type: Hessian minimizer
  step: full Newton
  minimum iterations: 1
  maximum iterations: 32
  relative tolerance: 1.0e-10
  absolute tolerance: 1.0e-06
  linear solver relative tolerance: 1.0e-10
  linear solver absolute tolerance: 1.0e-06
""")
    end
    Norma.run("cube.yaml")
    exo = ExodusDatabase("cube.e", "r")
    last_step = Exodus.read_number_of_time_steps(exo)
    nodal_names = Exodus.read_names(exo, NodalVariable)
    @test "eqps_n" in nodal_names
    eqps_n = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, "eqps_n"))
    Exodus.close(exo)
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    # Uniaxial yields uniform plastic strain across the cube; ε_z = 0.01,
    # ε_y = σy/E = 0.005, so plastic part ≈ 0.005 → eqps roughly that order.
    @test all(eqps_n .> 0.0)
    @test maximum(eqps_n) - minimum(eqps_n) < 1.0e-9   # uniform across the cube
end
