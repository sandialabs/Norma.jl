# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Exodus

# Read the six Voigt stress components from the last time step of an Exodus
# file.  `qualifier` selects which projected field to read: "cons" or "lump"
# for the 'both' recovery mode (names like sigma_xx_cons_n), or the empty
# string (default) for single-type recovery (names like sigma_xx_n).  The
# returned dictionary is always keyed by the bare "sigma_*_n" component name.
function _read_recovered_stress(exo_path::String, qualifier::String = "")
    exo = ExodusDatabase(exo_path, "r")
    last_step = Exodus.read_number_of_time_steps(exo)
    components = Dict{String,Vector{Float64}}()
    for axes in ("xx", "yy", "zz", "yz", "xz", "xy")
        key = "sigma_" * axes * "_n"
        read_name = isempty(qualifier) ? key : "sigma_" * axes * "_" * qualifier * "_n"
        components[key] = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, read_name))
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
  nodal recovery:
    method: lumped
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
  nodal recovery:
    method: consistent
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

@testset "Recovery — both (lumped + consistent), linear elastic, uniaxial cube" begin
    cp("../examples/single/static-solid/cube/standard/cube.g", "cube.g"; force=true)
    open("cube.yaml", "w") do io
        write(io, """
type: single
input mesh file: cube.g
output mesh file: cube.e
CSV output interval: 0
model:
  type: solid mechanics
  nodal recovery:
    method: both
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
    σ_l = _read_recovered_stress("cube.e", "lump")
    σ_c = _read_recovered_stress("cube.e", "cons")
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    # Both projections should recover the uniform uniaxial stress field.
    @test all(isapprox.(σ_l["sigma_zz_n"], 1.0e9; rtol=1.0e-6))
    for c in ("sigma_xx_n", "sigma_yy_n", "sigma_xy_n", "sigma_yz_n", "sigma_xz_n")
        @test maximum(abs.(σ_l[c])) < 1.0e3
    end
    @test all(isapprox.(σ_c["sigma_zz_n"], 1.0e9; rtol=1.0e-9))
    for c in ("sigma_xx_n", "sigma_yy_n", "sigma_xy_n", "sigma_yz_n", "sigma_xz_n")
        @test maximum(abs.(σ_c[c])) < 1.0e1
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
  nodal recovery:
    method: lumped
    internal variables: true
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

@testset "Recovery — von Mises stress and deformation gradient (J2 cube)" begin
    cp("../examples/materials/j2/cube.g", "cube.g"; force=true)
    open("cube.yaml", "w") do io
        write(io, """
type: single
input mesh file: cube.g
output mesh file: cube.e
CSV output interval: 0
model:
  type: solid mechanics
  nodal recovery:
    method: lumped
    von mises stress: true
    deformation gradient: true
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
    @test "von_mises_stress_n" in nodal_names
    @test "F_xx_n" in nodal_names
    @test "F_yy_n" in nodal_names
    @test "F_zz_n" in nodal_names
    @test "F_xy_n" in nodal_names
    vm = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, "von_mises_stress_n"))
    Fxx = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, "F_xx_n"))
    Fyy = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, "F_yy_n"))
    Fzz = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, "F_zz_n"))
    Fxy = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, "F_xy_n"))
    Fyz = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, "F_yz_n"))
    Fxz = Vector{Float64}(Exodus.read_values(exo, NodalVariable, last_step, "F_xz_n"))
    Exodus.close(exo)
    rm("cube.yaml"; force=true)
    rm("cube.g"; force=true)
    rm("cube.e"; force=true)
    # 1D J2, linear isotropic hardening, ε_z = 0.01:
    #   σ_zz = (σy + H·ε_z) / (1 + H/E)
    #   with σy = 1e9, H = 20e9, E = 200e9  →  σ_zz ≈ 1.0909e9
    # Uniaxial keeps σ_xx = σ_yy = 0, so σ_vm = σ_zz at every QP, and the
    # projected nodal vM should reproduce that uniform value.
    σ_vm_expected = (1.0e9 + 20.0e9 * 0.01) / (1.0 + 20.0e9 / 200.0e9)
    @test all(isapprox.(vm, σ_vm_expected; rtol=1.0e-2))
    @test maximum(vm) - minimum(vm) < 1.0e5   # uniform across the cube
    # F_zz = 1 + ε_z = 1.01 (small-strain regime); off-diagonal F ≈ 0
    @test all(isapprox.(Fzz, 1.01; atol=1.0e-3))
    @test all(Fxx .< 1.0) && all(Fyy .< 1.0)   # lateral Poisson contraction
    @test maximum(abs.(Fxy)) < 1.0e-3
    @test maximum(abs.(Fyz)) < 1.0e-3
    @test maximum(abs.(Fxz)) < 1.0e-3
end
