# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Same-geometry mesh-to-mesh field transfer (issue #161).  Verifies that
# transfer_field is exact (~machine precision) when the source field lies in
# the destination basis: cross-mesh transfer of f(X) = (X, 2Y, 3Z) on
# discretizations of the same physical cube [-0.5, 0.5]^3.

using Exodus

function _block_name_of(path::String)::String
    exo = ExodusDatabase(path, "r")
    name = Exodus.read_names(exo, Block)[1]
    Exodus.close(exo)
    return name
end

function _make_transfer_model(mesh_path::String, scratch_name::String, npts_override)
    block = _block_name_of(mesh_path)
    npts_yaml = npts_override === nothing ? "" : """
  num integration points:
    $block: $npts_override
"""
    yaml_path = scratch_name * ".yaml"
    cp(mesh_path, scratch_name * ".g"; force=true)
    open(yaml_path, "w") do io
        write(io, """
type: single
input mesh file: $(scratch_name).g
output mesh file: $(scratch_name).e
CSV output interval: 0
model:
  type: solid mechanics
  stress recovery: consistent
$npts_yaml  material:
    blocks:
      $block: elastic
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
    return Norma.create_simulation(yaml_path).model
end

function _apply_linear!(model)
    for i in 1:size(model.reference, 2)
        X, Y, Z = model.reference[:, i]
        model.displacement[:, i] = [X, 2Y, 3Z]
    end
end

function _ref_linear(model)
    out = similar(model.displacement)
    for i in 1:size(model.reference, 2)
        X, Y, Z = model.reference[:, i]
        out[:, i] = [X, 2Y, 3Z]
    end
    return out
end

@testset "Volumetric Field Transfer — linear in any basis" begin
    sources = (
        (:hex8_27, "../examples/single/static-solid/cube/standard/cube.g",                       nothing),
        (:hex8_8,  "../examples/single/implicit-dynamic-solid/pressure-bc-expand/cube.g",        nothing),
        (:tet4,    "../examples/element-types/tet4/cube/cube.g",                                 nothing),
        (:tet10,   "../examples/element-types/tet10/cube/cube.g",                                14),
    )
    models = Dict{Symbol,Norma.SolidMechanics}()
    scratch_files = String[]
    for (k, path, npts) in sources
        scratch = "tx_" * String(k)
        models[k] = _make_transfer_model(path, scratch, npts)
        push!(scratch_files, scratch * ".yaml", scratch * ".g", scratch * ".e")
    end
    keys_to_test = [first(t) for t in sources]
    # f(X) = (X, 2Y, 3Z) is in every basis we exercise (hex8 trilinear, tet4
    # linear, tet10 quadratic); cross-mesh L2 projection should reproduce it
    # to roundoff via the consistent mass.
    for src_k in keys_to_test, dst_k in keys_to_test
        src, dst = models[src_k], models[dst_k]
        _apply_linear!(src)
        out = Norma.transfer_field(src, src.displacement, dst)
        ref = _ref_linear(dst)
        err = maximum(abs.(out .- ref))
        @test err < 1.0e-12
    end
    for f in scratch_files
        rm(f; force=true)
    end
end
