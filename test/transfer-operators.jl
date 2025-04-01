# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "transfer-operators" begin
    cp("../examples/contact/transfer_operators/transfer.yaml", "transfer.yaml"; force=true)
    cp("../examples/contact/transfer_operators/src.yaml", "src.yaml"; force=true)
    cp("../examples/contact/transfer_operators/dst.yaml", "dst.yaml"; force=true)
    cp("../examples/contact/transfer_operators/src.g", "src.g"; force=true)
    cp("../examples/contact/transfer_operators/dst.g", "dst.g"; force=true)
    input_file = "transfer.yaml"
    sim = Norma.create_simulation(input_file)
    src_sim = sim.subsims[1]
    dst_sim = sim.subsims[2]
    src_model = src_sim.model
    dst_model = dst_sim.model
    src_mesh = src_model.mesh
    dst_mesh = dst_model.mesh
    rm("transfer.yaml")
    rm("src.yaml")
    rm("dst.yaml")
    src_side_set_id = 5
    dst_side_set_id = 6
    src_T = get_boundary_traction_force(src_mesh, src_side_set_id)
    println("Source side set:         $(length(src_T)) nodes")
    dst_T_real = get_boundary_traction_force(dst_mesh, dst_side_set_id)
    println("Destination side set:    $(length(dst_T_real)) nodes")
    H = Norma.get_square_projection_matrix(src_model, src_side_set_id)
    L = Norma.get_rectangular_projection_matrix(src_model, src_side_set_id, dst_model, dst_side_set_id)
    dst_T = L * inv(H) * src_T
    rel_er_tr = norm(dst_T - dst_T_real) / norm(dst_T_real)
    @printf("Relative error (traction):     %.4e\n", rel_er_tr)
    @test norm(dst_T - dst_T_real) / norm(dst_T_real) ≈ 0.0 atol = 1.0e-08
    M = Norma.get_square_projection_matrix(dst_model, dst_side_set_id)
    src_u = ones(length(src_T))
    dst_u = inv(M) * L * src_u
    dst_u_real = ones(length(dst_T_real))
    rel_er_disp = norm(dst_u - dst_u_real) / norm(dst_u_real)
    @printf("Relative error (displacement): %.4e\n", rel_er_disp)
    @test norm(dst_u - dst_u_real) / norm(dst_u_real) ≈ 0.0 atol = 1.0e-08
    Exodus.close(src_sim.params["input_mesh"])
    Exodus.close(src_sim.params["output_mesh"])
    Exodus.close(dst_sim.params["input_mesh"])
    Exodus.close(dst_sim.params["output_mesh"])
    rm("src.g")
    rm("dst.g")
    rm("src.e")
    rm("dst.e")
end
