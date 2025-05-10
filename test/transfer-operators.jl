# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
@testset "Transfer Operators" begin
    cp("../examples/contact/transfer-operators/transfer.yaml", "transfer.yaml"; force=true)
    cp("../examples/contact/transfer-operators/src.yaml", "src.yaml"; force=true)
    cp("../examples/contact/transfer-operators/dst.yaml", "dst.yaml"; force=true)
    cp("../examples/contact/transfer-operators/src.g", "src.g"; force=true)
    cp("../examples/contact/transfer-operators/dst.g", "dst.g"; force=true)
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
    src_bc_index = 4
    dst_bc_index = 4
    src_bc = src_model.boundary_conditions[src_bc_index]
    dst_bc = dst_model.boundary_conditions[dst_bc_index]
    src_side_set_id = src_bc.side_set_id
    dst_side_set_id = dst_bc.side_set_id
    expression_force = "1.0 * t"
    t = 1.0
    src_T = create_force(expression_force, src_mesh, src_side_set_id, t)
    Norma.norma_log(0, :info, "Source side set:         $(length(src_T)) nodes")
    dst_T_real = create_force(expression_force, dst_mesh, dst_side_set_id, t)
    Norma.norma_log(0, :info, "Destination side set:    $(length(dst_T_real)) nodes")
    H = Norma.get_square_projection_matrix(src_model, src_bc)
    L = Norma.get_rectangular_projection_matrix(dst_model, dst_bc, src_model, src_bc)
    dst_T = L * (H \ I) * src_T
    rel_er_tr = norm(dst_T - dst_T_real) / norm(dst_T_real)
    Norma.norma_logf(0, :summary, "Relative error (traction):     %.4e", rel_er_tr)
    @test norm(dst_T - dst_T_real) / norm(dst_T_real) ≈ 0.0 atol = 1.0e-08
    expression_disp = "1.0 * t"
    src_U = create_displacement(expression_disp, src_mesh, src_side_set_id, t)
    Norma.norma_log(0, :info, "Source side set:         $(length(src_U)) nodes")
    dst_U_real = create_displacement(expression_disp, dst_mesh, dst_side_set_id, t)
    Norma.norma_log(0, :info, "Destination side set:    $(length(dst_U_real)) nodes")
    W = Norma.get_square_projection_matrix(dst_model, dst_bc)
    L = Norma.get_rectangular_projection_matrix(dst_model, dst_bc, src_model, src_bc)
    dst_U = (W \ I) * L * src_U
    rel_er_disp = norm(dst_U - dst_U_real) / norm(dst_U_real)
    Norma.norma_logf(0, :summary, "Relative error (displacement): %.4e", rel_er_disp)
    @test norm(dst_U - dst_U_real) / norm(dst_U_real) ≈ 0.0 atol = 1.0e-08
    Exodus.close(src_sim.params["input_mesh"])
    Exodus.close(src_sim.params["output_mesh"])
    Exodus.close(dst_sim.params["input_mesh"])
    Exodus.close(dst_sim.params["output_mesh"])
    rm("src.g")
    rm("dst.g")
    rm("src.e")
    rm("dst.e")
end
