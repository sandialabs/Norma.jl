# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using NPZ



function LinearOpInfRom(params::Parameters)
    params["mesh smoothing"] = false
    fom_model = SolidMechanics(params)
    reference = fom_model.reference
    opinf_model_file = params["model"]["model-file"]
    opinf_model = NPZ.npzread(opinf_model_file)
    basis = opinf_model["basis"]
    _, _, reduced_dim = size(basis)
    num_dofs = reduced_dim
    time = 0.0
    failed = false
    null_vec = zeros(num_dofs)

    reduced_state = zeros(num_dofs)
    reduced_velocity = zeros(num_dofs)
    reduced_boundary_forcing = zeros(num_dofs)
    free_dofs = trues(num_dofs)
    boundary_conditions = Vector{BoundaryCondition}()
    return LinearOpInfRom(
        opinf_model,
        basis,
        reduced_state,
        reduced_velocity,
        reduced_boundary_forcing,
        null_vec,
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
        false,
    )
end

function QuadraticOpInfRom(params::Parameters)
    params["mesh smoothing"] = false
    fom_model = SolidMechanics(params)
    reference = fom_model.reference
    opinf_model_file = params["model"]["model-file"]
    opinf_model = NPZ.npzread(opinf_model_file)
    basis = opinf_model["basis"]
    _, _, reduced_dim = size(basis)
    num_dofs = reduced_dim
    time = 0.0
    failed = false
    null_vec = zeros(num_dofs)

    reduced_state = zeros(num_dofs)
    reduced_velocity = zeros(num_dofs)
    reduced_boundary_forcing = zeros(num_dofs)
    free_dofs = trues(num_dofs)
    boundary_conditions = Vector{BoundaryCondition}()
    return QuadraticOpInfRom(
        opinf_model,
        basis,
        reduced_state,
        reduced_velocity,
        reduced_boundary_forcing,
        null_vec,
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
        false,
    )
end

function CubicOpInfRom(params::Parameters)
    params["mesh smoothing"] = false
    fom_model = SolidMechanics(params)
    reference = fom_model.reference
    opinf_model_file = params["model"]["model-file"]
    opinf_model = NPZ.npzread(opinf_model_file)
    basis = opinf_model["basis"]
    _, _, reduced_dim = size(basis)
    num_dofs = reduced_dim
    time = 0.0
    failed = false
    null_vec = zeros(num_dofs)

    reduced_state = zeros(num_dofs)
    reduced_velocity = zeros(num_dofs)
    reduced_boundary_forcing = zeros(num_dofs)
    free_dofs = trues(num_dofs)
    boundary_conditions = Vector{BoundaryCondition}()
    return CubicOpInfRom(
        opinf_model,
        basis,
        reduced_state,
        reduced_velocity,
        reduced_boundary_forcing,
        null_vec,
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
        false,
    )
end

