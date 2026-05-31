# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using NPZ
using PyCall

function torch_tensor_vector(tensor)
    values = tensor.detach().cpu().tolist()
    if ndims(values) == 1
        return Vector{Float64}(values)
    end
    return Vector{Float64}(values[1, :])
end

function torch_tensor_matrix(tensor)
    values = tensor.detach().cpu().tolist()
    if ndims(values) == 2
        return Matrix{Float64}(values)
    end
    return Matrix{Float64}(values[1, :, :])
end

function flush_pycall_finalizers()
    GC.gc()
    return nothing
end

function require_single_threaded_nn_opinf()
    num_threads = Base.Threads.nthreads()
    if num_threads != 1
        norma_abortf(
            "Neural network OpInf ROM uses PyCall/Torch and currently requires exactly one Julia thread; current thread count is %d. Run with `julia --threads 1 ...`.",
            num_threads,
        )
    end
    return nothing
end

function NeuralNetworkOpInfRom(params::Dict{String,Any})
    require_single_threaded_nn_opinf()
    params["mesh smoothing"] = false
    fom_model = SolidMechanics(params)
    reference = fom_model.reference
    opinf_model_directory = params["model"]["model-directory"]
    basis_file = opinf_model_directory * "/nn-opinf-basis.npz"
    basis = NPZ.npzread(basis_file)
    basis = basis["basis"]
    py"""
    import torch
    def get_model(model_file):
      return torch.load(model_file,weights_only=False)
    """
    ensemble_size = params["model"]["ensemble-size"]
    model = []
    for i in 1:ensemble_size
      tmp =  py"get_model"(opinf_model_directory * "/stiffness-" * string(i-1) * ".pt")
      push!(model,tmp)
    end
    num_dofs_per_node,num_nodes_basis,reduced_dim = size(basis)
    num_dofs = reduced_dim

    time = 0.0
    failed = false

    reduced_state = zeros(num_dofs)
    reduced_velocity = zeros(num_dofs)
    reduced_boundary_forcing = zeros(num_dofs)
    free_dofs = trues(num_dofs)
    boundary_conditions = Vector{BoundaryCondition}()
    NeuralNetworkOpInfRom(
        model,
        basis,
        reduced_state,
        reduced_velocity,
        reduced_boundary_forcing,
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
    )
end


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
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
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
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
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
        free_dofs,
        boundary_conditions,
        time,
        failed,
        fom_model,
        reference,
    )
end
