# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using NPZ

# --- Neural-network OpInf ROM: optional PyTorch backend ----------------------
# Loading the torch stiffness/BC models and running their forward passes is the
# only part of Norma that needs Python. That code lives in the NormaPyTorchExt
# package extension, which Julia loads automatically when PyCall is available in
# the active environment. The stub below is reached only when it is not, and
# gives a clear, actionable error instead of a cryptic failure.
const NN_OPINF_PYTORCH_MESSAGE =
    "Neural-network OpInf ROMs require the optional PyCall + PyTorch backend, " *
    "which is not installed. Add it with `import Pkg; Pkg.add(\"PyCall\")` and " *
    "make sure the Python it uses has PyTorch (torch) available."

# Load a torch model from each path and return them as a vector. The real
# implementation lives in ext/NormaPyTorchExt.jl; this broad `::Any` fallback is
# only used when the extension is not loaded.
load_torch_models(_paths) = norma_abort(NN_OPINF_PYTORCH_MESSAGE)

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
    ensemble_size = params["model"]["ensemble-size"]
    model = load_torch_models([
        opinf_model_directory * "/stiffness-" * string(i - 1) * ".pt" for i in 1:ensemble_size
    ])
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
