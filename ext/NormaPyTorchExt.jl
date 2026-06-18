# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

# Optional PyTorch backend for the neural-network OpInf ROM. Julia loads this
# extension automatically when PyCall is present in the active environment. It
# is the only part of Norma that depends on Python. All pure-Julia OpInf/ROM
# code (linear/quadratic/cubic OpInf and the RBF kROMs) lives in src/ and works
# without this extension.
module NormaPyTorchExt

using PyCall
using LinearAlgebra
import Norma:
    load_torch_models,
    apply_bc,
    apply_bc_detail,
    evaluate,
    NeuralNetworkOpInfRom,
    SolidMechanicsOpInfDirichletBC,
    SolidMechanicsOpInfOverlapSchwarzBoundaryCondition,
    RomNewmark,
    RomHessianMinimizer

# Convert a torch tensor to a Julia Vector{Float64}.
function torch_tensor_vector(tensor)
    values = tensor.detach().cpu().tolist()
    if ndims(values) == 1
        return Vector{Float64}(values)
    end
    return Vector{Float64}(values[1, :])
end

# Convert a torch tensor to a Julia Matrix{Float64}.
function torch_tensor_matrix(tensor)
    values = tensor.detach().cpu().tolist()
    if ndims(values) == 2
        return Matrix{Float64}(values)
    end
    return Matrix{Float64}(values[1, :, :])
end

# Real implementation of the stub declared in src/opinf/opinf_model.jl: load a
# torch model from each path and return them as a vector.
function load_torch_models(paths::Vector{String})
    py"""
    import torch, os
    def _norma_load_torch_model(model_file):
        assert os.path.isfile(model_file), model_file + " cannot be found"
        return torch.load(model_file, weights_only=False)
    """
    return Any[py"_norma_load_torch_model"(path) for path in paths]
end

function apply_bc(model::NeuralNetworkOpInfRom, bc::SolidMechanicsOpInfDirichletBC)
    model.fom_model.time = model.time
    apply_bc(model.fom_model, bc.fom_bc)
    bc_vector = zeros(0)
    for node_index in bc.fom_bc.node_set_node_indices
        disp_val = model.fom_model.displacement[bc.fom_bc.offset, node_index]
        push!(bc_vector, disp_val)
    end
    py"""
    import numpy as np
    import torch
    def setup_inputs(x):
        xi = np.zeros((1, x.size))
        xi[0] = x
        return torch.tensor(xi)
    """
    if bc.offset == 1
        offset_name = "x"
    end
    if bc.offset == 2
        offset_name = "y"
    end
    if bc.offset == 3
        offset_name = "z"
    end
    reduced_bc_vector = bc.basis[1, :, :]' * bc_vector
    model_inputs = py"setup_inputs"(reduced_bc_vector)
    name = "u-" * bc.name * "-" * offset_name
    jdict = PyDict(Dict(name => model_inputs))
    ensemble_size = size(bc.nn_model)[1]
    for i in 1:ensemble_size
        reduced_forcing = bc.nn_model[i].forward(jdict)
        reduced_forcing = torch_tensor_vector(reduced_forcing)
        model.reduced_boundary_forcing[:] += reduced_forcing
    end
    model.reduced_boundary_forcing[:] = model.reduced_boundary_forcing[:] ./ ensemble_size
end

function apply_bc_detail(model::NeuralNetworkOpInfRom, bc::SolidMechanicsOpInfOverlapSchwarzBoundaryCondition)
    if typeof(Norma.coupled_subsim_of(bc).model) == Norma.SolidMechanics
        ## Apply BC to the FOM vector
        apply_bc_detail(model.fom_model, bc.fom_bc)
        # populate our own BC vector
        unique_node_indices = unique(bc.fom_bc.side_set_node_indices)
        bc_vector = zeros(3, length(unique_node_indices))
        for i in 1:length(unique_node_indices)
            node_index = unique_node_indices[i]
            bc_vector[:, i] = model.fom_model.displacement[:, node_index]
        end
        py"""
        import numpy as np
        import torch
        def setup_inputs(x):
            xi = np.zeros((1, x.size))
            xi[0] = x
            inputs = torch.tensor(xi)
            return inputs
        """
        reduced_bc_vector = zeros(size(bc.basis)[3])
        for i in 1:3
            reduced_bc_vector[:] += bc.basis[i, :, :]' * bc_vector[i, :]
        end
        model_inputs = py"setup_inputs"(reduced_bc_vector)
        name = "u-" * bc.name
        jdict = PyDict(Dict(name => model_inputs))
        ensemble_size = size(bc.nn_model)[1]
        local_reduced_forcing = zeros(size(model.reduced_boundary_forcing)[1])
        for i in 1:ensemble_size
            reduced_forcing = bc.nn_model[i].forward(jdict)
            reduced_forcing = torch_tensor_vector(reduced_forcing)
            local_reduced_forcing[:] += reduced_forcing
        end
        local_reduced_forcing[:] = local_reduced_forcing[:] ./ ensemble_size
        model.reduced_boundary_forcing[:] += local_reduced_forcing
    else
        throw("ROM-ROM coupling not supported yet")
    end
end

function evaluate(integrator::RomNewmark, solver::RomHessianMinimizer, model::NeuralNetworkOpInfRom)
    beta = integrator.β
    dt = integrator.time_step
    num_dof, = size(model.free_dofs)
    I = Matrix{Float64}(LinearAlgebra.I, num_dof, num_dof)
    py"""
    import numpy as np
    import torch
    def setup_inputs(x):
        xi = np.zeros((1, x.size))
        xi[0] = x
        inputs = {'x': torch.tensor(xi)}
        return inputs
    """
    ensemble_size = size(model.nn_model)[1]
    stiffness = zeros(num_dof, num_dof)
    f = zeros(num_dof)
    model_inputs = py"setup_inputs"(solver.solution)
    for i in 1:ensemble_size
        ft, Kt = model.nn_model[i].forward(model_inputs, return_jacobian=true)
        f += torch_tensor_vector(ft)
        stiffness += torch_tensor_matrix(Kt)
    end
    stiffness = stiffness ./ ensemble_size
    f = f ./ ensemble_size
    LHS = I / (dt * dt * beta) - stiffness
    RHS = model.reduced_boundary_forcing + 1.0 / (dt * dt * beta) .* integrator.disp_pre
    residual = RHS - LHS * solver.solution
    solver.hessian[:, :] = LHS
    solver.gradient[:] = -residual
    if any(!isfinite, solver.gradient)
        model.failed = true
        Norma.norma_log(0, :error, "Non-finite values detected in ROM residual. This may indicate solution divergence.")
        return nothing
    end
    return nothing
end

end # module NormaPyTorchExt
