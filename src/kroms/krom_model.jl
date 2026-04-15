using NPZ

function RBFKernelROM(params::Parameters)
    params["mesh smoothing"] = false
    fom_model = SolidMechanics(params)
    reference = fom_model.reference
    rom_model_file = params["model"]["model-file"]
    rom_model = NPZ.npzread(rom_model_file)

    interpolant = RBFKernelInterpolant(
        get_rbf_kernel(
            params["model"]["type"], 
            rom_model["rbf-shape-parameter"]
        ),
        Matrix{Float64}(rom_model["rbf-coefficients"]),
        Matrix{Float64}(rom_model["input-data"]),
    )
    forcing = rom_model["f"]
    
    basis = rom_model["basis"]
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
    return RBFKernelROM(
        interpolant,
        forcing, 
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

function get_rbf_kernel(model_type::String, shape::Float64)
    if model_type == "gaussian kernel rom"
        return GaussianRBFKernel(shape)
    elseif model_type == "inverse-quadratic kernel rom"
        return InverseQuadraticRBFKernel(shape)
    elseif model_type == "inverse-multiquadric kernel rom"
        return InverseMultiquadricRBFKernel(shape)
    elseif model_type == "thin-plate-spline kernel rom"
        return ThinPlateSplineRBFKernel(shape)
    elseif model_type == "basic-matern kernel rom"
        return BasicMaternKernel(shape)
    elseif model_type == "linear-matern kernel rom"
        return LinearMaternKernel(shape)
    elseif model_type == "quadatic-matern kernel rom"
        return QuadraticMaternKernel(shape)
    end
end
