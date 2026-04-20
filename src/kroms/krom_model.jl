using NPZ

function RBFKernelROM(params::Parameters)
    params["mesh smoothing"] = false
    fom_model = SolidMechanics(params)
    reference = fom_model.reference
    rom_model_file = params["model"]["model-file"]
    rom_model = NPZ.npzread(rom_model_file)

    interpolant = RBFKernelInterpolant(
        get_rbf_kernel(params["model"]["rbf-type"], rom_model["rbf-shape-parameter"]),
        Matrix{Float64}(rom_model["rbf-coefficients"]),
        Matrix{Float64}(rom_model["input-data"]),
    )

    basis = rom_model["basis"]
    _, _, reduced_dim = size(basis)
    num_dofs = reduced_dim
    time = 0.0
    failed = false
    forcing = zeros(num_dofs)
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

function get_rbf_kernel(rbf_type::String, shape::Float64)
    if rbf_type == "gaussian"
        return GaussianRBFKernel(shape)
    elseif rbf_type == "inverse-quadratic"
        return InverseQuadraticRBFKernel(shape)
    elseif rbf_type == "inverse-multiquadric"
        return InverseMultiquadricRBFKernel(shape)
    elseif rbf_type == "thin-plate-spline"
        return ThinPlateSplineRBFKernel(shape)
    elseif rbf_type == "basic-matern"
        return BasicMaternKernel(shape)
    elseif rbf_type == "linear-matern"
        return LinearMaternKernel(shape)
    elseif rbf_type == "quadatic-matern"
        return QuadraticMaternKernel(shape)
    end
end
