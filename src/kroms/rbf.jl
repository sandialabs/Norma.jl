# using LinearAlgebra
# using Printf

include("rbf_types.jl")
include("helpers.jl")

function evaluate_rbf_interpolant(interpolant::RBFKernelInterpolant, y::Vector{Float64})
    val, jac = evaluate_kernel_with_jacobian(interpolant.rbf, interpolant.training_data, y)
    return -interpolant.coefficients*val, -interpolant.coefficients*jac
end

function evaluate_kernel_with_jacobian(
    rbf::GaussianRBFKernel,
    X::Matrix{Float64},
    y::Vector{Float64},
)
    val = exp.(-rbf.shape^2*compute_squared_distances(X, y))
    deriv_over_r = -2.0*val
    jac = compute_rbf_jacobian(rbf.shape, deriv_over_r, X, y)
    return val, jac
end

function evaluate_kernel_with_jacobian(
    rbf::InverseQuadraticRBFKernel,
    X::Matrix{Float64},
    y::Vector{Float64},
)
    val = 1.0 ./ (1.0 .+ rbf.shape^2*compute_squared_distances(X, y))
    deriv_over_r = -2.0 * val .^ 2
    jac = compute_rbf_jacobian(rbf.shape, deriv_over_r, X, y)
    return val, jac
end

function evaluate_kernel_with_jacobian(
    rbf::InverseMultiquadricRBFKernel,
    X::Matrix{Float64},
    y::Vector{Float64},
)
    val = 1.0 ./ sqrt.(1.0 .+ rbf.shape^2*compute_squared_distances(X, y))
    deriv_over_r = -1.0*val .^ 3
    jac = compute_rbf_jacobian(rbf.shape, deriv_over_r, X, y)
    return val, jac
end

function evaluate_kernel_with_jacobian(
    rbf::ThinPlateSplineRBFKernel,
    X::Matrix{Float64},
    y::Vector{Float64},
)
    shapeTimesDists2 = rbf.shape^2*compute_squared_distances(X, y)
    logr = log.(sqrt.(shapeTimesDists2))

    ind = shapeTimesDists2 .!= 0.0
    val = zeros(length(shapeTimesDists2))
    val[ind] += shapeTimesDists2[ind] .* logr[ind]

    deriv_over_r = 2.0*logr .+ 1.0
    jac = compute_rbf_jacobian(rbf.shape, deriv_over_r, X, y)
    return val, jac
end

function evaluate_kernel_with_jacobian(
    rbf::BasicMaternKernel,
    X::Matrix{Float64},
    y::Vector{Float64},
)
    shapeTimesDists = rbf.shape*compute_distances(X, y)
    val = exp.(-shapeTimesDists)

    shapeTimesDists[shapeTimesDists .== 0.0].+=1e-16
    deriv_over_r = -val ./ shapeTimesDists
    jac = compute_rbf_jacobian(rbf.shape, deriv_over_r, X, y)

    return val, jac
end

function evaluate_kernel_with_jacobian(
    rbf::LinearMaternKernel,
    X::Matrix{Float64},
    y::Vector{Float64},
)
    shapeTimesDists = rbf.shape*compute_distances(X, y)
    expMinusr = exp.(-shapeTimesDists)
    val = (1.0 .+ shapeTimesDists) .* expMinusr
    deriv_over_r = -expMinusr
    jac = compute_rbf_jacobian(rbf.shape, deriv_over_r, X, y)
    return val, jac
end

function evaluate_kernel_with_jacobian(
    rbf::QuadraticMaternKernel,
    X::Matrix{Float64},
    y::Vector{Float64},
)
    shapeTimesDists2 = rbf.shape^2*compute_squared_distances(X, y)
    shapeTimesDists = sqrt.(shapeTimesDists2)
    expMinusr = exp.(-shapeTimesDists)
    onePlusShapeTimesDists = 1.0 .+ shapeTimesDists

    val = (3.0*onePlusShapeTimesDists+shapeTimesDists2) .* expMinusr
    deriv_over_r = -onePlusShapeTimesDists .* expMinusr

    jac = compute_rbf_jacobian(rbf.shape, deriv_over_r, X, y)
    return val, jac
end

function compute_rbf_jacobian(
    shape::Float64,
    deriv_over_r::Vector{Float64},
    X::Matrix{Float64},
    y::Vector{Float64},
)
    return (-shape^2*deriv_over_r) .* transpose(X .- y)
end


# Testing ------------------------------------------------------------
function fd_jac(f, x::Vector{Float64}, h::Float64)
    fx = f(x)
    nf = size(fx)[1]
    nx = size(x)[1]
    jac = zeros(nf, nx)
    for i = 1:nx
        ei = zeros(nx)
        ei[i] += h
        jac[:, i] += (f(x+ei)-fx)/h
    end
    return jac
end

function test_rbf()
    coefficients = rand(2, 4)
    training_data = rand(3, 4)
    shape = 0.5
    x = rand(3)

    rbf_list = [
        GaussianRBFKernel,
        InverseQuadraticRBFKernel,
        InverseMultiquadricRBFKernel,
        ThinPlateSplineRBFKernel,
        BasicMaternKernel,
        LinearMaternKernel,
        QuadraticMaternKernel,
    ]

    for rbf_type in rbf_list
        rbf = rbf_type(shape)
        println("RBF type = ", typeof(rbf))

        rbf = RBFKernelInterpolant(rbf, coefficients, training_data)
        val, true_jac = evaluate_rbf_interpolant(rbf, x)

        f(y) = evaluate_rbf_interpolant(rbf, y)[1]
        _test_fd_jac_for_decreasing_h(f, x, true_jac)
    end
end

function _test_fd_jac_for_decreasing_h(f, x, true_jac)
    println("   h           error")
    println("----------------------")
    for i = 1:8
        h = 10.0^(-1.0*i)
        approx_jac = fd_jac(f, x, h)
        error = frobenius_norm(approx_jac-true_jac)

        @printf("%1.1e      %1.3e\n", h, error)
    end
    println()
end


function test_fd_jac()
    println("f(x) = x^2")
    f(x) = x .^ 2
    dfdx(x) = diagm(2*x)
    x = ones(2)

    _test_fd_jac_for_decreasing_h(f, x, dfdx(x))
end

# test_rbf()
# test_fd_jac()
