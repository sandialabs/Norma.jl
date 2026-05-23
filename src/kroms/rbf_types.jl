abstract type RBFKernel end
struct GaussianRBFKernel <: RBFKernel
    shape::Float64
end
struct InverseQuadraticRBFKernel <: RBFKernel
    shape::Float64
end
struct InverseMultiquadricRBFKernel <: RBFKernel
    shape::Float64
end
struct ThinPlateSplineRBFKernel <: RBFKernel
    shape::Float64
end
struct BasicMaternKernel <: RBFKernel
    shape::Float64
end
struct LinearMaternKernel <: RBFKernel
    shape::Float64
end
struct QuadraticMaternKernel <: RBFKernel
    shape::Float64
end

struct RBFKernelInterpolant
    rbf::RBFKernel
    coefficients::Matrix{Float64}
    training_data::Matrix{Float64}
end
