# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

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
