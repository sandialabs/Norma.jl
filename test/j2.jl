include("../src/constitutive_def.jl")
include("../src/constitutive.jl")

using LinearAlgebra

params = Dict{String,Any}()
params["elastic modulus"] = 200.0e+09
params["Poisson's ratio"] = 0.25
params["density"] = 7800.0
params["yield stress"] = 1.0e+09
material = J2(params)
F = [1.01 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
Fp = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
eqps = 0.0
dt = 1.0e-06
Fe, Fp, eqps, sigma = stress_update(material, F, Fp, eqps, dt)