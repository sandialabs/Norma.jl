# Norma.jl 1.0: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
include("../src/constitutive_types.jl")
include("../src/constitutive.jl")

using LinearAlgebra

params = Parameters()
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