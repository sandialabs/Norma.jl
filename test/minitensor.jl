# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.

using Test
using LinearAlgebra
using StaticArrays

@testset "Minitensor" begin
    @testset "Skew And Symm" begin
        A = randn(3, 3)
        S = Norma.skew(A)
        Sym = Norma.symm(A)
        @test isapprox(S + S', zeros(3, 3); atol=1e-14)
        @test isapprox(Sym - Sym', zeros(3, 3); atol=1e-14)
        @test isapprox(A, S + Sym; atol=1e-14)
    end

    @testset "Hat And Vee" begin
        v = @SVector randn(3)
        W = Norma.hat(v)
        v2 = Norma.vee(W)
        @test isapprox(v2, v; atol=1e-14)
        w = @SVector randn(3)
        @test isapprox(W * w, cross(v, w); atol=1e-14)
    end

    @testset "Quaternion Round Trip" begin
        w = @SVector [0.1, -0.2, 0.3]
        R = Norma.rt_of_rv(w)
        w2 = Norma.rv_of_rt(R)
        @test isapprox(w2, w; atol=1e-12)

        q = Norma.q_of_rv(w)
        w3 = Norma.rv_of_q(q)
        @test isapprox(w3, w; atol=1e-12)

        R2 = Norma.rt_of_q(q)
        @test isapprox(R2, R; atol=1e-12)

        q2 = Norma.q_of_rt(R)
        @test isapprox(q2, q; atol=1e-12) || isapprox(q2, -q; atol=1e-12)
    end

    @testset "Rv Continue" begin
        w = @SVector [π, 0.0, 0.0]
        w_prev = @SVector [-π, 0.0, 0.0]
        w_cont = Norma.rv_continue(w, w_prev)
        @test isapprox(Norma.rt_of_rv(w_cont), Norma.rt_of_rv(w); atol=1e-12)
        @test dot(w_cont, w_prev) < 0.0  # same direction as previous
    end

    @testset "BCH Identity (approximate Exp(log(exp(x)exp(y))))" begin
        xv = @SVector [π / 2, 0, 0]
        yv = @SVector [0, π / 2, 0]
        x = Norma.hat(xv)
        y = Norma.hat(yv)
        X = exp(x)
        Y = exp(y)
        Z = X * Y
        w = Norma.bch(x, y)
        W = exp(w)
        @test isapprox(W, Z; atol=0.12)
    end
end
