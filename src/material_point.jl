# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# Material point driver: integrates a material model along a prescribed
# history at a single point, Norma's counterpart of LCM's Material Point
# Simulator. Each component of the deformation gradient is either
# prescribed, F_ij(t), or free, with the matching component of the first
# Piola-Kirchhoff stress prescribed instead, P_ij(t) (zero unless given).
# Uniaxial stress along x, for example, prescribes F11 and the components
# F21, F31, F32 that remove the rigid rotation, and leaves the other five
# free with zero stress.
#
# The Taylor and Sachs aggregates of a polycrystal build on the same steps:
# the Taylor (iso-strain) estimate gives every grain the macroscopic
# deformation gradient, and the Sachs (iso-stress) estimate gives every grain
# the macroscopic stress.

struct MaterialPointHistory
    time::Vector{Float64}
    F::Vector{SMatrix{3,3,Float64,9}}
    P::Vector{SMatrix{3,3,Float64,9}}
    σ::Vector{SMatrix{3,3,Float64,9}}  # Cauchy stress
    W::Vector{Float64}
    state::Vector{Vector{Float64}}
end

function cauchy_matrix(voigt::SVector{6,Float64})
    return @SMatrix [
        voigt[1] voigt[6] voigt[5]
        voigt[6] voigt[2] voigt[4]
        voigt[5] voigt[4] voigt[3]
    ]
end

function cauchy_stress(material::Solid, P::SMatrix{3,3,Float64,9}, F::SMatrix{3,3,Float64,9})
    return cauchy_matrix(voigt_cauchy_from_stress(material, P, F, det(F)))
end

# The components of F marked free, in column-major order.
free_components(free::AbstractMatrix{Bool}) = [CartesianIndex(i, j) for j in 1:3 for i in 1:3 if free[i, j]]

# The matrix of ∂P_a/∂F_b over the free components a, b.
function free_jacobian(AA::SArray{Tuple{3,3,3,3},Float64,4,81}, free_indices::Vector{CartesianIndex{2}})
    return [AA[a[1], a[2], b[1], b[2]] for a in free_indices, b in free_indices]
end

# F with `δ` added to its free components.
function add_free(F::SMatrix{3,3,Float64,9}, free_indices::Vector{CartesianIndex{2}}, δ::AbstractVector)
    G = MMatrix{3,3,Float64,9}(F)
    for (k, c) in enumerate(free_indices)
        G[c] += δ[k]
    end
    return SMatrix(G)
end

# A Newton correction that a material cannot take, beyond the range of a
# crystal plasticity flow rule for example, is halved up to this many times.
const MAXIMUM_CORRECTION_HALVINGS = 10

error_message(error::Exception) = error isa DomainError ? error.msg : sprint(showerror, error)

# `f(α)` for α = 1, 1/2, 1/4, ... until it returns without a math error;
# returns α and the value, or nothing when every halving failed.
function halve_until_admissible(f::Function)
    α = 1.0
    for _ in 0:MAXIMUM_CORRECTION_HALVINGS
        value = try
            f(α)
        catch error
            error isa _MATH_ERRORS || rethrow()
            nothing
        end
        value === nothing || return α, value
        α /= 2
    end
    return nothing
end

# One step of a material point from `state_old`: the free components of F,
# starting from those of F_start, are solved by Newton's method so that the
# matching components of P equal those of `target`. Returns F, W, P, the
# tangent AA and the new state.
function material_point_step(
    material::Solid,
    F_start::SMatrix{3,3,Float64,9},
    orientation::SMatrix{3,3,Float64,9},
    state_old::Vector{Float64},
    Δt::Float64,
    free_indices::Vector{CartesianIndex{2}},
    target::SMatrix{3,3,Float64,9};
    tolerance::Float64=1.0e-12,
    maximum_iterations::Int=50,
)
    F = F_start
    need_tangent = !isempty(free_indices)
    W, P, AA, state_new = constitutive_update(material, F, orientation, state_old, Δt; need_tangent=need_tangent)
    converged = isempty(free_indices)
    for _ in 1:maximum_iterations
        converged && break
        residual = [P[c] - target[c] for c in free_indices]
        if norm(residual) <= tolerance * max(1.0, norm(P))
            converged = true
            break
        end
        correction = -(free_jacobian(AA, free_indices) \ residual)
        trial = halve_until_admissible() do α
            F_trial = add_free(F, free_indices, α .* correction)
            return F_trial, constitutive_update(material, F_trial, orientation, state_old, Δt; need_tangent=true)
        end
        trial === nothing && throw(
            DomainError(norm(correction), "Material point: the material cannot take the stress-controlled correction."),
        )
        α, (F, (W, P, AA, state_new)) = trial
        # F is dimensionless: a correction at round-off ends the iteration
        # even where the stress residual cannot drop below the round-off of a
        # stiff material.
        if α * norm(correction) <= 1.0e-14
            converged = true
            break
        end
    end
    converged || throw(DomainError(0.0, "Material point: the stress-controlled components did not converge."))
    return F, W, P, AA, state_new
end

# Response of `material` along `deformation(t)` at `times`. The first time is
# the initial state, evaluated with a zero time increment. `free[i, j]` marks
# the components of F solved for, so that P_ij = stress(t)[i, j]; those
# components of deformation(t) only serve as the initial guess at the first
# time. `orientation` is the lattice orientation Q of oriented materials.
function material_point_history(
    material::Solid,
    deformation::Function,
    times::AbstractVector{<:Real};
    orientation::SMatrix{3,3,Float64,9}=SMatrix{3,3,Float64,9}(I),
    free::AbstractMatrix{Bool}=falses(3, 3),
    stress::Function=t -> zero(SMatrix{3,3,Float64,9}),
    tolerance::Float64=1.0e-12,
    maximum_iterations::Int=50,
)
    free_indices = free_components(free)
    state = number_states(material) > 0 ? initial_state(material) : Float64[]
    history = MaterialPointHistory(
        Float64[],
        SMatrix{3,3,Float64,9}[],
        SMatrix{3,3,Float64,9}[],
        SMatrix{3,3,Float64,9}[],
        Float64[],
        Vector{Float64}[],
    )
    F_previous = SMatrix{3,3,Float64,9}(deformation(times[1]))
    for (n, t) in enumerate(times)
        Δt = n == 1 ? 0.0 : Float64(t - times[n - 1])
        F_start = MMatrix{3,3,Float64,9}(deformation(t))
        for c in free_indices
            F_start[c] = F_previous[c]
        end
        F, W, P, _, state = try
            material_point_step(
                material,
                SMatrix(F_start),
                orientation,
                state,
                Δt,
                free_indices,
                SMatrix{3,3,Float64,9}(stress(t));
                tolerance=tolerance,
                maximum_iterations=maximum_iterations,
            )
        catch error
            error isa _MATH_ERRORS || rethrow()
            # The material itself can fail on the first evaluation of the step,
            # outside the halving of the corrections.
            error isa DomainError && occursin("Material point", error.msg) || norma_abortf(
                "Material point: the material cannot take the step at t = %.6e: %s", t, error_message(error)
            )
            norma_abortf("Material point: the stress-controlled components did not converge at t = %.6e.", t)
        end
        push!(history.time, Float64(t))
        push!(history.F, F)
        push!(history.P, P)
        push!(history.σ, cauchy_stress(material, P, F))
        push!(history.W, W)
        push!(history.state, copy(state))
        F_previous = F
    end
    return history
end

# Response of a polycrystal estimated from its grains: the time, the
# macroscopic deformation gradient, the volume averages of P and of the
# Cauchy stress, and the deformation gradient and state of every grain.
struct AggregateHistory
    time::Vector{Float64}
    F::Vector{SMatrix{3,3,Float64,9}}
    P::Vector{SMatrix{3,3,Float64,9}}
    σ::Vector{SMatrix{3,3,Float64,9}}
    grain_F::Vector{Vector{SMatrix{3,3,Float64,9}}}
    grain_state::Vector{Vector{Vector{Float64}}}
end

function AggregateHistory()
    return AggregateHistory(
        Float64[],
        SMatrix{3,3,Float64,9}[],
        SMatrix{3,3,Float64,9}[],
        SMatrix{3,3,Float64,9}[],
        Vector{SMatrix{3,3,Float64,9}}[],
        Vector{Vector{Float64}}[],
    )
end

function check_aggregate(materials::AbstractVector, orientations::AbstractVector, weights::AbstractVector)
    length(materials) == length(orientations) == length(weights) ||
        norma_abort("Aggregate: one material, orientation and weight per grain are required.")
    all(>=(0), weights) && sum(weights) > 0 || norma_abort("Aggregate: the grain weights must be nonnegative.")
    return weights ./ sum(weights)
end

initial_grain_state(material::Solid) = number_states(material) > 0 ? initial_state(material) : Float64[]

# Taylor (iso-strain) estimate: every grain g, of material materials[g],
# orientation orientations[g] and volume fraction weights[g], has the
# macroscopic deformation gradient deformation(t). Its free components are
# solved so that the volume average of P meets stress(t), as in
# material_point_history.
function taylor_history(
    materials::AbstractVector,
    orientations::AbstractVector,
    weights::AbstractVector,
    deformation::Function,
    times::AbstractVector{<:Real};
    free::AbstractMatrix{Bool}=falses(3, 3),
    stress::Function=t -> zero(SMatrix{3,3,Float64,9}),
    tolerance::Float64=1.0e-12,
    maximum_iterations::Int=50,
)
    w = check_aggregate(materials, orientations, weights)
    free_indices = free_components(free)
    need_tangent = !isempty(free_indices)
    states = [initial_grain_state(m) for m in materials]
    history = AggregateHistory()
    F_previous = SMatrix{3,3,Float64,9}(deformation(times[1]))
    for (n, t) in enumerate(times)
        Δt = n == 1 ? 0.0 : Float64(t - times[n - 1])
        responses_at(F) = [
            constitutive_update(materials[g], F, orientations[g], states[g], Δt; need_tangent=need_tangent) for
            g in eachindex(materials)
        ]
        F = MMatrix{3,3,Float64,9}(deformation(t))
        for c in free_indices
            F[c] = F_previous[c]
        end
        F = SMatrix(F)
        target = SMatrix{3,3,Float64,9}(stress(t))
        responses = try
            responses_at(F)
        catch error
            error isa _MATH_ERRORS || rethrow()
            norma_abortf(
                "Taylor aggregate: a grain cannot take the deformation at t = %.6e: %s", t, error_message(error)
            )
        end
        converged = !need_tangent
        for _ in 1:maximum_iterations
            converged && break
            P = sum(w[g] * responses[g][2] for g in eachindex(materials))
            residual = [P[c] - target[c] for c in free_indices]
            if norm(residual) <= tolerance * max(1.0, norm(P))
                converged = true
                break
            end
            AA = sum(w[g] * responses[g][3] for g in eachindex(materials))
            correction = -(free_jacobian(AA, free_indices) \ residual)
            trial = halve_until_admissible() do α
                F_trial = add_free(F, free_indices, α .* correction)
                return F_trial, responses_at(F_trial)
            end
            trial === nothing && norma_abortf("Taylor aggregate: a grain cannot take the correction at t = %.6e.", t)
            α, (F, responses) = trial
            if α * norm(correction) <= 1.0e-14
                converged = true
                break
            end
        end
        converged || norma_abortf("Taylor aggregate: the free components did not converge at t = %.6e.", t)
        states = [responses[g][4] for g in eachindex(materials)]
        push!(history.time, Float64(t))
        push!(history.F, F)
        push!(history.P, sum(w[g] * responses[g][2] for g in eachindex(materials)))
        push!(history.σ, sum(w[g] * cauchy_stress(materials[g], responses[g][2], F) for g in eachindex(materials)))
        push!(history.grain_F, fill(F, length(materials)))
        push!(history.grain_state, copy.(states))
        F_previous = F
    end
    return history
end

# Sachs (iso-stress) estimate in uniaxial stress along x: every grain carries
# the same first Piola-Kirchhoff stress, P = s e1 ⊗ e1, with its rigid
# rotation fixed by F21 = F31 = F32 = 0, and s is solved at every time so that
# the volume average of F11 - 1 equals axial_strain(t). The iteration on s
# ends when its Newton correction is within the stress tolerance.
function sachs_history(
    materials::AbstractVector,
    orientations::AbstractVector,
    weights::AbstractVector,
    axial_strain::Function,
    times::AbstractVector{<:Real};
    tolerance::Float64=1.0e-12,
    maximum_iterations::Int=50,
)
    w = check_aggregate(materials, orientations, weights)
    free_indices = free_components([true true true; false true true; false false true])
    k11 = findfirst(==(CartesianIndex(1, 1)), free_indices)
    num_grains = length(materials)
    states = [initial_grain_state(m) for m in materials]
    grain_F = fill(SMatrix{3,3,Float64,9}(I), num_grains)
    history = AggregateHistory()
    s = 0.0
    for (n, t) in enumerate(times)
        Δt = n == 1 ? 0.0 : Float64(t - times[n - 1])
        ε = Float64(axial_strain(t))
        steps_at(s) = [
            material_point_step(
                materials[g],
                grain_F[g],
                orientations[g],
                states[g],
                Δt,
                free_indices,
                SMatrix{3,3,Float64,9}(s, 0, 0, 0, 0, 0, 0, 0, 0);
                tolerance=tolerance,
            ) for g in 1:num_grains
        ]
        steps = try
            steps_at(s)
        catch error
            error isa _MATH_ERRORS || rethrow()
            norma_abortf("Sachs aggregate: a grain cannot take the stress at t = %.6e: %s", t, error_message(error))
        end
        converged = false
        for _ in 1:maximum_iterations
            mean_strain = sum(w[g] * (steps[g][1][1, 1] - 1) for g in 1:num_grains)
            # The axial compliance of each grain under the free components.
            compliance = sum(w[g] * inv(free_jacobian(steps[g][4], free_indices))[k11, k11] for g in 1:num_grains)
            correction = (ε - mean_strain) / compliance
            if abs(correction) <= tolerance * max(1.0, abs(s))
                converged = true
                break
            end
            trial = halve_until_admissible(α -> steps_at(s + α * correction))
            trial === nothing && norma_abortf("Sachs aggregate: a grain cannot take the stress at t = %.6e.", t)
            α, steps = trial
            s += α * correction
        end
        converged || norma_abortf("Sachs aggregate: the axial stress did not converge at t = %.6e.", t)
        grain_F = [steps[g][1] for g in 1:num_grains]
        states = [steps[g][5] for g in 1:num_grains]
        push!(history.time, Float64(t))
        push!(history.F, sum(w[g] * grain_F[g] for g in 1:num_grains))
        push!(history.P, sum(w[g] * steps[g][3] for g in 1:num_grains))
        push!(history.σ, sum(w[g] * cauchy_stress(materials[g], steps[g][3], grain_F[g]) for g in 1:num_grains))
        push!(history.grain_F, copy(grain_F))
        push!(history.grain_state, copy.(states))
    end
    return history
end
