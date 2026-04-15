
function compute_distances(
    X::Union{Matrix{Float64},Vector{Float64}},
    Y::Union{Matrix{Float64},Vector{Float64}},
)
    return sqrt.(compute_squared_distances(X, Y))
end

function compute_squared_distances(x::Vector{Float64}, y::Vector{Float64})
    return sum((x-y) .^ 2)
end

function compute_squared_distances(X::Matrix{Float64}, y::Vector{Float64})
    return [compute_squared_distances(collect(x), y) for x in eachcol(X)]
end

function compute_squared_distances(x::Vector{Float64}, Y::Matrix{Float64})
    return compute_squared_distances(Y, x)
end

function compute_squared_distances(X::Matrix{Float64}, Y::Matrix{Float64})
    return stack([compute_squared_distances(X, collect(y)) for y in eachcol(Y)])
end

function frobenius_norm(X::Matrix{Float64})
    return sqrt(sum(X .^ 2))
end

function test_compute_squared_distances()
    X = [1.0 2.0 3.0; 1.0 0.0 2.0]
    Y = [1.0 0.0; 1.0 2.0]
    true_dists = [0.0 2.0; 2.0 8.0; 5.0 9.0]

    computed_dists = compute_squared_distances(X, Y)
    error = frobenius_norm(true_dists-computed_dists)

    println("squared distance error = ", error)
    if error > 1e-16
        return false
    else
        return true
    end
end
