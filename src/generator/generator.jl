include("../types/types.jl")
using Distances, JLD2, Plots, Random
using .Types

function generate_instance(size::String, i::Int; write=true)
    K = 5
    M = 3
    B, S, P = 0
    if size == "S"
        B = 25
        S = 10
        P = 8
    elseif size == "M"
        B = 80
        S = 40
        P = 25
    elseif size == "L"
        B = 150
        S = 80
        P = 50
    elseif size == "XL"
        B = 300
        S = 120
        P = 70
    end
    BU_coords, S_coords = generate_coords(B, S)
    dist_mat = generate_dist(BU_coords, S_coords, B, S)
    parameters = generate_params(B, S)
    instance = Types.Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat, parameters...)
    dir_path = "instances/instances_"* size * "/"
    file_name = "inst_" * string(i) * ".jld2"
    full_path = dir_path * file_name
    if write
        Types.write_instance(instance, full_path)
    end
    return instance
end

function generate_coords(B, S)
    BU_coords = rand(10:800, (B, 2))
    S_coords = rand(10:800, (S, 2))
    return BU_coords, S_coords
end

function generate_dist(BU_coords, S_coords, B, S)
    metrica = Euclidean()
    mat = zeros(S, B)
    for i in 1:S
        for j in 1:B
            distancia = metrica(BU_coords[j, :], S_coords[i, :])
            mat[i, j] = distancia
        end
    end
    @debug "Wrote distance matrix"
    return trunc.(Int, mat)
end

function generate_params(B::Int64, S::Int64)
    Sk, Lk, Uk = generate_k(B, S)
    @debug "Wrote k"
    V, μ, T = generate_activities(B, S)
    @debug "Wrote activities"
    R, β = generate_risk(B, S)
    @debug "Wrote risk"
    return Sk, Lk, Uk, V, μ, T, R, β
end

function generate_k(B::Int64, S::Int64)
    K = 5
    Ks = rand(1:K, S)
    Sk = []
    for i in 1:K
        Si = findall(x->x==i, Ks)
        push!(Sk, Si)
    end
    Sk = convert(Vector{Vector{Int64}}, Sk)
    counts_k = [length(s) for s in Sk]
    Lk_flt = [k - 0.4k for k in counts_k]
    Uk_flt = [k + 0.4k for k in counts_k]
    Lk = trunc.(Int, Lk_flt)
    Uk = trunc.(Int, Uk_flt)
    return Sk, Lk, Uk
end

function generate_activities(B::Int64, S::Int64)
    M = 3
    V = [zeros(Int64, B) for _ in 1:M]
    for i in 1:length(V)
        V[i] = rand(10:30, B)
    end
    μ = [zeros(Int64, S) for _ in 1:M]
    for i in 1:length(μ)
        sum_vals = sum(V[i])
        μ[i] = fill(trunc(Int, sum_vals/S), S)
    end
    T = fill(0.3, M)
    return V, μ, T
end

function generate_risk(B::Int64, S::Int64)
    R = rand(2:15, B)
    sum_R = sum(R)
    lower = trunc(Int, trunc(sum_R/S) - 0.01sum_R)
    upper = trunc(Int, trunc(sum_R/S) + 0.01sum_R)
    β = rand(lower:upper, S)
    return R, β
end

function main()
    size = ARGS[1], n = Int(ARGS[2])
    for i in 1:n
        generate_instance(size, i)
    end
    @info "Done"
end
