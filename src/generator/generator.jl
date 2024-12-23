using Distances, JLD2, Random
using Types

function generate_instance(size::String, i::Int; write=true)
    K = 5
    M = 3
    B = 0
    S = 0
    P = 0
    if size == "S"
        B = 15
        S = 8
        P = 4
    elseif size == "M"
        B = 200
        S = 40
        P = 20
    elseif size == "L"
        B = 1500
        S = 120
        P = 100
    elseif size == "XL"
        B = 500
        S = 80
        P = 50
    end
    BU_coords, S_coords = generate_coords(B, S)
    dist_mat = generate_dist(BU_coords, S_coords, B, S)
    parameters = generate_params(B, S, P)
    instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat, parameters...)
    dir_path = "instances/instances_" * size * "/"
    if !isdir(dir_path)
        mkdir(dir_path)
    end
    file_name = "inst_" * string(i) * ".jld2"
    plot_name = replace(file_name, ".jld2" => ".png")
    full_path = dir_path * file_name
    full_plot = dir_path * plot_name
    if write
        write_instance(instance, full_path)
        plot_instance(instance, full_plot)
    end
    return instance
end

function generate_coords(B, S)
    BU_coords = rand(5:10000, (B, 2))
    S_coords = rand(5:10000, (S, 2))
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

function generate_params(B::Int64, S::Int64, P)
    Sk, Lk, Uk = generate_k(B, S)
    @debug "Wrote k"
    V, μ, T = generate_activities(B, S, P)
    @debug "Wrote activities"
    R, β = generate_risk(B, S, P)
    @debug "Wrote risk"
    return Sk, Lk, Uk, V, μ, T, R, β
end

function generate_k(B::Int64, S::Int64)
    K = 5
    Ks = rand(1:K, S)
    Sk = []
    for i in 1:K
        Si = findall(x -> x == i, Ks)
        push!(Sk, Si)
    end
    Sk = convert(Vector{Vector{Int64}}, Sk)
    counts_k = [length(s) for s in Sk]
    Lk_flt = [k - 0.85k for k in counts_k]
    Uk_flt = [k - 0.05k for k in counts_k]
    Lk = trunc.(Int, Lk_flt)
    Uk = trunc.(Int, Uk_flt)
    return Sk, Lk, Uk
end

function generate_activities(B::Int64, S::Int64, P)
    M = 3
    V = [zeros(Int64, B) for _ in 1:M]
    for i in eachindex(V)
        V[i] = rand(10:30, B)
    end
    μ = [zeros(Int64, S) for _ in 1:M]
    for i in eachindex(V)
        sum_vals = sum(V[i])
        τ = 0.6
        upper = trunc(Int, (trunc(Int, (sum_vals / P)) * (1 + τ)))
        #println("upper $i $upper")
        μ[i] = fill(upper, S)
    end
    T = fill(0.4, M)
    return V, μ, T
end

function generate_risk(B::Int64, S::Int64, P)
    R = rand(40:60, B)
    sum_R = sum(R)
    τ = 0.8
    upper = trunc(Int, trunc(Int, (sum_R / P) * (1 + τ)))
    #println("beta: $upper")
    β = fill(upper, S)
    return R, β
end

function main()
    size = ARGS[1]
    n = parse(Int64, ARGS[2])
    for i in 1:n
        generate_instance(size, i)
    end
    @info "Done"
end

#main()
