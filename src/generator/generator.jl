include("../types/types.jl")
using Distances, JLD2, Plots, Random, .Types

function write_jld2(size::String)
    K = 5
    M = 3
    if size == "S"
        B = 8
        S = 5
        P = 3
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        parameters = generate_params(size, B, S)
        instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat,parameters...)
        dir_path = "instances/"
        file_name = "inst_" * size * ".jld2"
        full_path = dir_path * file_name
        jldsave(full_path; instance)
    elseif size == "M"
        B = 80
        S = 40
        P = 25
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        parameters = generate_params(size, B, S)
        instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat,parameters...)
        dir_path = "instances/"
        file_name = "inst_" * size * ".jld2"
        full_path = dir_path * file_name
        jldsave(full_path; instance)
    elseif size == "L"
        B = 150
        S = 80
        P = 50
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        parameters = generate_params(size, B, S)
        instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat,parameters...)
        dir_path = "instances/"
        file_name = "inst_" * size * ".jld2"
        full_path = dir_path * file_name
        jldsave(full_path; instance)
    elseif size == "XL"
        B = 300
        S = 120
        P = 70
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        parameters = generate_params(size, B, S)
        instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat,parameters...)
        dir_path = "instances/"
        file_name = "inst_" * size * ".jld2"
        full_path = dir_path * file_name
        jldsave(full_path; instance)
    else
        @error "Size not recognized"
        exit(1)
    end

end


function generate_coords(B, S)
    BU_coords = rand(10:500, (B, 2))
    S_coords = rand(10:500, (S, 2))

    Plots.scatter(
        BU_coords[:, 1],
        BU_coords[:, 2],
        label = "BUs",
        markershape = :circle,
        markercolor = :blue,
    )
    Plots.scatter!(
        S_coords[:, 1],
        S_coords[:, 2],
        label = "Branches",
        markershape = :square,
        markercolor = :white,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    png("out/plot_instance_bu$B" * "_suc" * "$S" * ".png")
    #@info "Wrote plot and coords"
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
    @info "Wrote distance matrix"
    return trunc.(Int, mat)
end

function generate_params(size::String, B::Int64, S::Int64)
    if size == "S"
        Sk, Lk, Uk = generate_k(size, B, S)
        #@info "Wrote k"
        V, μ, T = generate_activities(size, B, S)
        #@info "Wrote activities"
        R, β = generate_risk(size, B, S)
        #@info "Wrote risk"
        return Sk, Lk, Uk, V, μ, T, R, β
    elseif size == "M"
        Sk, Lk, Uk = generate_k(size, B, S)
        @info "Wrote k"
        V, μ, T = generate_activities(size, B, S)
        @info "Wrote activities"
        R, β = generate_risk(size, B, S)
        @info "Wrote risk"
        return Sk, Lk, Uk, V, μ, T, R, β
    elseif size == "L"
        Sk, Lk, Uk = generate_k(size, B, S)
        @info "Wrote k"
        V, μ, T = generate_activities(size, B, S)
        @info "Wrote activities"
        R, β = generate_risk(size, B, S)
        @info "Wrote risk"
        return Sk, Lk, Uk, V, μ, T, R, β
    elseif size == "XL"
        Sk, Lk, Uk = generate_k(size, B, S)
        @info "Wrote k"
        V, μ, T = generate_activities(size, B, S)
        @info "Wrote activities"
        R, β = generate_risk(size, B, S)
        @info "Wrote risk"
        return Sk, Lk, Uk, V, μ, T, R, β
    end
end

function generate_k(size::String, B::Int64, S::Int64)
    if size == "S"
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
    elseif size == "M"
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
    elseif size == "L"
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
    elseif size == "XL"
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
end

function generate_activities(size::String, B::Int64, S::Int64)
    M = 3
    if size == "S"
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
    elseif size == "M"
        P = 25
        V = [Array{Int64}(undef, B) for _ in 1:M]
        for v in V
            v = rand(10:20, B)
        end
        μ = [Array{Int64}(undef, S) for _ in 1:M]
        for (index, μₘ) in enumerate(μ)
            sum_vals = sum(V[index])
            deviation = trunc(Int, sum_vals/10)
            lower = trunc(sum_vals/S) - deviation
            upper = trunc(sum_vals/S) + deviation
            μₘ = rand(lower:upper, S)
        end
        T = fill(0.1, M)
        return V, μ, T
    elseif size == "L"
        P = 50
        V = [Array{Int64}(undef, B) for _ in 1:M]
        for v in V
            v = rand(10:20, B)
        end
        μ = [Array{Int64}(undef, S) for _ in 1:M]
        for (index, μₘ) in enumerate(μ)
            sum_vals = sum(V[index])
            deviation = trunc(Int, sum_vals/10)
            lower = trunc(sum_vals/S) - deviation
            upper = trunc(sum_vals/S) + deviation
            μₘ = rand(lower:upper, S)
        end
        T = fill(0.1, M)
        return V, μ, T
    elseif size == "XL"
        P = 70
        V = [Array{Int64}(undef, B) for _ in 1:M]
        for v in V
            v = rand(10:20, B)
        end
        μ = [Array{Int64}(undef, S) for _ in 1:M]
        for (index, μₘ) in enumerate(μ)
            sum_vals = sum(V[index])
            deviation = trunc(Int, sum_vals/10)
            lower = trunc(sum_vals/S) - deviation
            upper = trunc(sum_vals/S) + deviation
            μₘ = rand(lower:upper, S)
        end
        T = fill(0.1, M)
        return V, μ, T
    end
end

function generate_risk(size::String, B::Int64, S::Int64)
    if size == "S"
        R = rand(2:10, B)
        sum_R = sum(R)
        lower = trunc(Int, trunc(sum_R/S) - 0.01sum_R)
        upper = trunc(Int, trunc(sum_R/S) + 0.01sum_R)
        β = rand(lower:upper, S)
        return R, β
    elseif size == "M"
        R = rand(10:20, B)
        sum_R = sum(R)
        deviation = trunc(Int, (sum_R/80))
        lower = trunc(sum_R/S) - deviation
        upper = trunc(sum_R/S) + deviation
        β = rand(lower:upper, S)
        return R, β
    elseif size == "L"
        R = rand(15:25, B)
        sum_R = sum(R)
        deviation = trunc(Int, sum_R/80)
        lower = trunc(sum_R/S) - deviation
        upper = trunc(sum_R/S) + deviation
        β = rand(lower:upper, S)
        return R, β
    elseif size == "XL"
        R = rand(20:35, B)
        sum_R = sum(R)
        deviation = trunc(Int, sum_R/80)
        lower = trunc(sum_R/S) - deviation
        upper = trunc(sum_R/S) + deviation
        β = rand(lower:upper, S)
        return R, β
    end
end

function main()
    size = ARGS[1]
    write_jld2(size)
    @info "Done"
end

# main()
