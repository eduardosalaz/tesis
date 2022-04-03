include("../types/types.jl")
using Distances, JLD2, Plots, Random, .Types

function write_jld2(size::String)
    K = 5
    M = 3
    if size == "S"
        B = 25
        S = 15
        P = 10
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
    @info "Wrote plot and coords"
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
        V, μ, T = generate_activities(size, B, S)
        R, β = generate_risk(size, B, S)

        return Sk, Lk, Uk, V, μ, T, R, β
    elseif size == "M"
        Sk, Lk, Uk = generate_k(size, B, S)
        V, μ, T = generate_activities(size, B, S)
        R, β = generate_risk(size, B, S)
        return Sk, Lk, Uk, V, μ, T, R, β
    elseif size == "L"
        Sk, Lk, Uk = generate_k(size, B, S)
        V, μ, T = generate_activities(size, B, S)
        R, β = generate_risk(size, B, S)
        return Sk, Lk, Uk, V, μ, T, R, β
    elseif size == "XL"
        Sk, Lk, Uk = generate_k(size, B, S)
        V, μ, T = generate_activities(size, B, S)
        R, β = generate_risk(size, B, S)
        return Sk, Lk, Uk, V, μ, T, R, β
    end
end

chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)]

function generate_k(size::String, B::Int64, S::Int64)
    if size == "S"
        K = 5
        size_sk = S / K |> Int
        shuffled_s = shuffle(collect(1:S))
        Sk = chunk(shuffled_s, size_sk)
        tope_Uk = length(Sk[1]) - 1
        Uk = rand(1:tope_Uk, K)
        Lk = rand(0:1, K)
        @info "Wrote k"
        return Sk, Lk, Uk
    elseif size == "M"
        P = 25
        K = 5
        size_sk = S / K |> Int
        shuffled_s = shuffle(collect(1:S))
        Sk = chunk(shuffled_s, size_sk)
        tope_Uk = length(Sk[1]) - 1
        Uk = rand(6:tope_Uk, K)
        Lk = rand(2:4, K)
        sum_uk = true
        while sum_uk
            if sum(Uk) < S
                Uk = rand(6:tope_Uk, K)
            else
                break
            end
        end
        sum_lk = true
        while sum_lk
            if sum(Lk) > S
                Lk = rand(2:4, K)
            else
                break
            end
        end
        @info "Wrote k"
        return Sk, Lk, Uk
    elseif size == "L"
        P = 50
        K = 5
        size_sk = S / K |> Int
        shuffled_s = shuffle(collect(1:S))
        Sk = chunk(shuffled_s, size_sk)
        tope_Uk = length(Sk[1]) - 1
        Uk = rand(10:tope_Uk, K)
        Lk = rand(4:6, K)
        sum_uk = true
        while sum_uk
            if sum(Uk) < S
                Uk = rand(10:tope_Uk, K)
            else
                break
            end
        end
        sum_lk = true
        while sum_lk
            if sum(Lk) > S
                Lk = rand(4:6, K)
            else
                break
            end
        end
        @info "Wrote k"
        return Sk, Lk, Uk
    elseif size == "XL"
        P = 70
        K = 5
        size_sk = S / K |> Int
        shuffled_s = shuffle(collect(1:S))
        Sk = chunk(shuffled_s, size_sk)
        tope_Uk = length(Sk[1]) - 1
        Uk = rand(15:tope_Uk, K)
        Lk = rand(7:9, K)
        sum_uk = true
        while sum_uk
            if sum(Uk) < S
                Uk = rand(15:tope_Uk, K)
            else
                break
            end
        end
        sum_lk = true
        while sum_lk
            if sum(Lk) > S
                Lk = rand(7:9, K)
            else
                break
            end
        end
        @info "Wrote k"
        return Sk, Lk, Uk
    end
end

function generate_activities(size::String, B::Int64, S::Int64)
    M = 3
    if size == "S"
        V = [zeros(Int64, B) for _ in 1:M]
        for i in 1:length(V)
            V[i] = rand(10:20, B)
        end
        μ = [zeros(Int64, S) for _ in 1:M]
        for i in 1:length(μ)
            sum_vals = sum(V[i])
            deviation = trunc(Int, sum_vals/100)
            lower = trunc(sum_vals/S) - deviation
            upper = trunc(sum_vals/S) + deviation
            μ[i] = rand(lower:upper, S)
        end
        T = fill(0.4, M)
        @info "Wrote activities"
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
        @info "Wrote activities"
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
        @info "Wrote activities"
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
        @info "Wrote activities"
        return V, μ, T
    end
end

function generate_risk(size::String, B::Int64, S::Int64)
    if size == "S"
        R = rand(5:10, B)
        sum_R = sum(R)
        deviation = trunc(Int, sum_R/80)
        lower = trunc(sum_R/S) - deviation
        upper = trunc(sum_R/S) + deviation
        β = rand(lower:upper, S)
        @info "Wrote risk"
        return R, β
    elseif size == "M"
        R = rand(10:20, B)
        sum_R = sum(R)
        deviation = trunc(Int, (sum_R/80))
        lower = trunc(sum_R/S) - deviation
        upper = trunc(sum_R/S) + deviation
        β = rand(lower:upper, S)
        @info "Wrote risk"
        return R, β
    elseif size == "L"
        R = rand(15:25, B)
        sum_R = sum(R)
        deviation = trunc(Int, sum_R/80)
        lower = trunc(sum_R/S) - deviation
        upper = trunc(sum_R/S) + deviation
        β = rand(lower:upper, S)
        @info "Wrote risk"
        return R, β
    elseif size == "XL"
        R = rand(20:35, B)
        sum_R = sum(R)
        deviation = trunc(Int, sum_R/80)
        lower = trunc(sum_R/S) - deviation
        upper = trunc(sum_R/S) + deviation
        β = rand(lower:upper, S)
        @info "Wrote risk"
        return R, β
    end
end

function write_file(mat, B, S, num_BU, num_Suc)
    open("instance_numBu$num_BU" * "_num_Suc" * "$num_Suc" * ".txt", "w") do io
        write(io, "S\n")
        write(io, "1 3\n") # s1
        write(io, "2 5\n") # s2
        write(io, "4\n") # s3
        write(io, "6 7\n") #s4
        write(io, "8\n") # s5

        write(io, "U\n")
        write(io, "2\n") #u1
        write(io, "2\n")
        write(io, "1\n")
        write(io, "2\n")
        write(io, "1\n")

        write(io, "L\n")
        write(io, "1\n") # l1  rango de 1 ≤ 2 = 1 o 2
        write(io, "1\n") # l2 rango de 1 ≤ 2 = 1 o 2
        write(io, "1\n") # l3 rango de 1 ≤ 1 = 1
        write(io, "1\n") # l4 rango de 1 ≤ 2 = 1 o 2
        write(io, "0\n") # l5 rango de 0 ≤ 1 = 0 o 1

        write(io, "P\n")
        write(io, "8\n") # solo 8 centros

    end
end

function main()
    size = ARGS[1]
    write_jld2(size)
    @info "Done"
end

main()
