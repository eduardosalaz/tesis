using Distances, JLD2, Random
using Types

function generate_instance(size::String, i::Int; write=true)
    K = 4
    M = 3
    B = 0
    S = 0
    P = 0
    if size == "S"
        B = 625
        S = 200
        P = 30
        size = "625_300_30"
    elseif size == "M"
        B = 1250
        S = 200
        P = 40
        size = "1250_155_62"
    elseif size == "L"
        B = 2500
        S = 400
        P = 40
        size = "2500_325_125"
    elseif size == "XS"
        B = 312
        S = 40
        P = 16
    end
    BU_coords, S_coords = generate_coords(B, S)
    dist_mat = generate_dist(BU_coords, S_coords, B, S)
    parameters = generate_params(B, S, P)
    instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat, parameters...)
    dir_path = "instances_new/new_ranges_005_newk" * size * "/"
    if !isdir(dir_path)
        mkdir(dir_path)
    end
    file_name = "inst_" * string(i) * "_" * size * ".jld2"
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
    return round.(Int, mat)
end

function generate_params(B::Int64, S::Int64, P)
    Sk, Lk, Uk, percentages = generate_k(S, P)
    @debug "Wrote k"
    V, μ, T = generate_activities(B, S, P)
    @debug "Wrote activities"
    R, β = generate_risk(B, S, P)
    @debug "Wrote risk"
    return Sk, Lk, Uk, V, μ, T, R, β
end

function generate_k(S::Int64, p::Int64)  # Now we need p as input
    K = 4
    percentages = [0.4, 0.3, 0.2, 0.1]
    shuffle!(percentages)  # Randomly assign percentages to types
    
    # Calculate how many facilities should be of each type
    counts_per_type = round.(Int, S .* percentages)
    
    # Adjust for rounding errors to ensure sum equals S
    while sum(counts_per_type) != S
        if sum(counts_per_type) < S
            idx = argmax(percentages .- counts_per_type/S)
            counts_per_type[idx] += 1
        else
            idx = argmin(percentages .- counts_per_type/S)
            counts_per_type[idx] -= 1
        end
    end
    
    # Create the assignment vector
    Ks = Int64[]
    for k in 1:K
        append!(Ks, fill(k, counts_per_type[k]))
    end
    shuffle!(Ks)  # Randomize the order
    
    # Create Sk (facilities of each type)
    Sk = [findall(x -> x == k, Ks) for k in 1:K]
    
    # Calculate bounds based on p, not on counts_per_type
    # For each type k, we want between (percentage-0.05)*p and (percentage+0.05)*p centers
    Lk = round.(Int, (percentages .- 0.05) .* p)
    Uk = round.(Int, (percentages .+ 0.05) .* p)
    
    # Ensure Lk isn't negative and Uk doesn't exceed available facilities
    Lk = max.(0, Lk)
    Uk = min.(Uk, counts_per_type)
    
    return Sk, Lk, Uk, percentages
end

function generate_activities(B::Int64, S::Int64, P)
    M = 3
    # Initialize V with M arrays, each of length B
    V = [zeros(Int64, B) for _ in 1:M]
    
    # Fill each array in V with random numbers according to their ranges
    for i in 1:B
        V[1][i] = rand(1:10)
        V[2][i] = rand(1000:10000)
        V[3][i] = rand(1000:5000)
    end
    
    # Initialize μ with M arrays, each of length S
    μ = [zeros(Int64, S) for _ in 1:M]
    
    # Calculate μ values based on sum of corresponding V array
    for i in 1:M
        sum_vals = sum(V[i])
        base_value = round(Int, sum_vals / P)
        for s in 1:S
            μ[i][s] = base_value
        end
    end
    
    T = [0.05, 0.05, 0.05]
    return V, μ, T
end

function generate_risk(B::Int64, S::Int64, P)
    R = rand(30:60, B)
    sum_R = sum(R)
    # rmin + rmax  / 2 (multiplicado por unidades basicas) dividido entre P multiplicado por un valor adicional
    τ = 0.1
    β = zeros(Int64, S)
    for s in 1:S
        base_value = round(Int, (((30 + 60 / 2) * B) / P)  * (1 + τ))
        β[s] = round(Int, base_value)
    end
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

# main()
