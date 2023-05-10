using JuMP, Gurobi, Types, Distances, Plots

function pdisp_model(instance::Instance)
    B = instance.B
    S = instance.S
    D = instance.D
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    P = instance.P
    V = instance.V
    μ = instance.μ
    T = instance.T
    R = instance.R
    β = instance.β
    k = 5
    s_coords = instance.S_coords
    model = Model(Gurobi.Optimizer)
    set_time_limit_sec(model, 50)
    metric = Distances.Euclidean()
    s_distances = trunc.(Int, Distances.pairwise(metric, s_coords, dims=1))
    @variable(model, y[1:S], Bin)
    @constraint(model, cardinality, sum(y) == P)
    @constraint(
        model,
        low_k[K in 1:k],
        Lk[K] <= sum(y[i] for i in Sk[K]),
    )

    @constraint(
        model,
        upp_k[K in 1:k],
        sum(y[i] for i in Sk[K]) <= Uk[K],
    )

    # lₖ ≤ ∑i ∈ Sₖ Yᵢ ≤ uₖ, k = 1 … 5
    @objective(model, Max, sum(s_distances[i, j] * y[i] * y[j] for i in 1:S, j in i+1:S))
    

    optimize!(model)
    tiempo = MOI.get(model, MOI.SolveTimeSec())
    time_int = trunc(Int, tiempo)

    Y = value.(model[:y])
    Y = round.(Int, Y)
    return Y
end

function remove!(V, item)
    deleteat!(V, findall(x -> x == item, V))
end

function pdisp_simple(d, p, N)
    maxdist  = 0
    bestpair = (0, 1)
    for i in 1:N
        for j in i+1:N
            if d[i,j]>maxdist
                maxdist = d[i,j]
                bestpair = (i,j)
            end
        end
    end
    P = Set([])
    push!(P, bestpair[1])
    push!(P, bestpair[2])

    while length(P) < p
        maxdist = 0
        vbest = 0
        for v in 1:N
            if v in P
                continue
            end
            mindist = Inf
            for vprime in P
                if d[v,vprime] < mindist
                    mindist = d[v,vprime]
                end
            end
            if mindist > maxdist
                maxdist = mindist
                vbest = v
            end
        end
        if vbest != 0 && !(vbest in P)
            
            push!(P, vbest)
        end
    end
    collection = collect(P)
    return collection
end

function count_k(P, Sk)
    count = zeros(Int, length(Sk))
    for i in P
        k = node_type(i, Sk)
        count[k] += 1
    end
    return count
end

function node_type(i, Sk)
    for k in eachindex(Sk)
        if i in Sk[k]
            return k
        end
    end
    println("Node $i not found in Sk")
end


function pdisp_3(instance)
    B = instance.B
    S = instance.S
    D = instance.D
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    P = instance.P
    V = instance.V
    μ = instance.μ
    T = instance.T
    R = instance.R
    β = instance.β
    k = 5
    s_coords = instance.S_coords
    metric = Distances.Euclidean()
    d = trunc.(Int, Distances.pairwise(metric, s_coords, dims=1))
    coords_S1 = s_coords[Sk[1],:]
    coords_S2 = s_coords[Sk[2],:]
    coords_S3 = s_coords[Sk[3],:]
    coords_S4 = s_coords[Sk[4],:]
    coords_S5 = s_coords[Sk[5],:]
    
    d1 = trunc.(Int, Distances.pairwise(metric, coords_S1, dims=1))
    d2 = trunc.(Int, Distances.pairwise(metric, coords_S2, dims=1))
    d3 = trunc.(Int, Distances.pairwise(metric, coords_S3, dims=1))
    d4 = trunc.(Int, Distances.pairwise(metric, coords_S4, dims=1))
    d5 = trunc.(Int, Distances.pairwise(metric, coords_S5, dims=1))
    N1 = length(Sk[1]) 
    N2 = length(Sk[2])
    N3 = length(Sk[3])
    N4 = length(Sk[4])
    N5 = length(Sk[5])
    p1 = Lk[1]
    p2 = Lk[2]
    p3 = Lk[3]
    p4 = Lk[4]
    p5 = Lk[5]

    pdisp1 = pdisp_simple(d1, p1, N1)
    pdisp2 = pdisp_simple(d2, p2, N2)
    pdisp3 = pdisp_simple(d3, p3, N3)
    pdisp4 = pdisp_simple(d4, p4, N4)
    pdisp5 = pdisp_simple(d5, p5, N5)

    pdisp1_fixed = Sk[1][pdisp1]
    pdisp2_fixed = Sk[2][pdisp2]
    pdisp3_fixed = Sk[3][pdisp3]
    pdisp4_fixed = Sk[4][pdisp4]
    pdisp5_fixed = Sk[5][pdisp5]

    pdisp_ok = Set(vcat([pdisp1_fixed, pdisp2_fixed, pdisp3_fixed, pdisp4_fixed, pdisp5_fixed]...))
    if length(pdisp_ok) != P
        count = count_k(pdisp_ok, Sk)
        while length(pdisp_ok) < P
            # Find the node v that maximizes the distance to its closest neighbor in P
            maxdist = 0
            vbest = 0
            for v in 1:N
                if v in pdisp_ok
                    continue
                end
                k = node_type(v, Sk)
                if count[k] >= Uk[k]
                    continue
                end
                dist = minimum([d[v,vprime] for vprime in pdisp_ok])
                if dist > maxdist
                    maxdist = dist
                    vbest = v
                end
            end
            # If no such node exists, stop the algorithm
            if vbest == 0
                break
            end
            # Add the node vbest to the set P and update the counts
            k = node_type(vbest, Sk)
            count[k] += 1
            push!(pdisp_ok, vbest)
        end
    end
    collection = collect(pdisp_ok)
    Y_bool = zeros(Int, instance.S)
    for idx in collection
        Y_bool[idx] = 1
    end
    return Y_bool
end
#=
function pdisp(d, p, N, Sk, Lk, Uk)
    # Find the pair of nodes with maximum distance
    maxdist = 0
    bestpair = (0, 1)
    for i in 1:N
        for j in i+1:N
            if d[i,j] > maxdist
                maxdist = d[i,j]
                bestpair = (i,j)
            end
        end
    end
    P = Set([bestpair[1], bestpair[2]])
    while length(P) < p
        # Find the node v that maximizes the distance to its closest neighbor in P
        maxdist = 0
        vbest = 0
        for v in 1:N
            if v in P
                continue
            end
            k = node_type(v, Sk)
            if count[k] >= Uk[k]
                continue
            end
            dist = minimum([d[v,vprime] for vprime in P])
            if dist > maxdist
                maxdist = dist
                vbest = v
            end
        end
        # If no such node exists, stop the algorithm
        if vbest == 0
            println("ni peyo")
            break
        end
        # Add the node vbest to the set P and update the counts
        k = node_type(vbest, Sk)
        count[k] += 1
        push!(P, vbest)
    end
    # Check if the solution satisfies the lower bounds
    count = count_k(P, Sk)
    for k = 1:length(Sk)
        if count[k] < Lk[k]
            println("puta mamada")
            return Set([])
        end
    end
    return P
end
=#

function pdisp_2(instance)
    B = instance.B
    S = instance.S
    D = instance.D
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    p = instance.P
    V = instance.V
    μ = instance.μ
    T = instance.T
    R = instance.R
    β = instance.β
    k = 5
    s_coords = instance.S_coords
    metric = Distances.Euclidean()
    d = trunc.(Int, Distances.pairwise(metric, s_coords, dims=1))
    N = S
    maxdist  = 0
    bestpair = (0, 1)
    for i in 1:N
        for j in i+1:N
            if d[i,j]>maxdist
                maxdist = d[i,j]
                bestpair = (i,j)
            end
        end
    end

    P = Set([])
    push!(P, bestpair[1])
    push!(P, bestpair[2])

    while length(P) < p
        maxdist = 0
        vbest = 0
        for v in 1:N
            if v in P
                continue
            end
            for vprime in P
                if d[v,vprime]>maxdist
                    maxdist = d[v,vprime]
                    # indexar el arreglo de las Sk para saber que tipo es v
                    # if Lk  >  countk
                    vbest   = v
                end
            end
        push!(P, vbest)
        end
    end
    collection = collect(P)
    Y_bool = zeros(Int, instance.S)
    for idx in collection
        Y_bool[idx] = 1
    end
    return Y_bool
end

function pdisp_heuristic(instance)
    P = instance.P
    s_coords = instance.S_coords
    S = instance.S
    metric = Distances.Euclidean()
    s_distances = trunc.(Int, Distances.pairwise(metric, s_coords, dims=1))
    idx_1, idx_2 = Tuple(argmax(s_distances)) # two furthest facs
    S_sol = [idx_1, idx_2] # S is our solution of our p-disp problem
    T = collect(1:S) # nodes not in our solution
    # se puede usar la suma de xⱼ a S_sol
    # o usar la facility mas cercana de xⱼ
    for s in S_sol
        remove!(T, s)
    end

    while length(S_sol) < P
        distances = Tuple{Int64,Int64}[]
        idx_max_dist = 0
        for i in T
            distance = 0
            for j in S
                distance += s_distances[i, j]
            end
            tupla = (distance, i)
            push!(distances, tupla)
        end
        distances_max = [distance[1] for distance in distances]::Array{Int64}
        indexes = [distance[2] for distance in distances]::Array{Int64}
        idx_max = argmax(distances_max)::Int64
        idx_max_dist = indexes[idx_max]
        remove!(T, idx_max_dist)
        push!(S_sol, idx_max_dist)
    end
    Y_bool = zeros(Int, instance.S)
    for idx in S_sol
        Y_bool[idx] = 1
    end
    return Y_bool
end

function main()
    instance = read_instance("instances\\800_150_50\\inst_1_800_150_50.jld2")
    println(@time pdisp_heuristic(instance))
    # println(@time pdisp_model(instance))
    println(@time pdisp_2(instance))
    println(@time pdisp_3(instance))
    Y_heur_1 = pdisp_heuristic(instance)
    Y_heur_2 = pdisp_2(instance)
    Y_heur_3 = pdisp_3(instance)
    Y_model =  pdisp_model(instance)
    Y_heur1_bool = Bool.(Y_heur_1)
    Y_heur2_bool = Bool.(Y_heur_2)
    Y_heur3_bool = Bool.(Y_heur_3)
    Y_model_bool = Bool.(Y_model)
    S = instance.S
    plot_font = "Computer Modern"
    default(fontfamily=plot_font,framestyle=:box, grid=false, tickfontsize=7)
    S_coords = instance.S_coords
    
    p = Plots.scatter(
        S_coords[:,1],
        S_coords[:,2],
        markershape = :square,
        markercolor = [(Y_heur1_bool[i] ? :red : :white ) for i in 1:S],
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        label = nothing,
    )

    Plots.pdf(p, "plot_heur1.pdf")

    p = Plots.scatter(
        S_coords[:,1],
        S_coords[:,2],
        markershape = :square,
        markercolor = [(Y_heur2_bool[i] ? :red : :white ) for i in 1:S],
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        label = nothing,
    )

    Plots.pdf(p, "plot_heur2.pdf")

    p = Plots.scatter(
        S_coords[:,1],
        S_coords[:,2],
        markershape = :square,
        markercolor = [(Y_heur3_bool[i] ? :red : :white ) for i in 1:S],
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        label = nothing,
    )

    Plots.pdf(p, "plot_heur3.pdf")

    p2 = Plots.scatter(
        S_coords[:,1],
        S_coords[:,2],
        markershape = :square,
        markercolor = [(Y_model_bool[i] ? :red : :white) for i in 1:S],
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        label = nothing,
    )

    Plots.pdf(p2, "plot_model.pdf")
end

main()