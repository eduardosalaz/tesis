using Distances
include("../types/types.jl")
using .Types

function remove!(V, item)
    deleteat!(V, findall(x->x==item, V))
end

function constructive(instance)
    instancia = Types.read_instance(instance)
    X = Matrix{Int64}[]
    Y = Vector{Int64}[]
    D = instancia.D
    Weight = 0
    Y = localize_facs(instancia)
    X = naive_assign_bu(instancia, Y)
    indices = findall(x->x==1, X)
    for indice in indices
        Weight += D[indice]
    end
    Y_bool = zeros(instancia.S)
    for idx in Y
        Y_bool[idx] = 1
    end
    solution = Types.Solution(instancia, X, Y_bool, Weight)
    Types.plot_solution(solution, "plot_solucion_algo.png")
    Types.write_solution(solution, "solucion_algo.jld2")
end

function localize_facs(instance; method=:random)
    # k_type = findall(x->x==1, idx_candidate .∈ Sk) # get k type of fac
    if method == :pdisp
        return pdisp(instance)
    elseif method == :random
        return random_init(instance)
    end

end

function naive_assign_bu(instance, Y)
    B = instance.B
    S = instance.S
    X = zeros(Int, S, B)
    D = instance.D
    Sk = instance.Sk

    for j in 1:B
        Dⱼ = D[:, j] # get each column as a vector
        sort!(Dⱼ) # sort to get min distances in order
        candidate = Dⱼ[1] # first should be the minimal distance of the BU j to any fac
        idx_candidate = findfirst((x->x==candidate), D[:, j]) # we dont search in Dⱼ as it has mixed idxs
        X[idx_candidate, j] = 1
    end
    return X


    # una fila representa una facility
    # una columna representa una BU
    # quiero encontrar el minimo de cada columna (la distancia minima de cada BU)
    # minimo de columna 1, minimo de distancia de BU1 a cualquier facility
    # irme columna por oclumna
    # agarrar el minimo
    # calcular si se respetan las tolerancias
    # si si, marcar en X 
    # si no, agarrar el segundo minimo hasta que haya uno
    # tengo que quitar las filas que no escogi en Y primero D[Y, :]
    # each column represents a facility


end

function pdisp(instance)
    P = instance.P
    s_coords = instance.S_coords
    D = instance.D
    S = instance.S
    metric = Euclidean()
    s_distances = pairwise(metric, s_coords, dims=1)
    idx_1, idx_2 = Tuple(argmax(s_distances)) # two furthest facs
    S_sol = Set([idx_1, idx_2]) # S is our solution of our p-disp problem
    T = Set(collect(1:S)) # nodes not in our solution
    for s in S_sol
        delete!(T, s)
    end

    while length(S_sol) < P # extremely slow
        min_distances = NamedTuple[]
        for node_sol in S_sol
            dist = []
            for node_candidate in T
                push!(dist, s_distances[node_sol, node_candidate])
            end
            min_dist = minimum(dist)
            min_idx = argmin(dist)
            push!(min_distances, (idx = min_idx, distance = min_dist))
        end
        distances = [min_dist.distance for min_dist in min_distances]
        indexes = [min_dist.idx for min_dist in min_distances]
        idx_max = argmax(distances)
        idx_max_dist = indexes[idx_max] # naming could be better
        push!(S_sol, idx_max_dist)
        delete!(T, idx_max_dist)
    end
    
end


function random_init(instance)
    P = instance.P
    S_coords = instance.S_coords
    S = instance.S
    S_vec = collect(1:S)
    Y = Int64[]
    for p in 1:P
        rand_elem = S_vec[rand(1:end)]
        push!(Y, rand_elem)
        remove!(S_vec, rand_elem)
    end
    return Y
end


constructive(ARGS[1])