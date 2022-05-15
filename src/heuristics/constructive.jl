using Distances
using Types

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
    println(Y)
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
    Types.plot_solution(solution, "plot_solucion_3algo2_1.png")
    Types.write_solution(solution, "solucion_3algo2_1.jld2")
end

function localize_facs(instance; method=:pdisp)
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

    # no empezar con la columna 1, de manera lexicografica
    # puede sesgar 
    # hay una manera mas inteligente de seleccionar el orden de las BUs
    # escoger la BU que tenga la menor distancia O la mayor distancia
    # esto se hace para ir descartando los peores nodos para al final dejar los mejores
    # agarrar el minimo de los minimos
    # o agarrar el maximo de los minimos para asegurar que ya asigne el peor
    # otro criterio: el costo de oportunidad
    # diferencia de la segunda con la primera y eliges con eso
    # al principio las restricciones no se violan
    # al final, las restricciones se violan y la asigno a una facility muy lejana
    # poner un tope de decir: no asignar a la 5ta mas lejana aunque me viole mi restriccion
    # tope = P/3, P/4 o P/2

    # otra heuristica:
        # Yi entera
        # resuelveme una relajacion del problema
        # Xij, tope 0-1, no es entera, es continua
        # se va a resolver MUY rapido
        # PERO YA ME UBICO LOS CENTROS
        # en base a eso ya los agarro

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
    # meter al pseudocodigo
end

function pdisp(instance)
    P = instance.P
    s_coords = instance.S_coords
    D = instance.D
    S = instance.S
    metric = Euclidean()
    s_distances = trunc.(Int, pairwise(metric, s_coords, dims=1))
    idx_1, idx_2 = Tuple(argmax(s_distances)) # two furthest facs
    S_sol = [idx_1, idx_2] # S is our solution of our p-disp problem
    T = collect(1:S) # nodes not in our solution
    # se puede usar la suma de xⱼ a S_sol
    # o usar la facility mas cercana de xⱼ
    for s in S_sol
        remove!(T, s)
    end

    while length(S_sol) < P # extremely slow
        distances = NamedTuple[]
        idx_max_dist = 0
        for i in T
            distance = 0
            for j in S
                distance += s_distances[i,j]
            end
            tupla = (dist=distance, idx=i)
            push!(distances, tupla)
        end
        distances_max = [distance[1] for distance in distances]
        indexes = [distance[2] for distance in distances]
        idx_max = argmax(distances_max)
        idx_max_dist = indexes[idx_max]
        remove!(T, idx_max_dist)
        push!(S_sol, idx_max_dist)
    end
    return S_sol
end


function random_init(instance)
    # tengo que agarrar los parametros para asignar una factible
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