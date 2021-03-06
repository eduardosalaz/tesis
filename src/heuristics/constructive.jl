using Distances
using Types
using JuMP, Gurobi

function remove!(V, item)
    deleteat!(V, findall(x -> x == item, V))
end

function constructive(instance)
    instancia = read_instance(instance)
    X = Matrix{Int64}[]
    Y = Vector{Int64}[]
    D = instancia.D
    method = :pdisp
    Weight = 0
    Y = localize_facs(instancia, method)
    println(Y)
    X = naive_assign_bu(instancia, Y)
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += D[indice]
    end
    if method ≠ :relax
        Y_bool = zeros(instancia.S)
        for idx in Y
            Y_bool[idx] = 1
        end
    else
        Y_bool = Y
    end
    solution = Solution(instancia, X, Y_bool, Weight)
    plot_solution(solution, "pruebas/plot_solucion_1algo31_pdisp.png")
    write_solution(solution, "pruebas/solucion_1algo31_pdisp.jld2")
end

function localize_facs(instance,method)
    # k_type = findall(x->x==1, idx_candidate .∈ Sk) # get k type of fac
    if method == :pdisp
        println("P-DISP")
        return pdisp(instance)
    elseif method == :random
        println("RANDOM")
        return random_init(instance)
    elseif method == :relax
        println("RELAXATION")
        return relax_init(instance)
    end

end

function minimize_assign(instance, Y)
    B = instance.B
    S = instance.S
    X = zeros(Int, S, B)
    D = instance.D

    branches_not_used = findall(x->x==0, Y)
    diff = [branches_not_used[i] - branches_not_used[i-1] for i in 2:length(branches_not_used)]
    signos = sign.(diff)
    D = instance.D
    D_copy = copy(D)
    pushfirst!(signos, 0) # first offset is 0 obviously as there is no previous row to be removed
    for i in eachindex(signos)
        D = D[1:end .≠ branches_not_used[i] - signos[i], :] # flawless
    end

    for j in 1:B
        min_row, min_col = argmin(D) # get the global minimum of the distance matrix
        min_val = D[min_row, min_col]
        original_row, original_col = findall(x->x==min_val, D_copy)
        # BUT, the entries in X must reflect the original D matrix
        # this breaks if there aren't unique entries in D
        X[original_row, original_col] = 1
        D = D[:, 1:end .≠ min_col] # remove from D the column which has that minimal value
        # in theory every column should be served?
    end
    return X
end

function naive_assign_bu(instance, Y)
    # i haven't been using Y
    branches_used = findall(x->x==1, Y)
    # BUT, the entries in X must reflect the original D matrix
    B = instance.B
    S = instance.S
    X = zeros(Int, S, B)

    minimums = Tuple{Int, Tuple{Int, Int}}[]

    for j in 1:B
        minimum = 1e9
        i_exported = 0
        for i in 1:S
            if i ∉ branches_used
                continue
            end
            Dij = D[i,j]
            if Dij < minimum
                minimum = Dij
                i_exported = i
            end
        end
        X[i_exported,j] = 1
        # min_and_index = (minimum, (i_exported, j))
        # push!(minimums, min_and_index)
    end
    
    Sk = instance.Sk



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
                distance += s_distances[i, j]
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


function relax_init(instance)

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
    m = instance.M
    k = instance.K

    model = Model(Gurobi.Optimizer) # THIS IS WHERE THE FUN BEGINS

    @variable(model, x[1:S, 1:B], lower_bound=0, upper_bound=1)
    # num suc and num bu, Xᵢⱼ
    @variable(model, y[1:S], Bin)
    # Yᵢ

    @objective(model, Min, sum(D .* x))
    # Xᵢⱼ * Dᵢⱼ

    @constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)

    # ∑ᵢ∈S Xᵢⱼ = 1, ∀ j ∈ B

    @constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y[i])

    # Xᵢⱼ ≤ Yᵢ , ∀ i ∈ S, j ∈ B

    @constraint(model, cardinality, sum(y) == P)

    # ∑ i ∈ S Yᵢ = p

    @constraint(model, risk[j in 1:B, i in 1:S], x[i, j] * R[j] <= β[i])

    # ∑ j ∈ B Xᵢⱼ Rⱼ ≤ βᵢ, ∀ i ∈ S

    @constraint(
        model,
        tol_l[i in 1:S, M in 1:m],
        y[i] * μ[M][i] * (1 - T[M]) <= sum(x[i, j] * V[m][j] for j in 1:B),
    )

    @constraint(
        model,
        tol_u[i in 1:S, M in 1:m],
        sum(x[i, j] * V[m][j] for j in 1:B) <= y[i] * μ[M][i] * (1 + T[M]),
    )

    # Yᵢμₘⁱ(1-tᵐ) ≤ ∑i∈S Xᵢⱼ vⱼᵐ ≤ Yᵢμₘʲ(1+tᵐ) ∀ j ∈ B, m = 1 … 3

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

    set_silent(model)
    optimize!(model)
    Y = trunc.(Int, value.(model[:y]))
    println(value.(model[:x]))
    return Y
end


constructive(ARGS[1])