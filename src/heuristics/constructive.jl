using Types
using JuMP
using Distances
using DelimitedFiles
using Gurobi
using Dates
using TimerOutputs
using MathOptInterface
const MOI = MathOptInterface
using Dates
using Random
function remove!(V, item)
    deleteat!(V, findall(x -> x == item, V))
end

function constructive(instance, id, init_method, assign_method; withdir=false, dir="")

    instancia = instance
    B = instancia.B
    S = instancia.S
    P = instancia.P
    X = Matrix{Int64}[]
    Y = Vector{Int64}[]
    D = instancia.D
    Weight = 0

    before1 = Dates.now()
    time_init1 = 0
    if init_method == "relax"
        Y, time_init1 = localize_facs(instancia, init_method)
        before1 = Dates.now()
    else
        Y = localize_facs(instancia, init_method)
    end

    if init_method ≠ "relax"
        Y_bool = zeros(Int, instancia.S)
        for idx in Y
            Y_bool[idx] = 1
        end
    else
        Y_bool = Y
    end
    after_init = Dates.now()
    delta_init1 = after_init - before1
    if init_method == "relax"
        delta_init1 = time_init1
    end
    #println("Y done with time: ", delta_init1)
    if assign_method == "naive"
        #println("NAIVE")
        X = naive_assign_bu(instancia, Y_bool)
    elseif assign_method == "opp"
        #println("OPP COST")
        X = oppCostAssignment(Y_bool, instancia)
    end
    after1 = Dates.now()
    delta_assign = after1 - after_init
    println("X done with time: ", delta_assign)
    delta1 = after1 - before1
    secs1 = round(delta1, Second)
    time1 = secs1.value
    if init_method == "relax"
        time1 += time_init1
    end

    before = Dates.now()
    time_init = 0
    if init_method == "relax"
        Y, time_init = localize_facs(instancia, init_method)
        before = Dates.now()
    else
        Y = localize_facs(instancia, init_method)
    end
    if init_method ≠ "relax"
        Y_bool = zeros(Int, instancia.S)
        for idx in Y
            Y_bool[idx] = 1
        end
    else
        Y_bool = Y
    end
    after_init2 = Dates.now()
    delta_init2 = after_init2 - before
    if init_method == "relax"
        delta_init2 = time_init
    end
    println("Y done with time: ", delta_init2)
    if assign_method == "naive"
        #println("NAIVE")
        X = naive_assign_bu(instancia, Y_bool)
    elseif assign_method == "opp"
        println("OPP COST")
        X = oppCostAssignment(Y_bool, instancia)
    end
    after = Dates.now()
    delta_assign2 = after - after_init2
    println("X done with time: ", delta_assign2)
    delta = after - before
    secs = round(delta, Second)
    time = secs.value
    if init_method == "relax"
        time += time_init
    end
    println("X done")
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += D[indice]
    end
    if time > time1
        time = time1
    end

    str_path = "sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * init_method * "_" * assign_method
    plot_str_path = str_path * ".png"
    solution_str_path = str_path * ".jld2"

    solution = Types.Solution(instancia, X, Y_bool, Weight, time)
    Types.plot_solution(solution, plot_str_path)
    Types.write_solution(solution, solution_str_path)
    # println(isFactible(solution, false))
    return solution
end

function localize_facs(instance, method)
    # k_type = findall(x->x==1, idx_candidate .∈ Sk) # get k type of fac
    if method == "pdisp"
        #println("P-DISP")
        return pdisp(instance)
    elseif method == "random"
        #println("RANDOM")
        return random_init(instance)
    elseif method == "relax"
        #println("RELAXATION")
        return relax_init(instance)
    end
end

function smarter_assign_bu(instance, Y)
    branches_used = findall(x -> x == 1, Y)
    B = instance.B
    S = instance.S
    X = zeros(Int, S, B)

    minimums = Tuple{Int,Tuple{Int,Int}}[]

    for j in 1:B
        minimum = 1e9
        i_exported = 0
        for i in 1:S
            if i ∉ branches_used
                continue
            end
            Dij = D[i, j]
            if Dij < minimum
                minimum = Dij
                i_exported = i
            end
        end
        min_and_index = (minimum, (i_exported, j))
        push!(minimums, min_and_index)
    end
    return X
end

function naive_assign_bu(instance, Y)
    # i haven't been using Y
    branches_used = findall(x -> x == 1, Y)
    # BUT, the entries in X must reflect the original D matrix
    B = instance.B
    S = instance.S
    K = instance.K
    D = instance.D
    X = zeros(Int, S, B)
    centers_used = 0
    for j in 1:B
        minimum = 1e9
        second_minimum = 1e9
        i_exported = 0
        for i in 1:S
            if i ∉ branches_used
                continue
            end
            Dij = D[i, j]
            if Dij < minimum
                second_minimum = minimum
                minimum = Dij
                i_exported = i
            end
        end
        X[i_exported, j] = 1
        centers_used += 1
    end
    return X
end

function minimums(matrix::Matrix, n)::Tuple{Vector{Int64},Vector{CartesianIndex{2}}}
    type = eltype(matrix)
    vals = fill(1e12, n)
    arr = Array{Tuple{type,CartesianIndex}}(undef, n)
    indices = Array{Int64}(undef, n)
    @inbounds for i ∈ axes(matrix, 1), j ∈ axes(matrix, 2)
        biggest, index = findmax(vals)
        if matrix[i, j] < biggest
            arr[index] = matrix[i, j], CartesianIndex(i, j)
            vals[index] = matrix[i, j]
        end
    end
    arr = sort(arr, by=x -> x[1])
    vals = [x[1] for x in arr]
    indices = [x[2] for x in arr]
    return vals, indices
end

function minimums(vec::Vector, n)::Tuple{Vector{Int64},Vector{Int64}}
    type = eltype(vec)
    vals = fill(1e9, n)
    arr = Array{Tuple{type,Int64}}(undef, n)
    indices = Array{Int64}(undef, n)
    @inbounds for i ∈ eachindex(vec)
        if vec[i] > 0
            biggest, index = findmax(vals)
            if vec[i] < biggest
                arr[index] = vec[i], i
                vals[index] = vec[i]
            end
        end
    end
    arr = sort(arr, by=x -> x[1])
    vals = [x[1] for x in arr]
    indices = [x[2] for x in arr]
    return vals, indices
end
#132407

function maximums(matrix, n)::Tuple{Vector{Int64},Vector{CartesianIndex{2}}}
    type = eltype(matrix)
    vals = zeros(type, n)
    indices = Array{Int64}(undef, n)
    arr = Array{Tuple{type,CartesianIndex}}(undef, n)
    @inbounds for i ∈ axes(matrix, 1), j ∈ axes(matrix, 2)
        smallest, index = findmin(vals)
        if matrix[i, j] > smallest
            arr[index] = matrix[i, j], CartesianIndex(i, j)
            vals[index] = matrix[i, j]
        end
    end
    arr = sort(arr, by=x -> x[1], rev=true)
    vals = [x[1] for x in arr]
    indices = [x[2] for x in arr]
    return vals, indices
end

function maximums2(M, n)
    v = vec(M)
    l = length(v)
    ix = [1:l;]
    partialsortperm!(ix, v, (l-n+1):l, initialized=true)
    vals = v[ix[(l-n+1):l]]
    indices = CartesianIndices(M)[ix[(l-n+1):l]]
    return vals, indices
end


function oppCostAssignment(Y, instance::Types.Instance)
    D = copy(instance.D)
    P = instance.P
    S = instance.S
    B = instance.B
    X = zeros(Int64, S, B)
    count = 0
    not_assigned_y = findall(y -> y == 0, Y)
    #println(not_assigned_y)
    for j in not_assigned_y
        D[j, :] .= -1
    end
    diff = copy(D)
    for i in 1:B
        minimal = 0
        minimals, _ = minimums(D[:, i], 1)
        minimal = minimals[1]
        diff[:, i] .= D[:, i] .- minimal
    end
    todos = false
    n = P - 5 # falta tweakear
    while !todos
        _, indices::Vector{CartesianIndex{2}} = maximums2(diff, n)
        constraints = Int64[]
        for indice::CartesianIndex in indices
            X_copy = Matrix{Int64}(undef, S, B)
            unsafe_copyto!(X_copy, 1, X, 1, S * B)
            col = indice[2]
            row = findfirst(x -> x == 0, diff[:, col])
            X_copy[row, col] = 1
            constraints_v = restricciones(X_copy, Y, instance; verbose=false)
            push!(constraints, constraints_v)
        end
        picked = CartesianIndex(1, 1)::CartesianIndex{2}
        if all(x -> x == 0, constraints) # si no se violan constraints, agarra el 0 de la columna del costo maximo
            indice = indices[1]::CartesianIndex{2}
            col = indice[2]::Int64
            row = findfirst(x -> x == 0, diff[:, col])
            picked = CartesianIndex(row, col)
        else
            if all(x -> x ≠ 0, constraints) # si todas violan constraints
                original_cons, idx = findmin(constraints) # agarra el que viole menos constraints
                indice = indices[idx]
                col = indice[2]
                original_row = findfirst(x -> x == 0, diff[:, col])::Int64
                # queremos buscar en esa columna el siguiente valor más cercano a 0 en diff
                # al hacerlo, nos acercamos al minimo valor de la matriz de distancias
                busqueda = trunc(Int, (P - 1)) # a lo mejor cambiarlo despues
                _, idxs_inner::Array{Int64} = minimums(diff[:, col], busqueda)
                # el primer minimo es el 0 de nuestra localizacion optima
                # lo desechamos para darle variedad a la busqueda
                idxs_inner2::Array{Int64} = idxs_inner[2:end]
                for row in idxs_inner2
                    X_copy = Matrix{Int64}(undef, S, B)
                    unsafe_copyto!(X_copy, 1, X, 1, S * B)
                    try
                        X_copy[row, col] = 1
                        constraints_v2 = restricciones(X_copy, Y, instance; verbose=false)
                        if constraints_v2 < original_cons
                            original_cons = constraints_v2
                            original_row = row
                        end
                    catch
                    end
                end
                picked = CartesianIndex(original_row, col)
            else # si hay una que no viola constraints
                for idx in eachindex(constraints)
                    if constraints[idx] == 0 # agarrala
                        indice = indices[idx]
                        col = indice[2]
                        row = findfirst(x -> x == 0, diff[:, col])
                        picked = CartesianIndex(row, col)
                        break
                    end
                end
            end
        end
        X[picked] = 1
        column = picked[2]::Int64
        diff[:, column] .= -1 # "apagamos" esa columna
        todos = true
        count += 1
        for col in eachcol(X)
            if all(x -> x == 0, col)
                todos = false
                break
            end
        end
    end
    return X
end

function restricciones(X_copy::Matrix{Int64}, Y_copy::Vector{Int64}, instance::Types.Instance; verbose=true)
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    Y = Y_copy
    X = X_copy
    Lk = instance.Lk
    Uk = instance.Uk
    Sk = instance.Sk
    μ = instance.μ
    T = instance.T
    M = instance.M
    S = instance.S
    P = instance.S
    B = instance.B
    β = instance.β
    R = instance.R
    K = instance.K
    V = instance.V
    number_constraints_violated = 0

    counts_k = Array{Int64}(undef, K)
    if sum(Y) > P
        if verbose
            println("Violando número de centros asignados ", sum(Y))
        end
    end
    for k_type in 1:K
        indices_k_type = Sk[k_type]
        count_k_type = 0
        for indice in indices_k_type
            if Y[indice] == 1
                count_k_type += 1
            end
        end
        counts_k[k_type] = count_k_type
    end

    for k in 1:K
        if counts_k[k] > Uk[k]
            if verbose
                println("Violando Uk en $k")
            end
            number_constraints_violated += 1
        end
    end
    for i in 1:S
        for m in 1:M
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= trunc(Int, (Y[i] * μ[m][i] * (1 + T[m]))))
                if verbose
                    println(i)
                    println(any(x -> x ≠ 0, X[i, :]))
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", trunc(Int, (Y[i] * μ[m][i] * (1 + T[m]))))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
        end
    end
    for i in 1:S
        if sum(X[i, j] * R[j] for j in 1:B) > β[i]
            if verbose
                println("violando riesgo en $i")
                println("β: ", β[i])
                println("R: ", sum(X[i, j] * R[j] for j in 1:B))
            end
            number_constraints_violated += 2
        end
    end
    return number_constraints_violated
end

function pdisp(instance)
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
    return S_sol
end

function random_init(instance)
    # tengo que agarrar los parametros para asignar una factible
    P = instance.P
    S = instance.S
    Y = Random.shuffle(collect(1:S))
    Y = Y[1:P]
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

    model = JuMP.Model(Gurobi.Optimizer) # THIS IS WHERE THE FUN BEGINS

    JuMP.@variable(model, x[1:S, 1:B], lower_bound = 0, upper_bound = 1)
    # num suc and num bu, Xᵢⱼ
    JuMP.@variable(model, y[1:S], Bin)
    # Yᵢ

    JuMP.@objective(model, Min, sum(D .* x))
    # Xᵢⱼ * Dᵢⱼ

    JuMP.@constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)

    # ∑ᵢ∈S Xᵢⱼ = 1, ∀ j ∈ B

    JuMP.@constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y[i])

    # Xᵢⱼ ≤ Yᵢ , ∀ i ∈ S, j ∈ B

    JuMP.@constraint(model, cardinality, sum(y) == P)

    # ∑ i ∈ S Yᵢ = p

    JuMP.@constraint(model, risk[j in 1:B, i in 1:S], x[i, j] * R[j] <= β[i])

    # ∑ j ∈ B Xᵢⱼ Rⱼ ≤ βᵢ, ∀ i ∈ S

    JuMP.@constraint(
        model,
        tol_l[i in 1:S, M in 1:m],
        y[i] * μ[M][i] * (1 - T[M]) <= sum(x[i, j] * V[m][j] for j in 1:B),
    )

    JuMP.@constraint(
        model,
        tol_u[i in 1:S, M in 1:m],
        sum(x[i, j] * V[m][j] for j in 1:B) <= y[i] * μ[M][i] * (1 + T[M]),
    )

    # Yᵢμₘⁱ(1-tᵐ) ≤ ∑i∈S Xᵢⱼ vⱼᵐ ≤ Yᵢμₘʲ(1+tᵐ) ∀ j ∈ B, m = 1 … 3

    JuMP.@constraint(
        model,
        low_k[K in 1:k],
        Lk[K] <= sum(y[i] for i in Sk[K]),
    )

    JuMP.@constraint(
        model,
        upp_k[K in 1:k],
        sum(y[i] for i in Sk[K]) <= Uk[K],
    )

    JuMP.set_silent(model)
    JuMP.optimize!(model)
    tiempo = MOI.get(model, MOI.SolveTimeSec())
    time_int = round(Int, tiempo)
    Y = round.(Int, JuMP.value.(model[:y]))
    return Y, time_int
end

function isFactible(solution::Types.Solution, verbose=true)
    number_constraints_violated = 0
    instance = solution.Instance
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    Y = solution.Y
    X = solution.X
    Lk = instance.Lk
    Uk = instance.Uk
    Sk = instance.Sk
    μ = instance.μ
    T = instance.T
    M = instance.M
    S = instance.S
    P = instance.S
    B = instance.B
    β = instance.β
    R = instance.R
    K = instance.K
    V = instance.V

    counts_k = []

    for y in eachindex(Y)
        if Y[y] == 1
            assignments_y = X[y, :]
            if !any(x -> x == 1, assignments_y)
                if verbose
                    println("Violando asignación de Y en: $y")
                end
                number_constraints_violated += 1
                # penalizar más aquí?
            end
        end
    end

    if sum(Y) > P
        if verbose
            println("Violando número de centros asignados ", sum(Y))
        end
        number_constraints_violated += 1
    end

    for j in 1:B
        if sum(X[i, j] for i in 1:S) ≠ 1
            if verbose
                println("Violando servicio a la BU: ", j)
                number_constraints_violated += 1
            end

        end
    end

    for k_type in 1:K
        indices_k_type = Sk[k_type]
        count_k_type = 0
        for indice in indices_k_type
            if Y[indice] == 1
                count_k_type += 1
            end
        end
        push!(counts_k, count_k_type)
    end

    for k in 1:K
        if counts_k[k] < Lk[k]
            if verbose
                println("Violando Lk en $k")
            end
            number_constraints_violated += 1
        end
        if counts_k[k] > Uk[k]
            if verbose
                println("Violando Uk en $k")
            end
            number_constraints_violated += 1
        end
    end

    for i in 1:S
        for m in 1:M
            if !(trunc(Int, Y[i] * μ[m][i] * (1 - 0.7)) <= sum(X[i, j] * V[m][j] for j in 1:B))
                if verbose
                    println("violando V inferior en i: $i y m: $m")
                    println("μ: ", trunc(Int, Y[i] * μ[m][i] * (1 - T[m])))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= trunc(Int, Y[i] * μ[m][i] * (1 + T[m])))
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", trunc(Int, Y[i] * μ[m][i] * (1 + T[m])))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
        end
    end

    for i in 1:S
        if sum(X[i, j] * R[j] for j in 1:B) > β[i]
            if verbose
                println("violando riesgo en $i")
                println("β: ", β[i])
                println("R: ", sum(X[i, j] * R[j] for j in 1:B))
            end
            number_constraints_violated += 1
        end
    end
    # println("Restricciones violadas: $number_constraints_violated")
    # intercambiar nodo i por nodo j 
    # mover nodo de territorio i a territorio j
    # apagar un branch y prender otro branch
    # checar ILS ??
    if number_constraints_violated ≠ 0
        return false, number_constraints_violated
    else
        return true, number_constraints_violated
    end
end

function main_constructive(init_method, assign_method; path="inst", read_file=true, instance_obj=nothing, id=0)
    instance = 1 # para traerlo al scope
    if read_file
        pattern = Regex("[t][_]\\d{1,3}")
        index = findfirst(pattern, path)
        almost_number = path[index]
        _, number = split(almost_number, "_")
        instance = Types.read_instance(path)
    else
        instance = instance_obj
        number = id
    end

    solution = constructive(instance, number, init_method, assign_method)
    write_solution(solution, "sol_300_60_20_1_viejomaximums.jld2")
    return solution
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_constructive(ARGS[2], ARGS[3]; path=ARGS[1])
end
