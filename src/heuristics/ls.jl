using Types, BenchmarkTools, DelimitedFiles
using TimerOutputs
using LinearAlgebra
using Dates
using StaticArrays

function find_one_in_column_unrolled(X::Matrix{Int64}, col::Int)
    rows = size(X, 1)
    @inbounds begin
        for i in 1:4:rows-3
            X[i, col] == 1 && return i
            X[i+1, col] == 1 && return i + 1
            X[i+2, col] == 1 && return i + 2
            X[i+3, col] == 1 && return i + 3
        end
        for i in (rows&-4+1):rows
            X[i, col] == 1 && return i
        end
    end
    return 0
end

# Highly optimized function to find all 1s in a row
function find_ones_in_row_optimized(X::Matrix{Int64}, row::Int)
    cols = size(X, 2)
    indices = Vector{Int}(undef, cols)  # Pre-allocate for worst case
    count = 0

    # Process 4 elements at a time
    @inbounds for j in 1:4:cols-3
        chunk = (X[row, j] == 1) | (X[row, j+1] == 1) << 1 | (X[row, j+2] == 1) << 2 | (X[row, j+3] == 1) << 3
        count += count_ones(chunk)

        chunk == 0 && continue

        chunk & 1 != 0 && (indices[count-count_ones(chunk)+1] = j)
        chunk & 2 != 0 && (indices[count-count_ones(chunk & 0b1110)+1] = j + 1)
        chunk & 4 != 0 && (indices[count-count_ones(chunk & 0b1100)+1] = j + 2)
        chunk & 8 != 0 && (indices[count] = j + 3)
    end

    # Handle remaining elements
    @inbounds for j in (cols&-4+1):cols
        if X[row, j] == 1
            count += 1
            indices[count] = j
        end
    end

    return resize!(indices, count)
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
    P = instance.P
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

    if sum(Y) != P
        if verbose
            println("Violando número de centros asignados ", sum(Y), " contra ", P)
        end
        number_constraints_violated += 1
    end

    for j in 1:B
        if sum(X[i, j] for i in 1:S) ≠ 1
            if verbose
                println("Violando servicio a la BU: ", j)
                println(sum(X[i, j] for i in 1:S))
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
            if !(round(Int, Y[i] * μ[m][i] * (1 - T[m])) <= sum(X[i, j] * V[m][j] for j in 1:B))
                if verbose
                    println("violando V inferior en i: $i y m: $m")
                    println("μ: ", round(Int, (Y[i] * μ[m][i] * (1 - T[m]))))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= round(Int, (Y[i] * μ[m][i] * (1 + T[m]))))
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", round(Int, (Y[i] * μ[m][i] * (1 + T[m]))))
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

function calculate_targets(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    S = instance.S
    targets_lower = Matrix{Int}(undef, S, M)
    targets_upper = Matrix{Int}(undef, S, M)
    for s in 1:S
        for m in 1:M
            targets_lower[s, m] = round(Int, (1 * μ[m][s] * (1 - T[m])))
            targets_upper[s, m] = round(Int, (1 * μ[m][s] * (1 + T[m])))
        end
    end
    return targets_lower, targets_upper
end

function calculate_targets_upper(instance)
    M = instance.M
    S = instance.S
    μ = instance.μ
    T = instance.T
    # Now we need targets for each center and each activity type
    targets = Matrix{Int}(undef, S, M)
    for s in 1:S
        for m in 1:M
            targets[s, m] = round(Int, (μ[m][s] * (1 + T[m])))
        end
    end
    return targets
end

function calculate_targets_lower(instance)
    M = instance.M
    S = instance.S
    μ = instance.μ
    T = instance.T
    targets = Matrix{Int}(undef, S, M)
    for s in 1:S
        for m in 1:M
            targets[s, m] = round(Int, (μ[m][s] * (1 - T[m])))
        end
    end
    return targets
end

function calculate_targets_optimized(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    S = instance.S

    # Create matrices to store lower and upper bounds for each center and activity
    targets_lower = Matrix{Int}(undef, S, M)
    targets_upper = Matrix{Int}(undef, S, M)

    # Calculate bounds for each center i and activity m
    for i in 1:S
        for m in 1:M
            targets_lower[i, m] = round(Int, μ[m][i] * (1 - T[m]))
            targets_upper[i, m] = round(Int, μ[m][i] * (1 + T[m]))
        end
    end

    # Convert to static matrices
    targets_lower = SMatrix{S,M,Int}(targets_lower)
    targets_upper = SMatrix{S,M,Int}(targets_upper)

    return targets_lower, targets_upper
end



function start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    for m in eachindex(V)
        mul!(view(values_matrix, :, m), X, V[m])
    end
    mul!(risk_vec, X, R)
    return values_matrix, risk_vec
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

function isFactible4(solution::Types.Solution, targets_lower, targets_upper, verbose=false)
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
    M = instance.M
    S = instance.S
    P = instance.P
    B = instance.B
    β = instance.β
    β = instance.β
    R = instance.R
    K = instance.K
    V = instance.V

    counts_k = []
    usables_i = findall(==(1), Y)
    remove = Set{Int64}()
    add = Set{Int64}()

    for y in usables_i
        assignments_y = X[y, :]
        if !any(x -> x == 1, assignments_y)
            if verbose
                println("Violando asignación de Y en: $y")
            end
            number_constraints_violated += 1
            # penalizar más aquí?
        end
    end

    if sum(Y) != P
        println(sum(Y))
        println(P)
        @error "P NO ES Y"
        number_constraints_violated += 1
    end

    for j in 1:B
        if sum(X[i, j] for i in usables_i) ≠ 1
            if verbose
                println("Violando servicio a la BU: ", j)
                number_constraints_violated += 1
            end
            #@error "X NO SERVIDA"
        end
    end

    for k_type in 1:K
        indices_k_type = Sk[k_type]
        count_k_type = 0
        for indice in indices_k_type
            if indice ∈ usables_i
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

    for i in usables_i
        for m in 1:M
            if (targets_lower[i, m] > sum(X[i, j] * V[m][j] for j in 1:B))
                if verbose
                    println("violando V inferior en i: $i y m: $m")
                    println("μ: ", targets_lower[m, i])
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                push!(add, i)
                number_constraints_violated += 1
            end
            if (sum(X[i, j] * V[m][j] for j in 1:B) > targets_upper[i, m])
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", targets_upper[m, i])
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                push!(remove, i)
                number_constraints_violated += 1
            end
        end
    end

    for i in usables_i
        if sum(X[i, j] * R[j] for j in 1:B) > β[i]
            if verbose
                println("violando riesgo en $i")
                println("β: ", β[i])
                println("R: ", sum(X[i, j] * R[j] for j in 1:B))
            end
            push!(remove, i) # por eso es set arriba, para evitar duplicados si viola las dos res de upper
            number_constraints_violated += 1
        end
    end
    if number_constraints_violated ≠ 0
        return false, number_constraints_violated, remove, add
    else
        return true, number_constraints_violated, remove, add
    end
end

# seguir agregando o quitando hasta que sea factible+

function minimums2(v, n)
    l = length(v)
    ix = [1:l;]
    partialsortperm!(ix, v, (l-n+1):l, rev=true, initialized=true)
    vals = reverse!(v[ix[(l-n+1):l]])
    indices = reverse!(CartesianIndices(v)[ix[(l-n+1):l]])
    return indices
end

function maximums(vec, n)::Tuple{Vector{Int64},Vector{CartesianIndex{1}}}
    type = eltype(vec)
    vals = zeros(type, n)
    indices = Array{Int64}(undef, n)
    arr = Array{Tuple{type,CartesianIndex}}(undef, n)
    smallest = 0
    index = 0
    @inbounds for i ∈ eachindex(vec)
        smallest, index = findmin(vals)
        if vec[i] > smallest
            arr[index] = vec[i], CartesianIndex(i)
            vals[index] = vec[i]
        end
    end
    arr = sort(arr, by=x -> x[1], rev=true)
    vals = [x[1] for x in arr]
    indices = [x[2] for x in arr]
    return vals, indices
end


function remove!(V, item)
    deleteat!(V, findall(x -> x == item, V))
end

function repair_solution1(solution, cons, targets_lower, targets_upper, remove, add)
    X = copy(solution.X)
    Y = copy(solution.Y)
    instance = solution.Instance
    B = instance.B
    S = instance.S
    D = copy(instance.D)
    M = instance.M
    V = instance.V
    R = instance.R
    β = instance.β
    P = instance.P
    n = round(Int, (P - 1))
    usables_i = findall(==(1), Y)
    not_usables_i = Set(findall(==(0), Y))
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)

    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    cant_fix = false
    for ĩ in remove
        if cant_fix
            break
        end
        values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
        # aqui el chiste es quitar las bus ASIGNADAS MAS LEJANAS
        # para asegurar que solo son las asignadas, multiplicamos la D por X
        all_asigneds = findall(==(1), X[ĩ, :])
        Dᵢ = copy(D[ĩ, :]) .* X[ĩ, :] # al hacer 0s los no asignados, no los agarra maximums
        vals, candidates_bus = maximums(Dᵢ, length(all_asigneds))
        factible_yet = false
        while !factible_yet
            factible_yet = true
            if length(candidates_bus) == 0
                cant_fix = true
                break
            end
            j = popfirst!(candidates_bus)[1] # arreglar el indice
            # cuales is podemos agarrar? las n mas cercanas? probamos todas? 
            # agarramos las n iₛ mas cercanas a j
            distance_to_bu = copy(D[:, j])
            for not_usable_i in not_usables_i
                distance_to_bu[not_usable_i] = 999999999 # a las is en Yi = 0, hacemos la distancia super grande para no fomentar su uso
            end
            j_assigned = false
            n = P
            candidates_is = minimums2(distance_to_bu, n)
            while !j_assigned
                i = 0
                if length(candidates_is) == 1
                    j_assigned = true # brincate la j   
                    break
                else
                    i = popfirst!(candidates_is)
                    i = i[1]
                end
                can_do_move = true
                if i == ĩ
                    can_do_move = false
                end
                #println("probando cambio en $i por $ĩ")

                for m in 1:M
                    values_matrix[ĩ, m] -= V[m][j] # restale a ĩ, previo
                    values_matrix[i, m] += V[m][j] # sumale a i, nuevo
                    if values_matrix[i, m] > targets_upper[i, m]
                        #println("no es factible por que se pasa el otro upper: ", values_matrix[i, m])
                        factible_yet = false
                    end
                    if values_matrix[ĩ, m] < targets_lower[i, m]
                        factible_yet = false
                        #println("no es factible por que se baja el lower: ", values_matrix[ĩ, m], " para ", targets_lower[m])
                        can_do_move = false
                        # no deberia de pasar porque entonces ĩ es infactible ahora
                    end
                end
                risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
                risk_vec[i] += R[j] # sumale a i el nuevo
                if risk_vec[i] > β[i]
                    # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                    can_do_move = false
                    #println("no es factible porque se pasa el otro risk: ", risk_vec[i])
                    factible_yet = false
                end
                if risk_vec[ĩ] > β[ĩ]
                    #println("no es factible porque se baja el risk")
                    factible_yet = false
                end
                if can_do_move
                    X[ĩ, j] = 0
                    X[i, j] = 1
                    j_assigned = true
                else
                    # si no puedo hacer el movimiento, restaura el valor de la ev parcial
                    for m in 1:M
                        values_matrix[ĩ, m] += V[m][j] # corrige el valor de la NO asignacion
                        values_matrix[i, m] -= V[m][j]
                    end
                    risk_vec[ĩ] += R[j]
                    risk_vec[i] -= R[j]
                end
            end
        end
    end
    newSol3 = Solution(instance, X, Y, solution.Weight, solution.Time)
    cant_fix = false

    for i in add
        #println("intentando $i")
        if cant_fix
            break
        end
        prev_assigned_bus = findall(==(1), X[i, :]) #indices de las bus asignadas previamente
        Dᵢ = D[i, :]
        for prev in prev_assigned_bus
            Dᵢ[prev] = 999999999
        end
        # las previamente asignadas las ponemos muy grandes, buscamos los minimoss
        candidates_bus = minimums2(Dᵢ, B) #indices de las bus mas cercanas al centro i
        factible_yet = false
        while !factible_yet
            solmamada = Solution(instance, X, Y, solution.Weight, solution.Time)
            algo, cons = isFactible(solmamada, false)
            inicializar = []
            if length(candidates_bus) == 0
                cant_fix = true
                break
            end
            j = popfirst!(candidates_bus)[1]
            ĩ = findfirst(==(1), X[:, j]) # obten la asignacion previa de cada j
            factible_yet = true
            can_do_move = true
            for m in 1:M
                if typeof(j) === Nothing
                    #println("uh con la j")
                    break
                end
                if typeof(ĩ) === Nothing
                    #println("uhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh")
                    break
                end
                values_matrix[ĩ, m] -= V[m][j] # restale a ĩ, previo
                values_matrix[i, m] += V[m][j] # sumale a i, nuevo
                if values_matrix[i, m] > targets_upper[i, m]
                    # no deberia de pasar porque entonces la infactibilidad cambia de razon
                    factible_yet = false
                    can_do_move = false
                end
                if values_matrix[i, m] < targets_lower[i, m]
                    factible_yet = false
                end
                if values_matrix[ĩ, m] < targets_lower[ĩ, m]
                    factible_yet = false
                    can_do_move = false
                    # no deberia de pasar porque entonces ĩ es infactible ahora
                end
            end
            if typeof(ĩ) === Nothing
                #println("uhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh 2")
                break
            end
            risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
            risk_vec[i] += R[j] # sumale a i el nuevo
            if risk_vec[i] > β[i]
                # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                can_do_move = false
                factible_yet = false
            end
            if can_do_move
                X[ĩ, j] = 0
                X[i, j] = 1
                #println("reasigne de $ĩ a $i")
            else
                # si no puedo hacer el movimiento, restaura el valor de la ev parcial
                for m in 1:M
                    values_matrix[ĩ, m] += V[m][j] # corrige el valor de la NO asignacion
                    values_matrix[i, m] -= V[m][j]
                end
                risk_vec[ĩ] += R[j]
                risk_vec[i] -= R[j]
            end
        end
    end

    Weight = 0
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += instance.D[indice]
    end
    newSol = Solution(instance, X, Y, Weight, solution.Time)
    #println(isFactible(newSol))
    return newSol
end


function repair_solution2(solution, cons, targets_lower, targets_upper, remove, add)
    X = copy(solution.X)
    Y = copy(solution.Y)
    instance = solution.Instance
    B = instance.B
    S = instance.S
    D = copy(instance.D)
    M = instance.M
    V = instance.V
    R = instance.R
    β = instance.β
    P = instance.P
    n = round(Int, (P - 1))
    not_usables_i = Set(findall(==(0), Y))
    usables_i = Set(findall(==(1), Y))
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    remove_vec = collect(remove) |> sort!
    remove_vec_fixed = collect(remove) |> sort!
    add_vec = collect(add) |> sort!
    add_vec_fixed = collect(add) |> sort!
    remove_factible = [false for i in remove_vec]
    add_factible = [false for i in add_vec]

    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    while !all(remove_factible) # mientras siga habiendo un falso en remover
        distances_fixed = Matrix{Int64}(undef, S, B)
        for i in 1:S
            if i ∈ remove_vec
                distances_fixed[i, :] = (X[i, :] .* D[i, :])
            else
                distances_fixed[i, :] .= 0 # cualquier BU a quitar que no pertenezca a un cliente en remove, ignora
            end
        end
        incumbent, coords = findmax(distances_fixed)
        col = coords[2]
        useful_rows = setdiff(usables_i, remove)
        distance_to_col_fixed = [i ∈ useful_rows ? D[i, col] : 10000000000000 for i in 1:S]
        minimum_possible_new_assignments, row = findmin(distance_to_col_fixed)
        can_do_move = true
        ĩ = coords[1] # vieja asignacion
        i = row
        j = col
        factible_yet = true
        for m in 1:M
            values_matrix[ĩ, m] -= V[m][j] # restale a ĩ, previo
            values_matrix[i, m] += V[m][j] # sumale a i, nuevo
            if values_matrix[i, m] > targets_upper[m]
                #println("no es factible por que se pasa el otro upper: ", values_matrix[i, m])
                factible_yet = false
            end
        end
        risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
        risk_vec[i] += R[j] # sumale a i el nuevo
        if risk_vec[i] > β[i]
            # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
            can_do_move = false
            #println("no es factible porque se pasa el otro risk: ", risk_vec[i])
            factible_yet = false
        end
        if risk_vec[ĩ] > β[ĩ]
            #println("no es factible porque se baja el risk")
            factible_yet = false
        end
        if can_do_move
            X[ĩ, j] = 0
            X[i, j] = 1
        else
            D[i, j] = 100000000000000000
            # si no puedo hacer el movimiento, restaura el valor de la ev parcial
            for m in 1:M
                values_matrix[ĩ, m] += V[m][j] # corrige el valor de la NO asignacion
                values_matrix[i, m] -= V[m][j]
            end
            risk_vec[ĩ] += R[j]
            risk_vec[i] -= R[j]
        end
        if factible_yet
            fixed_idx = findfirst(==(ĩ), remove_vec_fixed)
            remove!(remove_vec, ĩ)
            remove_factible[fixed_idx] = true
        end
    end

    sol_removed = Solution(instance, X, Y, solution.Weight, solution.Time)
    factible, constraints, remove, add = isFactible4(sol_removed, targets_lower, targets_upper)
    if factible
        return sol_removed
    end

    D = copy(instance.D)
    cant_fix = false
    escrito = false
    while !all(add_factible) # mientras siga habiendo un falso en agregr
        if cant_fix
            break
        end
        distances_fixed = Matrix{Int64}(undef, S, B)
        for i in 1:S
            if i ∈ add_vec
                prev_assigned_bus = findall(==(1), X[i, :]) #indices de las bus asignadas previamente
                distances_fixed[i, :] = D[i, :]
                distances_fixed[i, prev_assigned_bus] .= 100000000
            else
                distances_fixed[i, :] .= 1000000000 # cualquier BU a quitar que pertenezca a un centro en add, ignora
            end
        end
        incumbents_vec, coords_vec = minimums(distances_fixed, B)
        assigned = false
        while !assigned
            if length(coords_vec) == 0
                cant_fix = true
                break
            end
            coords = popfirst!(coords_vec)
            col = coords[2]
            ĩ = findfirst(==(1), X[:, col])
            row = coords[1]
            can_do_move = true
            i = row
            j = col
            factible_yet = true
            razon = 0
            for m in 1:M
                values_matrix[ĩ, m] -= V[m][j] # restale a ĩ, previo
                values_matrix[i, m] += V[m][j] # sumale a i, nuevo
                if values_matrix[i, m] > targets_upper[i, m]
                    #println("no es factible por que se pasa el otro upper: ", values_matrix[i, m])
                    factible_yet = false
                end
                if values_matrix[ĩ, m] < targets_lower[ĩ, m]
                    factible_yet = false
                    #println("no es factible por que se baja el lower: ", values_matrix[ĩ, m], " para ", targets_lower[m])
                    can_do_move = false
                    razon = 1
                    # no deberia de pasar porque entonces ĩ es infactible ahora
                end
                if values_matrix[i, m] < targets_lower[i, m]
                    factible_yet = false
                end
            end
            risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
            risk_vec[i] += R[j] # sumale a i el nuevo
            if risk_vec[i] > β[i]
                # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                can_do_move = false
                razon = 2
                #println("no es factible porque se pasa el otro risk: ", risk_vec[i])
                factible_yet = false
            end
            if risk_vec[ĩ] > β[ĩ]
                #println("no es factible porque se baja el risk")
                factible_yet = false
            end
            if can_do_move
                X[ĩ, j] = 0
                X[i, j] = 1
                assigned = true
            else
                # si no puedo hacer el movimiento, restaura el valor de la ev parcial
                for m in 1:M
                    values_matrix[ĩ, m] += V[m][j] # corrige el valor de la NO asignacion
                    values_matrix[i, m] -= V[m][j]
                end
                risk_vec[ĩ] += R[j]
                risk_vec[i] -= R[j]
            end
            if factible_yet
                fixed_idx = findfirst(==(i), add_vec_fixed)
                remove!(add_vec, i)
                add_factible[fixed_idx] = true
            end
        end
    end

    Weight = 0
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += instance.D[indice]
    end
    newSol = Solution(instance, X, Y, Weight, solution.Time)
    # println(isFactible(newSol, true))
    return newSol
end

function repair_solution4(solution, cons, targets_lower, targets_upper, remove, add)
    X = copy(solution.X)
    Y = copy(solution.Y)
    instance = solution.Instance
    S = instance.S
    M = 3
    B = instance.B
    V = instance.V
    R = instance.R
    β = instance.β
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    
    # Handle facilities that need to remove BUs
    for ĩ in remove
        assigned_bus = findall(==(1), X[ĩ, :])
        
        for j in assigned_bus  # Try every assigned BU
            # Try every possible facility (ignoring distance)
            for i in 1:S
                i == ĩ && continue
                Y[i] == 0 && continue
                
                # Test move feasibility
                can_do_move = true
                
                # Simulate move
                for m in 1:M
                    values_matrix[ĩ, m] -= V[m][j]
                    values_matrix[i, m] += V[m][j]
                    
                    # Check bounds
                    if values_matrix[i, m] > targets_upper[i, m] ||
                        values_matrix[i, m] < targets_lower[i, m] ||
                        values_matrix[ĩ, m] > targets_upper[ĩ, m] ||
                        values_matrix[ĩ, m] < targets_lower[ĩ, m]
                         can_do_move = false
                         break
                     end
                end
                
                # Check risk
                temp_risk_i = risk_vec[i] + R[j]
                temp_risk_ĩ = risk_vec[ĩ] - R[j]
                
                if temp_risk_i > β[i] || temp_risk_ĩ > β[ĩ]
                    can_do_move = false
                end
                
                if can_do_move
                    # Execute move
                    X[ĩ, j] = 0
                    X[i, j] = 1
                    risk_vec[ĩ] -= R[j]
                    risk_vec[i] += R[j]
                    break  # Move to next BU
                else
                    # Restore values if move not executed
                    for m in 1:M
                        values_matrix[ĩ, m] += V[m][j]
                        values_matrix[i, m] -= V[m][j]
                    end
                end
            end
        end
    end
    
    # Handle facilities that need BUs
    for i in add
        min_progress = minimum(values_matrix[i, m] / targets_lower[i, m] for m in 1:M)
        
        # Find critical measures (those furthest from target)
        critical_measures = findall(m -> values_matrix[i, m] / targets_lower[i, m] <= min_progress + 0.05, 1:M)
        
        # Sort BUs by their contribution to critical measures
        all_bus = collect(1:B)
        sort!(all_bus, by=j -> begin
            # Calculate weighted contribution to critical measures
            contribution = sum(V[m][j] for m in critical_measures) / length(critical_measures)
            return contribution
        end, rev=true)  # Highest contributors first
        
        for j in all_bus
            ĩ = findfirst(==(1), X[:, j])  # Current assignment
            ĩ === nothing && continue
            if ĩ === i
                continue
            end
            
            # First check if removing from current center would maintain feasibility
            would_maintain_feasibility = true
            for m in 1:M
                new_value_ĩ = values_matrix[ĩ, m] - V[m][j]
                if new_value_ĩ < targets_lower[ĩ, m]
                    would_maintain_feasibility = false
                    break
                end
            end
            
            if !would_maintain_feasibility
                continue  # Skip this BU if it would make current center infeasible
            end
            
            # Now test move feasibility for target center
            can_do_move = true
            
            # Temporary arrays to track changes
            temp_values_i = copy(values_matrix[i, :])
            temp_values_ĩ = copy(values_matrix[ĩ, :])
            
            # Simulate move
            for m in 1:M
                temp_values_ĩ[m] -= V[m][j]
                temp_values_i[m] += V[m][j]
                
                if temp_values_i[m] > targets_upper[i, m] ||
                   temp_values_i[m] < targets_lower[i, m] ||
                   temp_values_ĩ[m] > targets_upper[ĩ, m] ||
                   temp_values_ĩ[m] < targets_lower[ĩ, m]
                    can_do_move = false
                    break
                end
            end
            
            # Check risk constraints
            temp_risk_i = risk_vec[i] + R[j]
            temp_risk_ĩ = risk_vec[ĩ] - R[j]
            
            if temp_risk_i > β[i] || temp_risk_ĩ > β[ĩ]
                can_do_move = false
            end
            
            if can_do_move
                # Execute move and update all matrices
                X[ĩ, j] = 0
                X[i, j] = 1
                risk_vec[ĩ] -= R[j]
                risk_vec[i] += R[j]
                for m in 1:M
                    values_matrix[ĩ, m] = temp_values_ĩ[m]
                    values_matrix[i, m] = temp_values_i[m]
                end
                
                # Check if we've met the minimum targets
                if all(values_matrix[i, m] >= targets_lower[i, m] for m in 1:M)
                    break  # Stop adding BUs if we've met targets
                end
            end
        end
    end
    
    Weight = sum(instance.D[i,j] for i in 1:S, j in 1:B if X[i,j] == 1)

    #swap_problematic_centers!(X, Y, values_matrix, risk_vec, targets_lower, targets_upper, instance)
    return Solution(instance, X, Y, Weight, solution.Time)
end


function repair_solution_simplified(solution, cons, targets_lower, targets_upper, remove, add)
    X = copy(solution.X)
    Y = copy(solution.Y)
    instance = solution.Instance
    S = instance.S
    M = instance.M
    B = instance.B
    V = instance.V
    R = instance.R
    β = instance.β
    D = instance.D
    
    # Initialize and maintain values_matrix for efficient constraint checking
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    
    # Handle facilities that need to remove BUs
    for ĩ in remove
        assigned_bus = findall(==(1), X[ĩ, :])
        for j in assigned_bus
            # Try every other active facility
            candidates = collect(1:S)
            filter!(i -> i != ĩ && Y[i] == 1, candidates)
            
            for i in candidates
                # Test move feasibility using values_matrix
                can_do_move = true
                
                # Try the move
                for m in 1:M
                    values_matrix[ĩ, m] -= V[m][j]
                    values_matrix[i, m] += V[m][j]
                    
                    if values_matrix[i, m] > targets_upper[i, m] ||
                       values_matrix[i, m] < targets_lower[i, m] ||
                       values_matrix[ĩ, m] > targets_upper[ĩ, m] ||
                       values_matrix[ĩ, m] < targets_lower[ĩ, m]
                        can_do_move = false
                        break
                    end
                end
                
                # Check risk
                temp_risk_i = risk_vec[i] + R[j]
                temp_risk_ĩ = risk_vec[ĩ] - R[j]
                
                if temp_risk_i > β[i] || temp_risk_ĩ > β[ĩ]
                    can_do_move = false
                end
                
                if can_do_move
                    # Execute move
                    X[ĩ, j] = 0
                    X[i, j] = 1
                    risk_vec[ĩ] -= R[j]
                    risk_vec[i] += R[j]
                    break  # Found a valid move, go to next BU
                else
                    # Restore values if move not possible
                    for m in 1:M
                        values_matrix[ĩ, m] += V[m][j]
                        values_matrix[i, m] -= V[m][j]
                    end
                end
            end
        end
    end
    
    # Handle facilities that need BUs
    for i in add
        # Try all BUs
        for j in 1:B
            ĩ = findfirst(==(1), X[:, j])
            ĩ === nothing && continue
            ĩ == i && continue
            
            # Test move feasibility
            can_do_move = true
            
            # Try the move
            for m in 1:M
                values_matrix[ĩ, m] -= V[m][j]
                values_matrix[i, m] += V[m][j]
                
                if values_matrix[i, m] > targets_upper[i, m] ||
                   values_matrix[i, m] < targets_lower[i, m] ||
                   values_matrix[ĩ, m] > targets_upper[ĩ, m] ||
                   values_matrix[ĩ, m] < targets_lower[ĩ, m]
                    can_do_move = false
                    break
                end
            end
            
            # Check risk
            temp_risk_i = risk_vec[i] + R[j]
            temp_risk_ĩ = risk_vec[ĩ] - R[j]
            
            if temp_risk_i > β[i] || temp_risk_ĩ > β[ĩ]
                can_do_move = false
            end
            
            if can_do_move
                # Execute move
                X[ĩ, j] = 0
                X[i, j] = 1
                risk_vec[ĩ] -= R[j]
                risk_vec[i] += R[j]
                
                # Check if we've met all constraints
                if all(values_matrix[i, m] >= targets_lower[i, m] for m in 1:M)
                    break
                end
            else
                # Restore values if move not possible
                for m in 1:M
                    values_matrix[ĩ, m] += V[m][j]
                    values_matrix[i, m] -= V[m][j]
                end
            end
        end
    end
    
    Weight = sum(D[i,j] for i in 1:S, j in 1:B if X[i,j] == 1)
    return Solution(instance, X, Y, Weight, solution.Time)
end

# Optimized Interchange Move Implementation

function can_do_interchange_optimized(values_matrix::Matrix{Int}, V::Vector{Vector{Int}}, ĩ::Int, i✶::Int, j₁::Int, j₂::Int, risk_vec::Vector{Int}, R::Vector{Int}, β::Vector{Int}, targets_upper, targets_lower)
    @inbounds for m in 1:3
        value_ĩ = values_matrix[ĩ, m] - V[m][j₁] + V[m][j₂]
        value_i✶ = values_matrix[i✶, m] + V[m][j₁] - V[m][j₂]
        if value_ĩ > targets_upper[ĩ, m] || value_ĩ < targets_lower[ĩ, m] ||
           value_i✶ > targets_upper[i✶, m] || value_i✶ < targets_lower[i✶, m]
            return false
        end
    end
    risk_ĩ = risk_vec[ĩ] - R[j₁] + R[j₂]
    risk_i✶ = risk_vec[i✶] + R[j₁] - R[j₂]
    return risk_ĩ <= β[ĩ] && risk_i✶ <= β[i✶]
end

function find_best_interchange_optimized(X::Matrix{Int}, D::Matrix{Int}, V::Vector{Vector{Int}}, R::Vector{Int}, B::Int, M::Int, values_matrix::Matrix{Int}, risk_vec::Vector{Int}, targets_lower, targets_upper, β::Vector{Int}, strategy::Symbol)
    best_move = nothing
    best_weight_diff = 0.0

    @inbounds for j₁ in 1:B-1
        ĩ = find_one_in_column_unrolled(X, j₁)
        for j₂ in j₁+1:B
            i✶ = find_one_in_column_unrolled(X, j₂)
            weight_diff = D[i✶, j₁] + D[ĩ, j₂] - D[ĩ, j₁] - D[i✶, j₂]

            if weight_diff < 0 && can_do_interchange_optimized(values_matrix, V, ĩ, i✶, j₁, j₂, risk_vec, R, β, targets_upper, targets_lower)
                if strategy == :ff
                    return (ĩ=ĩ, i✶=i✶, j₁=j₁, j₂=j₂, weight_diff=weight_diff)
                elseif strategy == :bf && (best_move === nothing || weight_diff < best_weight_diff)
                    best_move = (ĩ=ĩ, i✶=i✶, j₁=j₁, j₂=j₂, weight_diff=weight_diff)
                    best_weight_diff = weight_diff
                end
            end
        end
        strategy == :ff && best_move !== nothing && break
    end
    return best_move
end

function apply_interchange_optimized!(X::Matrix{Int}, values_matrix::Matrix{Int}, risk_vec::Vector{Int}, move, V::Vector{Vector{Int}}, R::Vector{Int})
    ĩ, i✶, j₁, j₂ = move.ĩ, move.i✶, move.j₁, move.j₂
    @inbounds X[ĩ, j₁], X[i✶, j₂] = 0, 0
    @inbounds X[ĩ, j₂], X[i✶, j₁] = 1, 1
    @inbounds for m in 1:size(values_matrix, 2)
        values_matrix[ĩ, m] += V[m][j₂] - V[m][j₁]
        values_matrix[i✶, m] += V[m][j₁] - V[m][j₂]
    end
    @inbounds risk_vec[ĩ] += R[j₂] - R[j₁]
    @inbounds risk_vec[i✶] += R[j₁] - R[j₂]
end

function interchange_bu_improve_optimized(solution, targets_lower, targets_upper, strategy::Symbol)
    X, Y = solution.X, solution.Y
    instance = solution.Instance
    B, S, M = instance.B, instance.S, instance.M
    V, R, β = instance.V, instance.R, instance.β
    D = instance.D
    Weight = solution.Weight

    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)

    while true
        best_move = find_best_interchange_optimized(X, D, V, R, B, M, values_matrix, risk_vec, targets_lower, targets_upper, β, strategy)
        best_move === nothing && break
        apply_interchange_optimized!(X, values_matrix, risk_vec, best_move, V, R)
        Weight += best_move.weight_diff
    end

    Weight = dot(X, D)  # Recalculate weight to ensure accuracy
    #println("INTERCHANGE")
    #sol = Solution(instance, X, Y, Weight, solution.Time)
    #println(isFactible(sol))
    return Solution(instance, X, Y, Weight, solution.Time)
end

function minimums(matrix::Matrix, n)::Tuple{Vector{Int64},Vector{CartesianIndex{2}}}
    type = eltype(matrix)
    vals = fill(10000000000000, n)
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

function get_best_assignments(D, N)
    num_centers, num_clients = size(D)
    best_assignments = Dict{Int,Vector{Int64}}()
    for j in 1:num_clients
        # Use a temporary array to store the facility opportunity costs for this client
        costs = Tuple{Int64,Int}[]
        for i in 1:num_centers
            push!(costs, (D[i, j], i))
        end
        # Sort the costs
        sort!(costs)
        # Store the top N assignments for this client
        best_assignments[j] = best_assignments[j] = [cost[2] for cost in costs[1:N]] # extrae el indice nada mas
    end
    return best_assignments
end

function get_best_clients_for_centers(D::Matrix{Int64}, N::Int64)
    num_centers, num_clients = size(D)
    best_clients = Dict{Int,Vector{Int64}}()
    costs = Vector{Tuple{Int64,Int}}(undef, num_clients)
    @inbounds for i in 1:num_centers
        # Use a temporary array to store the client opportunity costs for this center
        # costs = Tuple{Int64,Int}[]
        @inbounds for j in 1:num_clients
            #push!(costs, (D[i, j], j))
            costs[j] = (D[i, j], j)
        end
        # Sort the costs
        sort!(costs)
        # Store the top N clients for this center
        best_clients[i] = Vector{Int64}(undef, N)
        @inbounds for k in 1:N
            best_clients[i][k] = costs[k][2]
        end
        # best_clients[i] = [cost[2] for cost in costs[1:N]] # Extract only the client index
    end
    return best_clients
end
#=
"""
Pre-Processing:

    Initialization: Start with a given solution (feasible or not).

    Compute Data Structures for Quick Lookups:
        Using the distance matrix DD, compute the get_best_assignments function which provides the top NN centers for each client.
        Similarly, calculate the get_best_clients_for_centers which lists the top NN clients for each center.

Local Search:

    While there is potential for improvement:

    a. Iterate Over Centers:
        For every currently active center ı~ı~:
            Consider deactivating this center, making its clients "orphaned".
            For every currently inactive center i✶i✶:
                Consider activating this center.

    b. Fulfilling New Centers:
        Use the get_best_clients_for_centers to get the nearest clients for i✶i✶.
        Assign clients to i✶i✶ until it becomes feasible or until there's no more benefit in adding more clients.
        Ensure that removing a client from another center doesn't make the original center unfeasible.

    c. Reassign Orphaned Clients:
        For each orphaned client (previously assigned to ı~ı~):
            Use the get_best_assignments data structure to find the most suitable center for this client.
            Reassign the client to a feasible center, ensuring the center remains feasible post-assignment.

    d. Evaluate Solution:
        After attempting to activate a new center and reassigning clients, evaluate the new solution.
        If the new solution is better (e.g., has a lower total distance or better meets other criteria), keep it; otherwise, revert to the previous solution.

    If a better solution is found during an iteration, continue the loop. Otherwise, terminate the local search.
# Pre-Processing
initial_solution = ...  # Your initial solution

# Compute quick lookup structures
best_assignments_for_clients = get_best_assignments(D, Y, N)
best_clients_for_centers = get_best_clients_for_centers(D, Y, N)

best_solution = initial_solution
improvement = true

# Local Search
while improvement
    improvement = false

    for each active_center ĩ in best_solution
        for each inactive_center i✶ in best_solution

            # Deactivate current center ĩ
            orphaned_clients = deactivate_center(ĩ)

            # Activate new center i✶ and assign its nearest clients
            assign_nearest_clients(i✶, best_clients_for_centers)

            # Reassign orphaned clients
            for client in orphaned_clients
                best_new_center = find_best_center_for_client(client, best_assignments_for_clients)
                if best_new_center is feasible after assignment
                    assign(client, best_new_center)
                end
            end

            # Evaluate the new solution
            new_solution = compute_new_solution(...)  # With the changes made
            if is_better(new_solution, best_solution)
                best_solution = new_solution
                improvement = true
            else
                # If the new solution isn't better, revert changes
                revert_changes()
            end
        end
    end
end

# Post-Processing
return best_solution
"""
=#

function deactivate_center_improve(solution, targets_lower, targets_upper, strategy=:bf)
    X, Y = solution.X, solution.Y
    instance = solution.Instance
    B, S, M, V, R, β, P = instance.B, instance.S, instance.M, instance.V, instance.R, instance.β, instance.P
    # D = copy(instance.D)
    D = instance.D
    Sk = instance.Sk
    K = 5
    Lk = instance.Lk
    Uk = instance.Uk
    Weight = solution.Weight
    not_usables_i = BitSet(findall(==(0), Y))
    usables_i = BitSet(findall(==(1), Y))
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    best_assignments_clients = get_best_assignments(D, P)
    best_clients_for_centers = get_best_clients_for_centers(D, B)
    count = count_k(usables_i, Sk)
    improvement = true
    total_time_spent = 0
    while improvement
        #while improvement
        improvement = false
        for ĩ in usables_i
            for i✶ in not_usables_i
                modified_X = Dict{Tuple{Int,Int},Int}()
                modified_values_matrix = Dict{Tuple{Int,Int},Int}()
                modified_risk = Dict{Int,Int}()
                useful = true
                ĩₖ = node_type(ĩ, Sk)
                i✶ₖ = node_type(i✶, Sk)
                count_ĩ = count[ĩₖ] - 1
                count_i✶ = count[i✶ₖ] + 1

                if count_ĩ <= Uk[ĩₖ] && count_ĩ >= Lk[ĩₖ] && count_i✶ <= Uk[i✶ₖ] && count_i✶ >= Lk[i✶ₖ]
                    #js_assigned = findall(==(1), @view X[ĩ, :])
                    js_assigned = find_ones_in_row_optimized(X, ĩ)
                    for j in js_assigned
                        if !haskey(modified_X, (ĩ, j))
                            modified_X[(ĩ, j)] = X[ĩ, j]
                        end
                        for m in 1:M
                            if !haskey(modified_values_matrix, (ĩ, m))
                                modified_values_matrix[(ĩ, m)] = values_matrix[ĩ, m]
                            end
                            values_matrix[ĩ, m] = 0
                        end
                        if !haskey(modified_risk, ĩ)
                            modified_risk[ĩ] = risk_vec[ĩ]
                        end
                        risk_vec[ĩ] = 0
                    end
                    X[ĩ, :] .= 0
                    # aqui se modifica la X entonces debemos de restar los vawlores de values_matrix y risk_vec a las js asignadas
                    js_assigned_set = Set(js_assigned)
                    weight_old_branch = sum(D[ĩ, js_assigned]) # total weight of old branch/center
                    weight_new_branch = 0
                    factible_yet = false
                    fulls_m = zeros(Int, M)
                    for client in best_clients_for_centers[i✶]
                        #for client in best_clients_for_centers[i✶]
                        potential_assignment_valid = true
                        #previous_i_client = findfirst(==(1), @views X[:, client])
                        previous_i_client = find_one_in_column_unrolled(X, client)
                        if previous_i_client !== 0 # si es nothing entonces estaba asignado a ĩ
                            for m in 1:M
                                if values_matrix[previous_i_client, m] - V[m][client] < targets_lower[previous_i_client, m]
                                    potential_assignment_valid = false
                                    break
                                end
                            end
                            if risk_vec[i✶] + R[client] > β[i✶]
                                potential_assignment_valid = false
                            end
                            for m in 1:M
                                if values_matrix[i✶, m] + V[m][client] > targets_upper[i✶, m]
                                    potential_assignment_valid = false
                                    break
                                end
                            end
                            if potential_assignment_valid
                                for m in 1:M
                                    ##=
                                    if !haskey(modified_values_matrix, (previous_i_client, m))
                                        modified_values_matrix[(previous_i_client, m)] = values_matrix[previous_i_client, m]
                                    end
                                    if !haskey(modified_values_matrix, (i✶, m))
                                        modified_values_matrix[(i✶, m)] = values_matrix[i✶, m]
                                    end
                                    #
                                    values_matrix[previous_i_client, m] -= V[m][client]
                                    values_matrix[i✶, m] += V[m][client]
                                    if values_matrix[i✶, m] > targets_lower[i✶, m]
                                        fulls_m[m] = 1
                                    end
                                end
                                ##=
                                if !haskey(modified_risk, previous_i_client)
                                    modified_risk[previous_i_client] = risk_vec[previous_i_client]
                                end

                                if !haskey(modified_risk, i✶)
                                    modified_risk[i✶] = risk_vec[i✶]
                                end
                                #
                                risk_vec[previous_i_client] -= R[client]
                                risk_vec[i✶] += R[client]
                            end
                        else
                            for m in 1:M
                                if values_matrix[i✶, m] + V[m][client] > targets_upper[i✶, m]
                                    potential_assignment_valid = false
                                    break
                                end
                            end
                            if risk_vec[i✶] + R[client] > β[i✶]
                                potential_assignment_valid = false
                            end
                            if potential_assignment_valid
                                for m in 1:M
                                    ##=
                                    if !haskey(modified_values_matrix, (i✶, m))
                                        modified_values_matrix[(i✶, m)] = values_matrix[i✶, m]
                                    end
                                    #
                                    values_matrix[i✶, m] += V[m][client]
                                    if values_matrix[i✶, m] > targets_lower[i✶, m]
                                        fulls_m[m] = 1
                                    end
                                end
                            end
                            ##=
                            if !haskey(modified_risk, i✶)
                                modified_risk[i✶] = risk_vec[i✶]
                            end
                            #
                            risk_vec[i✶] += R[client]
                        end
                        if potential_assignment_valid
                            if !haskey(modified_X, (i✶, client))
                                modified_X[(i✶, client)] = X[i✶, client]
                            end
                            if previous_i_client !== 0
                                if !haskey(modified_X, (previous_i_client, client))
                                    modified_X[(previous_i_client, client)] = X[previous_i_client, client]
                                end
                                X[previous_i_client, client] = 0
                            end
                            X[i✶, client] = 1
                            if client in js_assigned_set
                                delete!(js_assigned_set, client)
                            end
                            weight_new_branch += D[i✶, client]
                            if all(x -> x == 1, fulls_m)
                                factible_yet = true # ya llenamos el centro i✶
                                break
                            end
                        end
                    end
                    useful = true
                    if !factible_yet
                        useful = false
                    end
                    for orphaned_client in js_assigned_set
                        assigned_yet = false
                        for center in best_assignments_clients[orphaned_client]
                            if (center ∈ usables_i && center ≠ ĩ) || center == i✶
                                potential_assignment_valid = true
                                for m in 1:M
                                    if values_matrix[center, m] + V[m][orphaned_client] > targets_upper[center, m]
                                        potential_assignment_valid = false
                                        break
                                    end
                                end
                                if risk_vec[center] + R[orphaned_client] > β[center]
                                    potential_assignment_valid = false
                                end
                                if potential_assignment_valid
                                    if !haskey(modified_X, (center, orphaned_client))
                                        modified_X[(center, orphaned_client)] = X[center, orphaned_client]
                                    end
                                    X[center, orphaned_client] = 1
                                    assigned_yet = true
                                    weight_new_branch += D[center, orphaned_client]
                                    for m in 1:M
                                        ##=
                                        if !haskey(modified_values_matrix, (center, m))
                                            modified_values_matrix[(center, m)] = values_matrix[center, m]
                                        end
                                        #
                                        values_matrix[center, m] += V[m][orphaned_client]
                                    end
                                    ##=
                                    if !haskey(modified_risk, center)
                                        modified_risk[center] = risk_vec[center]
                                    end
                                    #
                                    risk_vec[center] += R[orphaned_client]
                                    break
                                end
                            end
                        end
                        if !assigned_yet
                            useful = false
                            break
                        end
                    end
                    if useful && weight_new_branch < weight_old_branch
                        Y[ĩ] = 0
                        Y[i✶] = 1
                        improvement = true
                        delete!(not_usables_i, i✶)
                        delete!(usables_i, ĩ)
                        push!(usables_i, i✶)
                        push!(not_usables_i, ĩ)
                        Weight = sum(X .* D)
                        count[ĩₖ] -= 1
                        count[i✶ₖ] += 1
                    else
                        for item in modified_X
                            i, j = item[1]
                            val = item[2]
                            X[i, j] = val
                        end
                        for (center, m) in keys(modified_values_matrix)
                            values_matrix[center, m] = modified_values_matrix[(center, m)]
                        end
                        for center in keys(modified_risk)
                            risk_vec[center] = modified_risk[center]
                        end
                        #
                        # Rollback de todos los cambios
                    end
                end
            end
        end
    end
    Weight = dot(X, D)
    #@show total_time_spent
    #println("DEACTIVATE")
    #sol = Solution(instance, X, Y, Weight, solution.Time)
    #println(isFactible(sol))
    return Solution(instance, X, Y, Weight, solution.Time)
end


function can_do_move_simple_optimized(values_matrix::Matrix{Int}, V::Vector{Vector{Int64}}, i::Int, ĩ::Int, j::Int, risk_vec::Vector{Int}, R::Vector{Int64}, β::Vector{Int}, targets_upper, targets_lower)
    @inbounds for m in 1:3
        valueĩ = values_matrix[ĩ, m] - V[m][j]
        value_i = values_matrix[i, m] + V[m][j]
        if valueĩ > targets_upper[ĩ, m] || valueĩ < targets_lower[ĩ, m] || value_i > targets_upper[i, m] || value_i < targets_lower[i, m]
            return false
        end
    end
    riskĩ = risk_vec[ĩ] - R[j]
    risk_i = risk_vec[i] + R[j]
    return riskĩ <= β[ĩ] && risk_i <= β[i]
end
function simple_bu_improve_optimized(solution, targets_lower, targets_upper, strategy::Symbol)
    X, Y = solution.X, solution.Y
    instance = solution.Instance
    B, S, M = instance.B, instance.S, instance.M
    V, R, β = instance.V, instance.R, instance.β
    D = instance.D
    Weight = solution.Weight
    usables_i = findall(==(1), Y)
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    while true
        best_move = find_best_move_simple_optimized(X, D, V, R, B, M, usables_i, values_matrix, risk_vec, targets_lower, targets_upper, β, strategy)
        best_move === nothing && break
        apply_move_simple_optimized!(X, values_matrix, risk_vec, best_move, V, R)
        Weight += best_move.weight_diff
    end
    Weight = dot(X, D)
    #println("SIMPLE ")
    #sol = Solution(instance, X, Y, Weight, solution.Time)
    #println(isFactible(sol))
    return Solution(instance, X, Y, Weight, solution.Time)
end
function find_best_move_simple_optimized(X::Matrix{Int}, D::Matrix{Int}, V::Vector{Vector{Int}}, R::Vector{Int}, B::Int, M::Int, usables_i::Vector{Int}, values_matrix::Matrix{Int}, risk_vec::Vector{Int}, targets_lower, targets_upper, β::Vector{Int}, strategy::Symbol)
    best_move = nothing
    best_weight_diff = 0.0

    @inbounds for j in 1:B
        ĩ = find_one_in_column_unrolled(X, j)
        min_index = argmin(@view D[:, j])
        ĩ == min_index && continue

        for i in usables_i
            i == ĩ && continue
            weight_diff = D[i, j] - D[ĩ, j]

            if weight_diff < 0 && can_do_move_simple_optimized(values_matrix, V, i, ĩ, j, risk_vec, R, β, targets_upper, targets_lower)
                if strategy == :ff
                    return (new_i=i, old_i=ĩ, j=j, weight_diff=weight_diff)
                elseif strategy == :bf && (best_move === nothing || weight_diff < best_weight_diff)
                    best_move = (new_i=i, old_i=ĩ, j=j, weight_diff=weight_diff)
                    best_weight_diff = weight_diff
                end
            end
        end
        strategy == :ff && best_move !== nothing && break
    end
    return best_move
end
function apply_move_simple_optimized!(X::Matrix{Int}, values_matrix::Matrix{Int}, risk_vec::Vector{Int}, move, V::Vector{Vector{Int}}, R::Vector{Int})
    new_i, old_i, j = move.new_i, move.old_i, move.j
    @inbounds X[old_i, j] = 0
    @inbounds X[new_i, j] = 1
    @inbounds for m in 1:size(values_matrix, 2)
        values_matrix[old_i, m] -= V[m][j]
        values_matrix[new_i, m] += V[m][j]
    end
    @inbounds risk_vec[old_i] -= R[j]
    @inbounds risk_vec[new_i] += R[j]
end

# Chain Move Implementation: Moves k clients in a chain between open facilities
function can_do_chain_move_optimized(values_matrix::Matrix{Int}, V::Vector{Vector{Int}},
    moves::Vector{NamedTuple{(:from_facility, :to_facility, :client),
        Tuple{Int,Int,Int}}}, risk_vec::Vector{Int},
    R::Vector{Int}, β::Int, targets_upper::SVector{3,Float32},
    targets_lower::SVector{3,Float32})
    # Create temporary copies to check feasibility
    temp_values = copy(values_matrix)
    temp_risks = copy(risk_vec)

    # Apply all moves to temporary matrices
    @inbounds for move in moves
        from_facility = move.from_facility
        to_facility = move.to_facility
        client = move.client

        # Update facility values
        for m in 1:3
            temp_values[from_facility, m] -= V[m][client]
            temp_values[to_facility, m] += V[m][client]

            # Check constraints after each update
            if temp_values[from_facility, m] > targets_upper[m] ||
               temp_values[from_facility, m] < targets_lower[m] ||
               temp_values[to_facility, m] > targets_upper[m] ||
               temp_values[to_facility, m] < targets_lower[m]
                return false
            end
        end

        # Update risks
        temp_risks[from_facility] -= R[client]
        temp_risks[to_facility] += R[client]

        # Check risk constraints
        if temp_risks[from_facility] > β || temp_risks[to_facility] > β
            return false
        end
    end
    return true
end

function find_best_chain_move_optimized(X::Matrix{Int}, D::Matrix{Int}, V::Vector{Vector{Int}},
    R::Vector{Int}, B::Int, M::Int, usables_i::Vector{Int},
    values_matrix::Matrix{Int}, risk_vec::Vector{Int},
    targets_lower::SVector{3,Float32}, targets_upper::SVector{3,Float32},
    β::Int, chain_length::Int, strategy::Symbol)
    best_move = nothing
    best_weight_diff = 0.0

    # Pre-allocate vector for building chains
    current_chain = Vector{NamedTuple{(:from_facility, :to_facility, :client),
        Tuple{Int,Int,Int}}}(undef, chain_length)

    # Set to track used facilities and clients in current chain
    used_clients = BitSet(B)

    function try_extend_chain(chain_pos::Int, weight_so_far::Float64)
        # If chain is complete, check if it's the best so far
        if chain_pos > chain_length
            if can_do_chain_move_optimized(values_matrix, V, current_chain, risk_vec,
                R, β, targets_upper, targets_lower)
                if strategy == :ff
                    best_move = copy(current_chain)
                    return true
                elseif strategy == :bf && (best_move === nothing ||
                                           weight_so_far < best_weight_diff)
                    best_move = copy(current_chain)
                    best_weight_diff = weight_so_far
                end
            end
            return false
        end

        # Get previous facility in chain (or all facilities if starting)
        prev_facility = chain_pos == 1 ? nothing : current_chain[chain_pos-1].to_facility

        # Try each client
        @inbounds for client in 1:B
            client in used_clients && continue

            # Get current facility of client
            from_facility = find_one_in_column_unrolled(X, client)

            # Try each possible destination facility
            for to_facility in usables_i
                # Skip if same facility or invalid connection
                (to_facility == from_facility ||
                 (prev_facility !== nothing && from_facility != prev_facility)) && continue

                # Calculate weight difference for this move
                move_weight = D[to_facility, client] - D[from_facility, client]
                new_weight = weight_so_far + move_weight

                # Only continue if potentially improving
                if new_weight < 0
                    # Add move to chain
                    current_chain[chain_pos] = (from_facility=from_facility,
                        to_facility=to_facility, client=client)
                    push!(used_clients, client)

                    # Recursively try to complete chain
                    if try_extend_chain(chain_pos + 1, new_weight)
                        return true
                    end

                    # Backtrack
                    delete!(used_clients, client)
                end
            end
        end
        return false
    end

    # Try to build chain
    try_extend_chain(1, 0.0)
    return best_move
end

function apply_chain_move_optimized!(X::Matrix{Int}, values_matrix::Matrix{Int},
    risk_vec::Vector{Int}, moves::Vector{NamedTuple{
        (:from_facility, :to_facility, :client),Tuple{Int,Int,Int}}},
    V::Vector{Vector{Int}}, R::Vector{Int})
    # First remove all clients from their original facilities
    @inbounds for move in moves
        X[move.from_facility, move.client] = 0
    end

    # Then assign all clients to their new facilities
    @inbounds for move in moves
        X[move.to_facility, move.client] = 1

        # Update facility values
        for m in 1:size(values_matrix, 2)
            values_matrix[move.from_facility, m] -= V[m][move.client]
            values_matrix[move.to_facility, m] += V[m][move.client]
        end

        # Update risks
        risk_vec[move.from_facility] -= R[move.client]
        risk_vec[move.to_facility] += R[move.client]
    end
end

function chain_move_improve_optimized(solution, targets_lower::SVector{3,Float32},
    targets_upper::SVector{3,Float32}, strategy::Symbol,
    chain_length::Int)
    X, Y = solution.X, solution.Y
    instance = solution.Instance
    B, S, M = instance.B, instance.S, instance.M
    V, R, β = instance.V, instance.R, instance.β[1]
    D = instance.D
    Weight = solution.Weight

    # Get list of open facilities
    usables_i = findall(==(1), Y)

    # Initialize matrices for constraint checking
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)

    while true
        best_move = find_best_chain_move_optimized(X, D, V, R, B, M, usables_i, values_matrix,
            risk_vec, targets_lower, targets_upper, β,
            chain_length, strategy)
        best_move === nothing && break

        # Calculate total weight difference
        weight_diff = 0.0
        @inbounds for move in best_move
            weight_diff += D[move.to_facility, move.client] - D[move.from_facility, move.client]
        end

        apply_chain_move_optimized!(X, values_matrix, risk_vec, best_move, V, R)
        Weight += weight_diff
    end

    Weight = dot(X, D)  # Recalculate weight to ensure accuracy
    return Solution(instance, X, Y, Weight, solution.Time)
end

# Ejection Chain Move: Temporarily remove clients and find optimal reassignments
struct EjectionMove
    client::Int
    from_facility::Int
    to_facility::Int
    weight_diff::Float64
end

function can_assign_to_facility(client::Int, facility::Int, values_matrix::Matrix{Int},
    V::Vector{Vector{Int}}, risk_vec::Vector{Int}, R::Vector{Int},
    β::Int, targets_upper::SVector{3,Float32},
    targets_lower::SVector{3,Float32})
    @inbounds begin
        # Check all resource constraints
        for m in 1:3
            new_value = values_matrix[facility, m] + V[m][client]
            if new_value > targets_upper[m] || new_value < targets_lower[m]
                return false
            end
        end

        # Check risk constraint
        new_risk = risk_vec[facility] + R[client]
        return new_risk <= β
    end
end

# Check if removing a set of clients from their facilities maintains lower target constraints
function can_remove_clients(values_matrix::Matrix{Int}, V::Vector{Vector{Int}},
    clients_to_eject::Vector{Int}, facility_clients::Dict{Int,Vector{Int}},
    targets_lower::SVector{3,Float32})
    # Create temporary values matrix
    temp_values = copy(values_matrix)

    # Group clients by their current facility
    facility_removals = Dict{Int,Vector{Int}}()
    for client in clients_to_eject
        for (facility, clients) in facility_clients
            if client in clients
                if haskey(facility_removals, facility)
                    push!(facility_removals[facility], client)
                else
                    facility_removals[facility] = [client]
                end
                break
            end
        end
    end

    # Check each facility's lower bounds after removals
    @inbounds for (facility, clients) in facility_removals
        for m in 1:3
            # Remove all clients' contributions
            for client in clients
                temp_values[facility, m] -= V[m][client]
            end
            # Check if still above lower bound
            if temp_values[facility, m] < targets_lower[m]
                return false
            end
        end
    end
    return true
end

function can_assign_to_facility(client::Int, facility::Int, values_matrix::Matrix{Int},
    V::Vector{Vector{Int}}, risk_vec::Vector{Int}, R::Vector{Int},
    β::Int, targets_upper::SVector{3,Float32})
    @inbounds begin
        # Check all resource constraints
        for m in 1:3
            new_value = values_matrix[facility, m] + V[m][client]
            if new_value > targets_upper[m]
                return false
            end
        end

        # Check risk constraint
        new_risk = risk_vec[facility] + R[client]
        return new_risk <= β
    end
end

function find_best_ejection_chain_optimized(X::Matrix{Int}, D::Matrix{Int},
    V::Vector{Vector{Int}}, R::Vector{Int},
    B::Int, M::Int, usables_i::Vector{Int},
    values_matrix::Matrix{Int},
    risk_vec::Vector{Int},
    targets_lower::SVector{3,Float32},
    targets_upper::SVector{3,Float32},
    β::Int, chain_length::Int, strategy::Symbol)
    best_chain = nothing
    best_total_diff = 0.0

    # Create mapping of facilities to their currently assigned clients
    facility_clients = Dict{Int,Vector{Int}}()
    for i in usables_i
        facility_clients[i] = findall(j -> X[i, j] == 1, 1:B)
    end

    # Temporary matrices for evaluating moves
    temp_values = copy(values_matrix)
    temp_risks = copy(risk_vec)

    # Track ejected clients and their original facilities
    ejected = Dict{Int,Int}()  # client => original_facility
    moves = Vector{EjectionMove}()

    function try_reassign_clients(current_weight::Float64)
        if length(moves) == length(ejected)
            if strategy == :ff && current_weight < 0
                best_chain = copy(moves)
                return true
            elseif strategy == :bf && current_weight < best_total_diff
                best_chain = copy(moves)
                best_total_diff = current_weight
            end
            return false
        end

        # Try each remaining ejected client
        for (client, orig_facility) in ejected
            client in [m.client for m in moves] && continue

            # Try each possible facility
            for new_facility in usables_i
                new_facility == orig_facility && continue

                # Calculate weight difference
                weight_diff = D[new_facility, client] - D[orig_facility, client]
                new_total_weight = current_weight + weight_diff

                # Early pruning if not improving
                if strategy == :ff && new_total_weight >= 0
                    continue
                end

                # Check if assignment is feasible
                if can_assign_to_facility(client, new_facility, temp_values, V, temp_risks,
                    R, β, targets_upper)
                    # Temporarily apply the assignment
                    for m in 1:3
                        temp_values[new_facility, m] += V[m][client]
                    end
                    temp_risks[new_facility] += R[client]

                    # Record the move
                    push!(moves, EjectionMove(client, orig_facility, new_facility, weight_diff))

                    # Recursively try to assign remaining clients
                    if try_reassign_clients(new_total_weight)
                        return true
                    end

                    # Backtrack
                    pop!(moves)
                    for m in 1:3
                        temp_values[new_facility, m] -= V[m][client]
                    end
                    temp_risks[new_facility] -= R[client]
                end
            end
        end
        return false
    end

    # Try different sets of clients to eject
    for eject_count in 1:chain_length
        for i in 1:B
            clients_subset = collect(max(1, i - eject_count + 1):i)

            # First check if we can remove these clients without violating lower bounds
            if !can_remove_clients(values_matrix, V, clients_subset, facility_clients, targets_lower)
                continue
            end

            # Reset temporary state
            copyto!(temp_values, values_matrix)
            copyto!(temp_risks, risk_vec)
            empty!(ejected)
            empty!(moves)

            # Record original facilities and remove clients
            for client in clients_subset
                orig_facility = find_one_in_column_unrolled(X, client)
                ejected[client] = orig_facility

                # Remove from temporary matrices
                for m in 1:3
                    temp_values[orig_facility, m] -= V[m][client]
                end
                temp_risks[orig_facility] -= R[client]
            end

            # Try to find good reassignment of ejected clients
            try_reassign_clients(0.0)

            # If using first-fit and found improving solution, stop
            if strategy == :ff && best_chain !== nothing
                break
            end
        end
    end

    return best_chain
end
function ejection_chain_improve_optimized(solution, targets_lower::SVector{3,Float32},
    targets_upper::SVector{3,Float32}, strategy::Symbol,
    chain_length::Int)
    X, Y = solution.X, solution.Y
    instance = solution.Instance
    B, S, M = instance.B, instance.S, instance.M
    V, R, β = instance.V, instance.R, instance.β[1]
    D = instance.D
    Weight = solution.Weight

    usables_i = findall(==(1), Y)
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)

    while true
        best_moves = find_best_ejection_chain_optimized(X, D, V, R, B, M, usables_i,
            values_matrix, risk_vec,
            targets_lower, targets_upper,
            β, chain_length, strategy)
        best_moves === nothing && break

        # Apply the moves
        total_weight_diff = sum(move.weight_diff for move in best_moves)
        apply_ejection_chain!(X, values_matrix, risk_vec, best_moves, V, R)
        Weight += total_weight_diff
    end

    Weight = dot(X, D)  # Recalculate weight to ensure accuracy
    return Solution(instance, X, Y, Weight, solution.Time)
end

function apply_ejection_chain!(X::Matrix{Int}, values_matrix::Matrix{Int},
    risk_vec::Vector{Int}, moves::Vector{EjectionMove},
    V::Vector{Vector{Int}}, R::Vector{Int})
    # First verify the entire chain maintains feasibility
    temp_values = copy(values_matrix)
    temp_risks = copy(risk_vec)

    # Remove all clients from their original facilities
    @inbounds for move in moves
        X[move.from_facility, move.client] = 0

        for m in 1:3
            values_matrix[move.from_facility, m] -= V[m][move.client]
        end
        risk_vec[move.from_facility] -= R[move.client]
    end

    # Then assign them to their new facilities
    @inbounds for move in moves
        X[move.to_facility, move.client] = 1

        for m in 1:3
            values_matrix[move.to_facility, m] += V[m][move.client]
        end
        risk_vec[move.to_facility] += R[move.client]
    end
end

# Generate combinations of k elements from collection 1:n
function combinations_optimized(n::Int, k::Int)
    result = Vector{Vector{Int}}()
    temp = Vector{Int}(undef, k)

    function generate_combinations(start::Int, pos::Int)
        if pos > k
            push!(result, copy(temp))
            return
        end

        # Try each possible element for this position
        @inbounds for i in start:n-k+pos
            temp[pos] = i
            generate_combinations(i + 1, pos + 1)
        end
    end

    generate_combinations(1, 1)
    return result
end

# Alternative iterative version if memory is a concern
struct CombinationIterator
    n::Int
    k::Int
end

function Base.iterate(iter::CombinationIterator, state=nothing)
    n, k = iter.n, iter.k

    # Initialize first combination
    if state === nothing
        k == 0 && return (Int[], nothing)
        k > n && return nothing
        return (collect(1:k), collect(1:k))
    end

    current = state

    # Find the rightmost element that can be incremented
    i = k
    while i > 0 && current[i] == n - k + i
        i -= 1
    end

    # If no element can be incremented, we're done
    i == 0 && return nothing

    # Increment element i and set subsequent elements accordingly
    current[i] += 1
    for j in (i+1):k
        current[j] = current[j-1] + 1
    end

    return (copy(current), current)
end

Base.length(iter::CombinationIterator) = binomial(iter.n, iter.k)
Base.eltype(::Type{CombinationIterator}) = Vector{Int}

# New types to support proper ejection chain structure
struct DisplacementMove
    client::Int
    from_facility::Int
    to_facility::Int
    displaced_client::Union{Nothing,Int}
    weight_diff::Float64
end

struct ChainStep
    move::DisplacementMove
    values_delta::Vector{Float32}  # Changes in resource values
    risk_delta::Float32           # Change in risk value
end

mutable struct EjectionChain
    steps::Vector{ChainStep}
    total_weight_diff::Float64
end

function find_best_classical_chain(X::Matrix{Int}, D::Matrix{Int},
    V::Vector{Vector{Int}}, R::Vector{Int},
    values_matrix::Matrix{Int}, risk_vec::Vector{Int},
    targets_lower::SVector{3,Float32}, targets_upper::SVector{3,Float32},
    β::Int, max_chain_length::Int, Y::Vector{Int})

    best_chain = nothing
    best_improvement = 0.0

    # Try each client as the initial move
    for initial_client in 1:size(X, 2)
        current_facility = find_one_in_column_unrolled(X, initial_client)

        # Try moving to each possible facility
        for new_facility in findall(==(1), Y)
            new_facility == current_facility && continue

            chain = try_build_chain(
                initial_client, current_facility, new_facility,
                X, D, V, R, values_matrix, risk_vec,
                targets_lower, targets_upper, β, max_chain_length, Y
            )

            if chain !== nothing && chain.total_weight_diff < best_improvement
                best_chain = chain
                best_improvement = chain.total_weight_diff
            end
        end
    end

    return best_chain
end

function try_build_chain(initial_client::Int, from_facility::Int, to_facility::Int,
    X::Matrix{Int}, D::Matrix{Int}, V::Vector{Vector{Int}}, R::Vector{Int},
    values_matrix::Matrix{Int}, risk_vec::Vector{Int},
    targets_lower::SVector{3,Float32}, targets_upper::SVector{3,Float32},
    β::Int, max_chain_length::Int, Y)

    # Initialize chain
    chain = EjectionChain(ChainStep[], 0.0)

    # Temporary matrices for evaluating moves
    temp_values = copy(values_matrix)
    temp_risks = copy(risk_vec)
    temp_X = copy(X)

    # Try to build chain starting with initial move
    current_client = initial_client
    current_from = from_facility
    current_to = to_facility

    while length(chain.steps) < max_chain_length
        # Create potential move
        weight_diff = D[current_to, current_client] - D[current_from, current_client]

        # Find who we would displace (if anyone)
        displaced_client = find_displaced_client_enhanced(
            temp_X, current_to, current_client,
            D, V, R, temp_values, temp_risks,
            targets_upper, targets_lower, β, Y
        )
        move = DisplacementMove(
            current_client,
            current_from,
            current_to,
            displaced_client,
            weight_diff
        )

        # Calculate resource changes
        values_delta = zeros(Float32, 3)
        for m in 1:3
            values_delta[m] = V[m][current_client]
        end
        risk_delta = R[current_client]

        # Check if move is feasible
        if !is_move_feasible(move, temp_values, temp_risks, values_delta, risk_delta,
            targets_lower, targets_upper, β)
            return nothing
        end

        # Apply move temporarily
        apply_move!(temp_X, temp_values, temp_risks, move, values_delta, risk_delta)

        # Record step
        push!(chain.steps, ChainStep(move, values_delta, risk_delta))
        chain.total_weight_diff += weight_diff

        # If no displacement, chain ends naturally
        if displaced_client === nothing
            return chain
        end

        # Set up next iteration with displaced client
        current_client = displaced_client
        current_from = current_to

        # Find best facility for displaced client
        best_facility = nothing
        best_facility_cost = Inf

        for candidate_facility in findall(==(1), Y)
            candidate_facility == current_from && continue

            if can_assign_to_facility(current_client, candidate_facility, temp_values, V, temp_risks,
                R, β, targets_upper, targets_lower)
                facility_cost = D[candidate_facility, current_client]
                if facility_cost < best_facility_cost
                    best_facility = candidate_facility
                    best_facility_cost = facility_cost
                end
            end
        end

        # If we can't place displaced client, chain fails
        if best_facility === nothing
            return nothing
        end

        current_to = best_facility
    end

    return chain
end

function is_move_feasible(move::DisplacementMove,
    values_matrix::Matrix{Int}, risk_vec::Vector{Int},
    values_delta::Vector{Float32}, risk_delta::Int64,
    targets_lower::SVector{3,Float32}, targets_upper::SVector{3,Float32}, β::Int)

    # Check resource constraints at both facilities
    for m in 1:3
        # Check removal from original facility
        new_value_from = values_matrix[move.from_facility, m] - values_delta[m]
        if new_value_from < targets_lower[m]
            return false
        end

        # Check addition to new facility
        new_value_to = values_matrix[move.to_facility, m] + values_delta[m]
        if new_value_to > targets_upper[m]
            return false
        end
    end

    # Check risk constraints
    new_risk_from = risk_vec[move.from_facility] - risk_delta
    new_risk_to = risk_vec[move.to_facility] + risk_delta

    return new_risk_to <= β
end

function apply_move!(X::Matrix{Int},
    values_matrix::Matrix{Int}, risk_vec::Vector{Int},
    move::DisplacementMove, values_delta::Vector{Float32}, risk_delta::Int64)

    # Update assignment matrix
    X[move.from_facility, move.client] = 0
    X[move.to_facility, move.client] = 1

    # Update resource values
    for m in 1:3
        values_matrix[move.from_facility, m] -= values_delta[m]
        values_matrix[move.to_facility, m] += values_delta[m]
    end

    # Update risk values
    risk_vec[move.from_facility] -= risk_delta
    risk_vec[move.to_facility] += risk_delta
end

# Replace your current ejection chain improve with this
function classical_ejection_chain_improve(solution, targets_lower::SVector{3,Float32},
    targets_upper::SVector{3,Float32}, max_chain_length::Int)

    X, Y = solution.X, solution.Y
    instance = solution.Instance
    B, S, M = instance.B, instance.S, instance.M
    V, R, β = instance.V, instance.R, instance.β[1]
    D = instance.D

    # Initialize resource tracking matrices
    values_matrix = Matrix{Float32}(undef, S, M)
    risk_vec = Vector{Float32}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)

    improved = true
    while improved
        improved = false

        best_chain = find_best_classical_chain(
            X, D, V, R, values_matrix, risk_vec,
            targets_lower, targets_upper, β, max_chain_length, Y
        )

        if best_chain !== nothing
            # Apply the chain
            for step in best_chain.steps
                apply_move!(X, values_matrix, risk_vec, step.move, step.values_delta, step.risk_delta)
            end
            improved = true
        end
    end

    Weight = dot(X, D)  # Recalculate final weight
    return Solution(instance, X, Y, Weight, solution.Time)
end

# Utility functions for finding current assignments and displaced clients

function find_displaced_client(X::Matrix{Int}, facility::Int, incoming_client::Int)
    """
    Find if any client will be displaced from the facility by the incoming client.
    Returns the displaced client or nothing if no displacement occurs.
    
    Parameters:
    - X: Current assignment matrix
    - facility: Target facility
    - incoming_client: Client being moved to the facility
    """
    # Get current assignments at the facility
    current_clients = findall(j -> X[facility, j] == 1, 1:size(X, 2))

    if isempty(current_clients)
        return nothing
    end

    # In your case, since there might be capacity/resource constraints,
    # we might need to displace a client even if there's theoretical room

    # Strategy 1: Displace the client with highest cost if needed
    if !isempty(current_clients)
        # You might want to choose based on various criteria:
        # - Highest cost client
        # - Client that violates constraints least when moved
        # - Client with most alternative placement options
        return current_clients[1]  # For now, just take the first one
    end

    return nothing
end

function find_current_facility(X::Matrix{Int}, client::Int)
    """
    Find the facility currently serving a client.
    
    Parameters:
    - X: Assignment matrix
    - client: Client to find
    """
    for i in 1:size(X, 1)
        if X[i, client] == 1
            return i
        end
    end
    error("Client $client not assigned to any facility")
end

function find_best_facility_for_displaced(client::Int, current_facility::Int,
    X::Matrix{Int}, D::Matrix{Int},
    values_matrix::Matrix{Int}, risk_vec::Vector{Int},
    V::Vector{Vector{Int}}, R::Vector{Int},
    targets_upper::SVector{3,Float32}, targets_lower::SVector{3,Float32},
    β::Int, Y::Vector{Int})
    """
    Find the best facility to place a displaced client.
    
    Returns (facility_index, cost) or (nothing, Inf) if no feasible placement found.
    """
    best_facility = nothing
    best_cost = Inf

    for facility in findall(==(1), Y)
        facility == current_facility && continue

        # Check if assignment is feasible
        if can_assign_to_facility(client, facility, values_matrix, V, risk_vec,
            R, β, targets_upper, targets_lower)
            cost = D[facility, client]
            if cost < best_cost
                best_facility = facility
                best_cost = cost
            end
        end
    end

    return best_facility, best_cost
end

# Enhanced version of find_displaced_client with more sophisticated selection
function find_displaced_client_enhanced(X::Matrix{Int}, facility::Int, incoming_client::Int,
    D::Matrix{Int}, V::Vector{Vector{Int}}, R::Vector{Int},
    values_matrix::Matrix{Int}, risk_vec::Vector{Int},
    targets_upper::SVector{3,Float32}, targets_lower::SVector{3,Float32},
    β::Int, Y::Vector{Int})
    """
    Enhanced version that considers multiple factors when choosing which client to displace.
    """
    current_clients = findall(j -> X[facility, j] == 1, 1:size(X, 2))
    isempty(current_clients) && return nothing

    # Calculate scores for each potential client to displace
    scores = Dict{Int,Float64}()
    for client in current_clients
        # 1. Cost of current assignment
        current_cost = D[facility, client]

        # 2. Number of alternative facilities available
        alternatives = count(i ->
                i != facility && Y[i] == 1 &&
                    can_assign_to_facility(client, i, values_matrix, V, risk_vec,
                        R, β, targets_upper, targets_lower),
            1:size(X, 1)
        )

        # 3. Resource usage
        resource_usage = sum(V[m][client] for m in 1:length(V))

        # 4. Risk contribution
        risk_contribution = R[client]

        # Combine factors into a score (lower is better for displacement)
        scores[client] = -alternatives * 100 +  # Weight heavily towards clients with more options
                         current_cost * 0.5 +    # Consider current cost
                         resource_usage * 0.3 +  # Consider resource usage
                         risk_contribution * 0.2 # Consider risk contribution
    end

    # Return client with lowest score (best to displace)
    return argmin(scores)
end

# Function to check if a move maintains feasibility
function move_maintains_feasibility(client::Int, from_facility::Int, to_facility::Int,
    X::Matrix{Int}, V::Vector{Vector{Int}}, R::Vector{Int},
    values_matrix::Matrix{Int}, risk_vec::Vector{Int},
    targets_upper::SVector{3,Float32}, targets_lower::SVector{3,Float32},
    β::Int)
    """
    Check if moving a client from one facility to another maintains feasibility
    """
    # Create temporary matrices to check the move
    temp_values = copy(values_matrix)
    temp_risks = copy(risk_vec)

    # Remove from original facility
    for m in 1:length(V)
        temp_values[from_facility, m] -= V[m][client]
    end
    temp_risks[from_facility] -= R[client]

    # Check if removal maintains lower bounds
    for m in 1:length(V)
        if temp_values[from_facility, m] < targets_lower[m]
            return false
        end
    end

    # Add to new facility
    for m in 1:length(V)
        temp_values[to_facility, m] += V[m][client]
    end
    temp_risks[to_facility] += R[client]

    # Check if addition maintains upper bounds and risk constraint
    for m in 1:length(V)
        if temp_values[to_facility, m] > targets_upper[m]
            return false
        end
    end

    return temp_risks[to_facility] <= β
end

function repair_solution_improved(solution, cons, targets_lower, targets_upper)
    X = copy(solution.X)
    Y = copy(solution.Y)
    instance = solution.Instance
    B = instance.B
    S = instance.S
    M = instance.M
    V = instance.V
    R = instance.R
    β = instance.β
    D = instance.D

    # Initialize constraint tracking matrices
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)

    # For each facility, identify activities that are over/under
    facility_issues = Dict{Int,Dict{String,Vector{Int}}}()
    for i in 1:S
        over_activities = Int[]
        under_activities = Int[]
        for m in 1:M
            if values_matrix[i, m] > targets_upper[i, m]
                push!(over_activities, m)
            elseif values_matrix[i, m] < targets_lower[i, m]
                push!(under_activities, m)
            end
        end

        if !isempty(over_activities) || !isempty(under_activities)
            facility_issues[i] = Dict(
                "over" => over_activities,
                "under" => under_activities
            )
        end
    end

    # Try direct moves first
    for (i, issues) in facility_issues
        assigned_buses = findall(==(1), X[i, :])

        # If we have over-activities, try to remove buses
        if !isempty(issues["over"])
            bus_scores = Tuple{Int,Float64}[]
            for j in assigned_buses
                score = 0.0
                # Score based on contribution to over/under activities
                for m in issues["over"]
                    score += V[m][j]  # Higher contribution to over = better to remove
                end
                for m in get(issues, "under", Int[])  # Use get for safe access
                    score -= 2 * V[m][j]  # Contribution to under = bad to remove
                end
                push!(bus_scores, (j, score))
            end

            sort!(bus_scores, by=x -> x[2], rev=true)

            # Try to reassign highest scoring buses
            for (j, _) in bus_scores
                # Find potential new facilities
                potential_facilities = findall(==(1), Y)
                filter!(ĩ -> ĩ != i, potential_facilities)

                # Score each potential facility
                facility_scores = Tuple{Int,Float64}[]
                for ĩ in potential_facilities
                    score = -D[ĩ, j]  # Base score on negative distance (closer is better)
                    feasible = true

                    # Check feasibility and calculate margin to constraints
                    for m in 1:M
                        new_value = values_matrix[ĩ, m] + V[m][j]
                        if new_value > targets_upper[ĩ, m]
                            feasible = false
                            break
                        end
                        # Add score for how far from upper bound
                        score += (targets_upper[ĩ, m] - new_value)
                    end

                    # Check risk constraint
                    if risk_vec[ĩ] + R[j] > β[ĩ]
                        feasible = false
                    end

                    if feasible
                        push!(facility_scores, (ĩ, score))
                    end
                end

                # Try best scoring facility
                if !isempty(facility_scores)
                    sort!(facility_scores, by=x -> x[2], rev=true)
                    ĩ = facility_scores[1][1]

                    # Make the move
                    X[i, j] = 0
                    X[ĩ, j] = 1

                    # Update tracking
                    for m in 1:M
                        values_matrix[i, m] -= V[m][j]
                        values_matrix[ĩ, m] += V[m][j]
                    end
                    risk_vec[i] -= R[j]
                    risk_vec[ĩ] += R[j]
                end
            end
        end
    end

    # Try swaps for remaining issues
    for (i, issues) in facility_issues
        if !isempty(issues["over"]) && !isempty(issues["under"])
            assigned_buses = findall(==(1), X[i, :])

            for j in assigned_buses
                best_swap = nothing
                best_improvement = -Inf

                # Look for beneficial swaps with other facilities
                for ĩ in 1:S
                    if ĩ == i || Y[ĩ] == 0
                        continue
                    end

                    for j̃ in findall(==(1), X[ĩ, :])
                        # Calculate potential improvement
                        feasible = true
                        improvement = 0.0

                        # Simulate swap
                        temp_values_i = copy(values_matrix[i, :])
                        temp_values_ĩ = copy(values_matrix[ĩ, :])
                        temp_risk_i = risk_vec[i]
                        temp_risk_ĩ = risk_vec[ĩ]

                        # Update temp values
                        for m in 1:M
                            temp_values_i[m] -= V[m][j]
                            temp_values_i[m] += V[m][j̃]
                            temp_values_ĩ[m] -= V[m][j̃]
                            temp_values_ĩ[m] += V[m][j]

                            # Check constraints
                            if temp_values_i[m] > targets_upper[i, m] ||
                               temp_values_i[m] < targets_lower[i, m] ||
                               temp_values_ĩ[m] > targets_upper[ĩ, m] ||
                               temp_values_ĩ[m] < targets_lower[ĩ, m]
                                feasible = false
                                break
                            end

                            # Calculate improvement
                            if m in issues["over"]
                                improvement += values_matrix[i, m] - temp_values_i[m]
                            end
                            if m in issues["under"]
                                improvement += temp_values_i[m] - values_matrix[i, m]
                            end
                        end

                        # Check risk constraints
                        temp_risk_i = risk_vec[i] - R[j] + R[j̃]
                        temp_risk_ĩ = risk_vec[ĩ] - R[j̃] + R[j]
                        if temp_risk_i > β[i] || temp_risk_ĩ > β[ĩ]
                            feasible = false
                        end

                        # Consider distance in improvement
                        distance_change = D[i, j̃] + D[ĩ, j] - D[i, j] - D[ĩ, j̃]
                        improvement -= distance_change

                        if feasible && improvement > best_improvement
                            best_improvement = improvement
                            best_swap = (ĩ, j̃)
                        end
                    end
                end

                # Make the best swap if found
                if best_swap !== nothing
                    ĩ, j̃ = best_swap
                    # Swap buses
                    X[i, j] = 0
                    X[i, j̃] = 1
                    X[ĩ, j̃] = 0
                    X[ĩ, j] = 1

                    # Update tracking
                    for m in 1:M
                        values_matrix[i, m] = values_matrix[i, m] - V[m][j] + V[m][j̃]
                        values_matrix[ĩ, m] = values_matrix[ĩ, m] - V[m][j̃] + V[m][j]
                    end
                    risk_vec[i] = risk_vec[i] - R[j] + R[j̃]
                    risk_vec[ĩ] = risk_vec[ĩ] - R[j̃] + R[j]
                end
            end
        end
    end

    # Calculate new solution weight
    Weight = sum(D[i, j] for i in 1:S, j in 1:B if X[i, j] == 1)
    sol = Solution(instance, X, Y, Weight, solution.Time)
    println("===========================================================")
    isFactible(sol, true)
    return sol
end


struct ActivityMetric
    avg::Float64
    min::Float64
    max::Float64
    range::Float64
end

mutable struct FacilityState
    id::Int
    current_values::Vector{Float64}
    min_thresholds::Vector{Float64}
    max_thresholds::Vector{Float64}
    normalized_gaps::Vector{Float64}
    current_risk::Float64
    risk_capacity::Float64
end

function get_activity_metrics(instance::Instance)
    activity_metrics = Vector{ActivityMetric}(undef, instance.M)

    for m in 1:instance.M
        values = Float64.(instance.V[m])  # Convert to Float64 for precision
        activity_metrics[m] = ActivityMetric(
            mean(values),
            minimum(values),
            maximum(values),
            maximum(values) - minimum(values)
        )
    end

    return activity_metrics
end

function repair_solution_v2(solution, targets_lower, targets_upper)
    instance = solution.Instance
    X = copy(solution.X)
    Y = copy(solution.Y)

    # Get activity metrics for normalization
    activity_metrics = get_activity_metrics(instance)

    # Initialize tracking
    values_matrix = zeros(Int, instance.S, instance.M)
    risk_vec = zeros(Int, instance.S)
    for i in 1:instance.S, j in 1:instance.B
        if X[i, j] == 1
            for m in 1:instance.M
                values_matrix[i, m] += instance.V[m][j]
            end
            risk_vec[i] += instance.R[j]
        end
    end

    function get_violation_severity(facility)
        severity = 0.0
        for m in 1:instance.M
            current = values_matrix[facility, m]
            under = max(0, targets_lower[facility, m] - current)
            over = max(0, current - targets_upper[facility, m])
            # Normalize by activity scale
            severity += under / activity_metrics[m].avg + over / activity_metrics[m].avg
        end

        # Add risk violation if any
        if risk_vec[facility] > instance.β[facility]
            severity += (risk_vec[facility] - instance.β[facility]) / instance.β[facility]
        end

        return severity
    end

    function simulate_swap(i, j, ĩ, j̃)
        # Calculate new states after potential swap
        new_values_i = copy(values_matrix[i, :])
        new_values_ĩ = copy(values_matrix[ĩ, :])
        new_risk_i = risk_vec[i] - instance.R[j] + instance.R[j̃]
        new_risk_ĩ = risk_vec[ĩ] - instance.R[j̃] + instance.R[j]

        for m in 1:instance.M
            new_values_i[m] = new_values_i[m] - instance.V[m][j] + instance.V[m][j̃]
            new_values_ĩ[m] = new_values_ĩ[m] - instance.V[m][j̃] + instance.V[m][j]
        end

        return new_values_i, new_values_ĩ, new_risk_i, new_risk_ĩ
    end

    function evaluate_state(values, risk, facility, tolerance)
        violation = 0.0

        for m in 1:instance.M
            lower = targets_lower[facility, m] * (1 - tolerance)
            upper = targets_upper[facility, m] * (1 + tolerance)

            under = max(0, lower - values[m])
            over = max(0, values[m] - upper)
            violation += under / activity_metrics[m].avg + over / activity_metrics[m].avg
        end

        if risk > instance.β[facility] * (1 + tolerance)
            violation += (risk - instance.β[facility]) / instance.β[facility]
        end

        return violation
    end

    # Main repair loop
    tolerance = 0.05  # Initial 5% violation allowance
    max_iterations = 100
    iteration = 0

    while tolerance > 0.001 && iteration < max_iterations
        println("Iteration $iteration, tolerance: $tolerance")
        made_improvement = false

        # Get and sort facilities by violation severity
        facility_violations = [(f=i, v=get_violation_severity(i)) for i in 1:instance.S if Y[i] == 1]
        sort!(facility_violations, by=x -> x.v, rev=true)

        # Try improvements for each violating facility
        for (i, _) in facility_violations
            current_violation_i = get_violation_severity(i)
            current_violation_i == 0 && continue

            assigned_to_i = findall(==(1), X[i, :])

            # Try swaps with other facilities
            for ĩ in (f.f for f in facility_violations if f.f != i)
                assigned_to_ĩ = findall(==(1), X[ĩ, :])

                for j in assigned_to_i, j̃ in assigned_to_ĩ
                    # Simulate swap
                    new_values_i, new_values_ĩ, new_risk_i, new_risk_ĩ =
                        simulate_swap(i, j, ĩ, j̃)

                    # Evaluate new states with current tolerance
                    new_violation_i = evaluate_state(new_values_i, new_risk_i, i, tolerance)
                    new_violation_ĩ = evaluate_state(new_values_ĩ, new_risk_ĩ, ĩ, tolerance)

                    current_violation_ĩ = get_violation_severity(ĩ)
                    total_current = current_violation_i + current_violation_ĩ
                    total_new = new_violation_i + new_violation_ĩ

                    # Accept if total violation decreases
                    if total_new < total_current
                        # Execute swap
                        X[i, j], X[i, j̃] = 0, 1
                        X[ĩ, j̃], X[ĩ, j] = 0, 1

                        # Update tracking
                        values_matrix[i, :], values_matrix[ĩ, :] = new_values_i, new_values_ĩ
                        risk_vec[i], risk_vec[ĩ] = new_risk_i, new_risk_ĩ

                        made_improvement = true
                        println("Swap improved violation from $total_current to $total_new")
                        break
                    end
                end
                made_improvement && break
            end
            made_improvement && break
        end

        if !made_improvement
            tolerance *= 0.5
            println("No improvement found, reducing tolerance to $tolerance")
        end

        iteration += 1
    end

    # Calculate final weight
    weight = sum(instance.D[i, j] for i in 1:instance.S, j in 1:instance.B if X[i, j] == 1)
    return Solution(instance, X, Y, weight, solution.Time)
end

function repair_solution_v3(solution, targets_lower, targets_upper)
    instance = solution.Instance
    X = copy(solution.X)
    Y = copy(solution.Y)

    # Initialize tracking matrices
    values_matrix = zeros(Float32, instance.S, instance.M)
    risk_vec = zeros(Float32, instance.S)

    # Initial state
    for i in 1:instance.S, j in 1:instance.B
        if X[i, j] == 1
            for m in 1:instance.M
                values_matrix[i, m] += instance.V[m][j]
            end
            risk_vec[i] += instance.R[j]
        end
    end

    function calculate_violation(values_matrix, risk_vec)
        total = 0.0
        for i in 1:instance.S
            if Y[i] == 0
                continue
            end

            for m in 1:instance.M
                target = (targets_upper[i, m] + targets_lower[i, m]) / 2
                current = values_matrix[i, m]

                # Relative violation calculation
                if current < targets_lower[i, m]
                    relative_gap = (targets_lower[i, m] - current) / target
                    total += relative_gap^2  # Square to emphasize larger violations
                elseif current > targets_upper[i, m]
                    relative_gap = (current - targets_upper[i, m]) / target
                    total += relative_gap^2
                end
            end
        end
        return total
    end

    function identify_facilities_to_repair()
        needs_repair = Dict{Int,Vector{Symbol}}()

        for i in 1:instance.S
            Y[i] == 0 && continue

            needs = Symbol[]
            for m in 1:instance.M
                if values_matrix[i, m] > targets_upper[i, m]
                    push!(needs, :remove)
                    break
                end
            end

            for m in 1:instance.M
                if values_matrix[i, m] < targets_lower[i, m]
                    push!(needs, :add)
                    break
                end
            end

            if risk_vec[i] > instance.β[i]
                push!(needs, :remove)
            end

            !isempty(needs) && (needs_repair[i] = unique(needs))
        end

        return needs_repair
    end

    function simulate_move(i, j, ĩ)
        new_values = copy(values_matrix)
        new_risks = copy(risk_vec)

        # Remove from i, add to ĩ
        for m in 1:instance.M
            new_values[i, m] -= instance.V[m][j]
            new_values[ĩ, m] += instance.V[m][j]
        end
        new_risks[i] -= instance.R[j]
        new_risks[ĩ] += instance.R[j]

        return new_values, new_risks
    end

    # Main repair loop
    max_iterations = 150
    iteration = 0
    best_violation = calculate_violation(values_matrix, risk_vec)
    println("Initial violation: $best_violation")

    # Track moves to avoid cycles
    tabu_moves = Set{Tuple{Int,Int,Int}}()  # (from_facility, to_facility, bu)
    tabu_tenure = 5

    while iteration < max_iterations && best_violation > 0.001
        needs_repair = identify_facilities_to_repair()
        isempty(needs_repair) && break

        made_improvement = false
        current_violation = calculate_violation(values_matrix, risk_vec)

        # Handle facilities needing repair
        for (facility, needs) in needs_repair
            if :remove in needs
                # Find candidate BUs to remove
                assigned_bus = findall(==(1), X[facility, :])
                # Sort by distance, furthest first
                sorted_bus = sort(assigned_bus,
                    by=j -> instance.D[facility, j],
                    rev=true)

                for j in sorted_bus
                    # Find candidate facilities to receive this BU
                    candidates = [i for i in 1:instance.S
                                  if Y[i] == 1 && i != facility]

                    # Sort by distance to BU
                    sort!(candidates, by=i -> instance.D[i, j])

                    for ĩ in candidates
                        # Skip if move is in tabu list
                        (facility, ĩ, j) in tabu_moves && continue

                        # Simulate move
                        new_values, new_risks = simulate_move(facility, j, ĩ)
                        new_violation = calculate_violation(new_values, new_risks)

                        # Accept if improves violation
                        if new_violation < current_violation
                            # Execute move
                            X[facility, j] = 0
                            X[ĩ, j] = 1
                            values_matrix = new_values
                            risk_vec = new_risks

                            # Update tracking
                            push!(tabu_moves, (facility, ĩ, j))
                            length(tabu_moves) > tabu_tenure && pop!(tabu_moves)

                            made_improvement = true
                            current_violation = new_violation
                            if current_violation < best_violation
                                best_violation = current_violation
                                println("Iteration $iteration: New best violation: $best_violation")
                            end
                            break
                        end
                    end
                    made_improvement && break
                end
            end

            if :add in needs && !made_improvement
                # Similar logic for facilities needing additions
                # Look for closest unassigned or poorly assigned BUs
                potential_bus = findall(j -> sum(X[:, j]) == 1, 1:instance.B)
                sort!(potential_bus, by=j -> instance.D[facility, j])

                for j in potential_bus
                    current_facility = findfirst(i -> X[i, j] == 1, 1:instance.S)
                    (current_facility, facility, j) in tabu_moves && continue

                    new_values, new_risks = simulate_move(current_facility, j, facility)
                    new_violation = calculate_violation(new_values, new_risks)

                    if new_violation < current_violation
                        X[current_facility, j] = 0
                        X[facility, j] = 1
                        values_matrix = new_values
                        risk_vec = new_risks

                        push!(tabu_moves, (current_facility, facility, j))
                        length(tabu_moves) > tabu_tenure && pop!(tabu_moves)

                        made_improvement = true
                        current_violation = new_violation
                        if current_violation < best_violation
                            best_violation = current_violation
                            println("Iteration $iteration: New best violation: $best_violation")
                        end
                        break
                    end
                end
            end

            made_improvement && break
        end

        !made_improvement && println("Iteration $iteration: No improvement found")
        iteration += 1
    end

    # Calculate final weight
    weight = sum(instance.D[i, j] for i in 1:instance.S, j in 1:instance.B if X[i, j] == 1)


    println("Repair finished after $iteration iterations")

    println("SOLUTION")
    sol3 = Solution(instance, X, Y, weight, solution.Time)
    println(isFactible(sol3))
    println("Final violation: $(calculate_violation(values_matrix, risk_vec))")
    new_sol = Solution(instance, X, Y, weight, solution.Time)
    fine_tune_assignments!(X, values_matrix, risk_vec, targets_lower, targets_upper, instance, new_sol)
    weight = sum(instance.D[i, j] for i in 1:instance.S, j in 1:instance.B if X[i, j] == 1)
    return Solution(instance, X, Y, weight, solution.Time)
end

function fine_tune_assignments!(X, values_matrix, risk_vec, targets_lower, targets_upper, instance, solution)
    instance = solution.Instance
    X = copy(solution.X)
    Y = copy(solution.Y)

    function get_split_violations()
        split_facilities = Dict{Int, Dict{Int,Tuple{Float64,Float64}}}()
        
        for i in 1:instance.S
            Y[i] == 0 && continue
            
            activity_gaps = Dict{Int,Tuple{Float64,Float64}}()
            has_under = false
            has_over = false
            
            for m in 1:instance.M
                current = values_matrix[i,m]
                if current < targets_lower[i,m]
                    gap = targets_lower[i,m] - current
                    activity_gaps[m] = (gap, 0.0)  # (need_to_add, need_to_reduce)
                    has_under = true
                elseif current > targets_upper[i,m]
                    gap = current - targets_upper[i,m]
                    activity_gaps[m] = (0.0, gap)
                    has_over = true
                end
            end
            
            # Only include if facility has both under and over violations
            if has_under && has_over
                split_facilities[i] = activity_gaps
            end
        end
        return split_facilities
    end

    function evaluate_swap_for_split(i, j, ĩ, j̃, needed_gaps)
        improvement_score = 0.0
        total_improvements = 0  # Count how many violations we improve
        creates_new_violation = false
        
        for (m, (need_add, need_reduce)) in needed_gaps
            current_change = -instance.V[m][j] + instance.V[m][j̃]  # Net change for facility i
            
            if need_add > 0  # Need to increase this activity
                if current_change > 0
                    improvement_score += current_change  # Use absolute improvement
                    total_improvements += 1
                elseif current_change < 0  # Moving in wrong direction
                    improvement_score -= abs(current_change)  # Penalize
                    creates_new_violation = true
                end
            elseif need_reduce > 0  # Need to decrease this activity
                if current_change < 0
                    improvement_score += abs(current_change)  # Use absolute improvement
                    total_improvements += 1
                elseif current_change > 0  # Moving in wrong direction
                    improvement_score -= current_change  # Penalize
                    creates_new_violation = true
                end
            end
        end
        
        # Only consider valid if we improve at least as many violations as we have
        if total_improvements < length(needed_gaps) || creates_new_violation
            return -Inf
        end
        
        return improvement_score
     end

    function try_fix_split_violations(split_facilities)
        made_improvement = false
        
        for (i, needed_gaps) in split_facilities
            assigned_to_i = findall(==(1), X[i, :])
            best_swap = nothing
            best_score = -1.0
            
            for ĩ in 1:instance.S
                ĩ == i && continue
                Y[ĩ] == 0 && continue
                
                assigned_to_ĩ = findall(==(1), X[ĩ, :])
                
                for j in assigned_to_i
                    for j̃ in assigned_to_ĩ
                        # Check if swap would be feasible for facility ĩ
                        would_violate = false
                        new_values = copy(values_matrix)
                        new_risks = copy(risk_vec)
                        
                        # Simulate swap
                        for m in 1:instance.M
                            new_values[i,m] = values_matrix[i,m] - instance.V[m][j] + instance.V[m][j̃]
                            new_values[ĩ,m] = values_matrix[ĩ,m] - instance.V[m][j̃] + instance.V[m][j]
                            
                            # Check if swap would create new violations in ĩ
                            if new_values[ĩ,m] < targets_lower[ĩ,m] || 
                               new_values[ĩ,m] > targets_upper[ĩ,m]
                                would_violate = true
                                break
                            end
                        end
                        
                        if !would_violate
                            score = evaluate_swap_for_split(i, j, ĩ, j̃, needed_gaps)
                            if score > best_score
                                best_score = score
                                best_swap = (j, j̃, ĩ, new_values, new_risks)
                            end
                        end
                    end
                end
            end
            
            # Execute best swap if found
            if best_swap !== nothing
                j, j̃, ĩ, new_values, new_risks = best_swap
                X[i,j], X[i,j̃] = 0, 1
                X[ĩ,j̃], X[ĩ,j] = 0, 1
                values_matrix = new_values
                risk_vec = new_risks
                made_improvement = true
                println("Fixed split violation in facility $i with score $best_score")
            end
        end
        
        return made_improvement
    end

    # Main loop
    iteration = 0
    max_iterations = 150
    
    while iteration < max_iterations
        split_facilities = get_split_violations()
        if isempty(split_facilities)
            println("No split violations found")
            break
        end
        
        println("Iteration $iteration: Found $(length(split_facilities)) facilities with split violations")
        if !try_fix_split_violations(split_facilities)
            println("Could not improve split violations further")
            break
        end
        
        iteration += 1
    end
    
    return X
end

function swap_problematic_centers!(X, Y, values_matrix, risk_vec, targets_lower, targets_upper, instance)
    B, S, M, V, R, β, P = instance.B, instance.S, instance.M, instance.V, instance.R, instance.β, instance.P
    D = instance.D
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk

    function identify_problematic_centers()
        problematic = Set{Int}()
        for i in 1:S
            Y[i] == 0 && continue
            
            has_over = false
            has_under = false
            for m in 1:M
                if values_matrix[i,m] > targets_upper[i,m]
                    has_over = true
                elseif values_matrix[i,m] < targets_lower[i,m]
                    has_under = true
                end
            end
            
            has_over || has_under && push!(problematic, i)
        end
        println(problematic)
        return problematic
    end

    function rollback_changes(modified_X, modified_values, modified_risk)
        # Restore X matrix
        for ((i,j), val) in modified_X
            X[i,j] = val
        end
        
        # Restore values matrix
        for ((i,m), val) in modified_values
            values_matrix[i,m] = val
        end
        
        # Restore risk vector
        for (i, val) in modified_risk
            risk_vec[i] = val
        end
        
        return false  # Indicate the swap was unsuccessful
    end

    
    function try_center_swap(ĩ, i✶)
        modified_X = Dict{Tuple{Int,Int},Int}()
        modified_values = Dict{Tuple{Int,Int},Int}()
        modified_risk = Dict{Int,Int}()
        
        # First, get the BUs that will become orphaned from ĩ
        js_assigned_to_old = findall(==(1), X[ĩ,:])
        js_assigned_set = Set(js_assigned_to_old)
        
        # Clear the old center
        for j in js_assigned_to_old
            modified_X[(ĩ, j)] = X[ĩ, j]
            for m in 1:M
                modified_values[(ĩ, m)] = values_matrix[ĩ, m]
                values_matrix[ĩ, m] = 0
            end
            modified_risk[ĩ] = risk_vec[ĩ]
            risk_vec[ĩ] = 0
        end
        X[ĩ,:] .= 0
        
        # Try to populate new center with any suitable BUs
        fulls_m = zeros(Int, M)
        factible_yet = false
        
        # Try assigning BUs to new center (considering all BUs, not just the orphaned ones)
        for j in 1:B
            potential_assignment_valid = true
            current_center = findfirst(==(1), X[:,j])
            
            if current_center !== nothing  # If BU is assigned somewhere
                # Check if we can remove from current center
                for m in 1:M
                    if values_matrix[current_center, m] - V[m][j] < targets_lower[current_center, m]
                        potential_assignment_valid = false
                        break
                    end
                end
            end
            
            # Check if we can add to new center
            if potential_assignment_valid
                for m in 1:M
                    if values_matrix[i✶, m] + V[m][j] > targets_upper[i✶, m]
                        potential_assignment_valid = false
                        break
                    end
                end
                
                if risk_vec[i✶] + R[j] > β[i✶]
                    potential_assignment_valid = false
                end
            end
            
            if potential_assignment_valid
                # Make the assignment
                if current_center !== nothing
                    modified_X[(current_center, j)] = X[current_center, j]
                    X[current_center, j] = 0
                    for m in 1:M
                        modified_values[(current_center, m)] = values_matrix[current_center, m]
                        values_matrix[current_center, m] -= V[m][j]
                    end
                    modified_risk[current_center] = risk_vec[current_center]
                    risk_vec[current_center] -= R[j]
                end
                
                modified_X[(i✶, j)] = X[i✶, j]
                X[i✶, j] = 1
                
                for m in 1:M
                    modified_values[(i✶, m)] = get(modified_values, (i✶, m), values_matrix[i✶, m])
                    values_matrix[i✶, m] += V[m][j]
                    if values_matrix[i✶, m] > targets_lower[i✶, m]
                        fulls_m[m] = 1
                    end
                end
                
                modified_risk[i✶] = get(modified_risk, i✶, risk_vec[i✶])
                risk_vec[i✶] += R[j]
                
                if j in js_assigned_set
                    delete!(js_assigned_set, j)
                end
                
                if all(==(1), fulls_m)
                    factible_yet = true
                    break
                end
            end
        end
        
        if !factible_yet
            # Couldn't populate new center adequately
            return rollback_changes(modified_X, modified_values, modified_risk)
        end
        
        # Now handle orphaned BUs
        for j in js_assigned_set
            assigned = false
            for i in 1:S
                (i == ĩ || Y[i] == 0) && continue  # Skip closed centers and the one we're closing
                
                # Check if assignment is feasible
                can_assign = true
                for m in 1:M
                    if values_matrix[i,m] + V[m][j] > targets_upper[i,m]
                        can_assign = false
                        break
                    end
                end
                
                if risk_vec[i] + R[j] > β[i]
                    can_assign = false
                end
                
                if can_assign
                    modified_X[(i,j)] = X[i,j]
                    X[i,j] = 1
                    for m in 1:M
                        modified_values[(i,m)] = get(modified_values, (i,m), values_matrix[i,m])
                        values_matrix[i,m] += V[m][j]
                    end
                    modified_risk[i] = get(modified_risk, i, risk_vec[i])
                    risk_vec[i] += R[j]
                    assigned = true
                    break
                end
            end
            
            if !assigned
                return rollback_changes(modified_X, modified_values, modified_risk)
            end
        end
        
        return true
    end


    problematic = identify_problematic_centers()
    improvement = false
    
    for ĩ in problematic
        # Try each potential replacement center
        for i✶ in 1:S
            Y[i✶] == 1 && continue  # Skip already open centers
            
            if try_center_swap(ĩ, i✶)
                Y[ĩ] = 0
                Y[i✶] = 1
                improvement = true
                println("Swapped problematic center $ĩ with center $i✶")
                break
            end
        end
    end
    
    return improvement
end

function mainLocal(; path="solucion_grasp_16_625_feas.jld2")
    solution = read_solution(path)
    newSolution = localSearch(solution)
    println(isFactible(newSolution))
    println(newSolution.Weight)
    write_solution(newSolution, "sol_ls_1_625_new003_newrange.jld2")
    #plot_solution(newSolution, "plot_sol_2_1250_viejo_relax_ls.png")
    return newSolution
end


function localSearch(solution)
    oldSol = solution
    oldTime = oldSol.Time
    instance = solution.Instance
    factible_after_repair = false
    targets_lower, targets_upper = calculate_targets(instance)
    targets_lower_op, targets_upper_op = calculate_targets_optimized(instance)
    println(isFactible(oldSol))
    factible, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper)
    println("is factible4 oldsol ", isFactible4(oldSol, targets_lower, targets_upper))
    repaired = oldSol
    before_time = Dates.now()
    original_weight = 10000000000000
    if factible
        println("Factible")
        original_weight = solution.Weight
        #println(original_weight)
        factible_after_repair = true
    else
        println("Reparando")
        repaired_1 = repair_solution1(oldSol, constraints, targets_lower, targets_upper, remove, add)
        fac_repaired_1, cons = isFactible(repaired_1, true)
        #println(fac_repaired_1)
        if !fac_repaired_1
            println("con la 4")
            #=
            repaired_2 = repair_solution4(oldSol, constraints, targets_lower, targets_upper, remove, add)
            fac_repaired_2, cons = isFactible(repaired_2, true)
            factible2, constraints2, remove2, add2 = isFactible4(repaired_2, targets_lower, targets_upper)
            =#
            println("========================================================== 3")
            repaired_3 = repair_solution_v3(oldSol, targets_lower, targets_upper)
            fac_repaired_3, cons3 = isFactible(repaired_3, true)
            println("=============================================================================== 4")
            factible3, constraints3, remove3, add3 = isFactible4(repaired_3, targets_lower, targets_upper)
            repaired_4 = repair_solution4(repaired_3, constraints3, targets_lower, targets_upper, remove3, add3)
            fac_repaired_4, cons4 = isFactible(repaired_4, true)
            if fac_repaired_4
                factible_after_repair = true
                repaired = repaired_4
            end
        else
            repaired = repaired_1
            factible_after_repair = true
        end
        #println("intercambiando")
        #oldSol = interchange_bu_repair(oldSol, targets_lower, targets_upper, :bf)
        #println("simple")
        #oldSol = simple_bu_repair(oldSol, targets_lower, targets_upper, :bf)
        if factible_after_repair
            original_weight = repaired.Weight
            #println(repaired.Weight)
            println("Reparada ")
        end
    end
    if !factible_after_repair
        @error "INFACTIBLE"
        return 0
    end
    oldSol = repaired
    D_original = oldSol.Instance.D
    improvement = true
    loop = 0

    while improvement
        loop += 1
        improvement = false  # Reset the flag at the start of each loop iteration
        prev_weight = oldSol.Weight

        # Array to keep track of individual improvements
        improvements = Bool[]
        # First improvement function
        println("simple optimized bf")
        #println(@benchmark simple_bu_improve_optimized($oldSol, $targets_lower_op, $targets_upper_op, :bf))
        sol_moved_bu = simple_bu_improve_optimized(oldSol, targets_lower_op, targets_upper_op, :bf)
        new_weight_moved = sol_moved_bu.Weight
        println(new_weight_moved)

        push!(improvements, new_weight_moved < prev_weight)
        if improvements[end]
            prev_weight = new_weight_moved
            #println("En el loop loop el movimiento simple mejora con un new_weight_moved")
            oldSol = sol_moved_bu  # Update oldSol if there was an improvement
        end

        # Second improvement function

        println("interchange optimized bf")
        #println(@benchmark interchange_bu_improve_optimized($oldSol, $targets_lower_op, $targets_upper_op, :bf))
        sol_interchanged_bu = interchange_bu_improve_optimized(oldSol, targets_lower_op, targets_upper_op, :bf)
        new_weight_moved = sol_interchanged_bu.Weight
        println(new_weight_moved)
        push!(improvements, new_weight_moved < prev_weight)
        if improvements[end]
            prev_weight = new_weight_moved
            #println("En el loop loop el movimiento intercambio mejora con un new_weight_moved")
            oldSol = sol_interchanged_bu  # Update oldSol if there was an improvement
        end
        println("---------------------")

        #println(isFactible(sol_interchanged_bu, true))

        # Third improvement function

        println("deactivate ")
        #println(@benchmark deactivate_center_improve($oldSol, $targets_lower, $targets_upper))
        sol_deactivated_center = deactivate_center_improve(oldSol, targets_lower, targets_upper)
        new_weight_moved = sol_deactivated_center.Weight
        println(new_weight_moved)
        push!(improvements, new_weight_moved < prev_weight)
        if improvements[end]
            prev_weight = new_weight_moved
            #println("En el loop $loop el movimiento desactivar mejora con un $new_weight_moved")
            oldSol = sol_deactivated_center  # Update oldSol if there was an improvement
        end
        println("---------------------")
        #=
                println(" chain ")
                #println(@benchmark deactivate_center_improve($oldSol, $targets_lower, $targets_upper))
                sol_chain = chain_move_improve_optimized(oldSol, targets_lower_op, targets_upper_op, :bf, 2)
                new_weight_moved = sol_chain.Weight
                println(new_weight_moved)
                push!(improvements, new_weight_moved < prev_weight)
                if improvements[end]
                    prev_weight = new_weight_moved
                    #println("En el loop $loop el movimiento desactivar mejora con un $new_weight_moved")
                    oldSol = sol_chain  # Update oldSol if there was an improvement
                end
                =#
        #=
        println("---------------------")
        println("ejection chain ")
        #println(@benchmark deactivate_center_improve($oldSol, $targets_lower, $targets_upper))
        sol_chain = classical_ejection_chain_improve(oldSol, targets_lower_op, targets_upper_op, 4)
        new_weight_moved = sol_chain.Weight
        println(new_weight_moved)
        push!(improvements, new_weight_moved < prev_weight)
        if improvements[end]
            prev_weight = new_weight_moved
            #println("En el loop $loop el movimiento desactivar mejora con un $new_weight_moved")
            oldSol = sol_chain  # Update oldSol if there was an improvement
        end
        println("---------------------")
        =#
        improvement = any(improvements)
    end
    after_time = Dates.now()
    println(after_time - before_time)
    println(oldSol.Weight)
    return oldSol
end

if abspath(PROGRAM_FILE) == @__FILE__
    mainLocal(; path=ARGS[1])
end