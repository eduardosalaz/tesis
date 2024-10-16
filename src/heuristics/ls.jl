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
            if !(Y[i] * μ[m][i] * (1 - T[m]) <= sum(X[i, j] * V[m][j] for j in 1:B))
                if verbose
                    println("violando V inferior en i: $i y m: $m")
                    println("μ: ", (Y[i] * μ[m][i] * (1 - T[m])))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= Y[i] * μ[m][i] * (1 + T[m]))
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", (Y[i] * μ[m][i] * (1 + T[m])))
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
    targets_lower = Vector{Float32}(undef, M)
    targets_upper = Vector{Float32}(undef, M)
    for m in 1:M
        targets_lower[m] = (1 * μ[m][1] * (1 - T[m]))
        targets_upper[m] = (1 * μ[m][1] * (1 + T[m]))
    end
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
    β₁ = β[1]
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
            @error "X NO SERVIDA"
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
            if (targets_lower[m] > sum(X[i, j] * V[m][j] for j in 1:B))
                if verbose
                    println("violando V inferior en i: $i y m: $m")
                    println("μ: ", targets_lower[m])
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                push!(add, i)
                number_constraints_violated += 1
            end
            if (sum(X[i, j] * V[m][j] for j in 1:B) > targets_upper[m])
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", targets_upper[m])
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                push!(remove, i)
                number_constraints_violated += 1
            end
        end
    end

    for i in usables_i
        if sum(X[i, j] * R[j] for j in 1:B) > β₁
            if verbose
                println("violando riesgo en $i")
                println("β: ", β₁)
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
    β = instance.β[1]
    P = instance.P
    n = round(Int, (P - 1))
    usables_i = findall(==(1), Y)
    not_usables_i = Set(findall(==(0), Y))
    values_matrix = Matrix{Float32}(undef, S, M)
    risk_vec = Vector{Float32}(undef, S)

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
                    if values_matrix[i, m] > targets_upper[m]
                        #println("no es factible por que se pasa el otro upper: ", values_matrix[i, m])
                        factible_yet = false
                    end
                    if values_matrix[ĩ, m] < targets_lower[m]
                        factible_yet = false
                        #println("no es factible por que se baja el lower: ", values_matrix[ĩ, m], " para ", targets_lower[m])
                        can_do_move = false
                        # no deberia de pasar porque entonces ĩ es infactible ahora
                    end
                end
                risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
                risk_vec[i] += R[j] # sumale a i el nuevo
                if risk_vec[i] > β
                    # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                    can_do_move = false
                    #println("no es factible porque se pasa el otro risk: ", risk_vec[i])
                    factible_yet = false
                end
                if risk_vec[ĩ] > β
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
                values_matrix[ĩ, m] -= V[m][j] # restale a ĩ, previo
                values_matrix[i, m] += V[m][j] # sumale a i, nuevo
                if values_matrix[i, m] > targets_upper[m]
                    # no deberia de pasar porque entonces la infactibilidad cambia de razon
                    factible_yet = false
                    can_do_move = false
                end
                if values_matrix[i, m] < targets_lower[m]
                    factible_yet = false
                end
                if values_matrix[ĩ, m] < targets_lower[m]
                    factible_yet = false
                    can_do_move = false
                    # no deberia de pasar porque entonces ĩ es infactible ahora
                end
            end
            risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
            risk_vec[i] += R[j] # sumale a i el nuevo
            if risk_vec[i] > β
                # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                can_do_move = false
                factible_yet = false
            end
            if can_do_move
                X[ĩ, j] = 0
                X[i, j] = 1
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
    β = instance.β[1]
    P = instance.P
    n = round(Int, (P - 1))
    not_usables_i = Set(findall(==(0), Y))
    usables_i = Set(findall(==(1), Y))
    values_matrix = Matrix{Float32}(undef, S, M)
    risk_vec = Vector{Float32}(undef, S)
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
        if risk_vec[i] > β
            # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
            can_do_move = false
            #println("no es factible porque se pasa el otro risk: ", risk_vec[i])
            factible_yet = false
        end
        if risk_vec[ĩ] > β
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
                if values_matrix[i, m] > targets_upper[m]
                    #println("no es factible por que se pasa el otro upper: ", values_matrix[i, m])
                    factible_yet = false
                end
                if values_matrix[ĩ, m] < targets_lower[m]
                    factible_yet = false
                    #println("no es factible por que se baja el lower: ", values_matrix[ĩ, m], " para ", targets_lower[m])
                    can_do_move = false
                    razon = 1
                    # no deberia de pasar porque entonces ĩ es infactible ahora
                end
                if values_matrix[i, m] < targets_lower[m]
                    factible_yet = false
                end
            end
            risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
            risk_vec[i] += R[j] # sumale a i el nuevo
            if risk_vec[i] > β
                # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                can_do_move = false
                razon = 2
                #println("no es factible porque se pasa el otro risk: ", risk_vec[i])
                factible_yet = false
            end
            if risk_vec[ĩ] > β
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

# Optimized Interchange Move Implementation

function can_do_interchange_optimized(values_matrix::Matrix{Float32}, V::Vector{Vector{Int}}, ĩ::Int, i✶::Int, j₁::Int, j₂::Int, risk_vec::Vector{Float32}, R::Vector{Int}, β::Int, targets_upper::SVector{3,Float32}, targets_lower::SVector{3,Float32})
    @inbounds for m in 1:3
        value_ĩ = values_matrix[ĩ, m] - V[m][j₁] + V[m][j₂]
        value_i✶ = values_matrix[i✶, m] + V[m][j₁] - V[m][j₂]
        if value_ĩ > targets_upper[m] || value_ĩ < targets_lower[m] ||
           value_i✶ > targets_upper[m] || value_i✶ < targets_lower[m]
            return false
        end
    end
    risk_ĩ = risk_vec[ĩ] - R[j₁] + R[j₂]
    risk_i✶ = risk_vec[i✶] + R[j₁] - R[j₂]
    return risk_ĩ <= β && risk_i✶ <= β
end

function find_best_interchange_optimized(X::Matrix{Int}, D::Matrix{Int}, V::Vector{Vector{Int}}, R::Vector{Int}, B::Int, M::Int, values_matrix::Matrix{Float32}, risk_vec::Vector{Float32}, targets_lower::SVector{3,Float32}, targets_upper::SVector{3,Float32}, β::Int, strategy::Symbol)
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

function apply_interchange_optimized!(X::Matrix{Int}, values_matrix::Matrix{Float32}, risk_vec::Vector{Float32}, move, V::Vector{Vector{Int}}, R::Vector{Int})
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

function interchange_bu_improve_optimized(solution, targets_lower::SVector{3,Float32}, targets_upper::SVector{3,Float32}, strategy::Symbol)
    X, Y = solution.X, solution.Y
    instance = solution.Instance
    B, S, M = instance.B, instance.S, instance.M
    V, R, β = instance.V, instance.R, instance.β[1]
    D = instance.D
    Weight = solution.Weight

    values_matrix = Matrix{Float32}(undef, S, M)
    risk_vec = Vector{Float32}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)

    while true
        best_move = find_best_interchange_optimized(X, D, V, R, B, M, values_matrix, risk_vec, targets_lower, targets_upper, β, strategy)
        best_move === nothing && break
        apply_interchange_optimized!(X, values_matrix, risk_vec, best_move, V, R)
        Weight += best_move.weight_diff
    end

    Weight = dot(X, D)  # Recalculate weight to ensure accuracy
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
    B, S, M, V, R, β, P = instance.B, instance.S, instance.M, instance.V, instance.R, instance.β[1], instance.P
    # D = copy(instance.D)
    D = instance.D
    Sk = instance.Sk
    K = 5
    Lk = instance.Lk
    Uk = instance.Uk
    Weight = solution.Weight
    not_usables_i = BitSet(findall(==(0), Y))
    usables_i = BitSet(findall(==(1), Y))
    values_matrix = Matrix{Float32}(undef, S, M)
    risk_vec = Vector{Float32}(undef, S)
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
                modified_values_matrix = Dict{Tuple{Int,Int},Float32}()
                modified_risk = Dict{Int,Float32}()
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
                                if values_matrix[previous_i_client, m] - V[m][client] < targets_lower[m]
                                    potential_assignment_valid = false
                                    break
                                end
                            end
                            if risk_vec[i✶] + R[client] > β
                                potential_assignment_valid = false
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
                                    if values_matrix[i✶, m] > targets_lower[m]
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
                                ##=
                                if !haskey(modified_values_matrix, (i✶, m))
                                    modified_values_matrix[(i✶, m)] = values_matrix[i✶, m]
                                end
                                #
                                values_matrix[i✶, m] += V[m][client]
                                if values_matrix[i✶, m] > targets_lower[m]
                                    fulls_m[m] = 1
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
                                    if values_matrix[center, m] + V[m][orphaned_client] > targets_upper[m]
                                        potential_assignment_valid = false
                                        break
                                    end
                                end
                                if risk_vec[center] + R[orphaned_client] > β
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
    return Solution(instance, X, Y, Weight, solution.Time)
end

function calculate_targets_optimized(instance)
    μ = instance.μ
    T = instance.T

    targets_lower = SVector{3,Float32}(
        1 * μ[1][1] * (1 - T[1]),
        1 * μ[2][1] * (1 - T[2]),
        1 * μ[3][1] * (1 - T[3])
    )

    targets_upper = SVector{3,Float32}(
        1 * μ[1][1] * (1 + T[1]),
        1 * μ[2][1] * (1 + T[2]),
        1 * μ[3][1] * (1 + T[3])
    )

    return targets_lower, targets_upper
end

function can_do_move_simple_optimized(values_matrix::Matrix{Float32}, V::Vector{Vector{Int64}}, i::Int, ĩ::Int, j::Int, risk_vec::Vector{Float32}, R::Vector{Int64}, β::Int64, targets_upper::SVector{3,Float32}, targets_lower::SVector{3,Float32})
    @inbounds for m in 1:3
        valueĩ = values_matrix[ĩ, m] - V[m][j]
        value_i = values_matrix[i, m] + V[m][j]
        if valueĩ > targets_upper[m] || valueĩ < targets_lower[m] || value_i > targets_upper[m] || value_i < targets_lower[m]
            return false
        end
    end
    riskĩ = risk_vec[ĩ] - R[j]
    risk_i = risk_vec[i] + R[j]
    return riskĩ <= β && risk_i <= β
end
function simple_bu_improve_optimized(solution, targets_lower::SVector{3,Float32}, targets_upper::SVector{3,Float32}, strategy::Symbol)
    X, Y = solution.X, solution.Y
    instance = solution.Instance
    B, S, M = instance.B, instance.S, instance.M
    V, R, β = instance.V, instance.R, instance.β[1]
    D = instance.D
    Weight = solution.Weight
    usables_i = findall(==(1), Y)
    values_matrix = Matrix{Float32}(undef, S, M)
    risk_vec = Vector{Float32}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    while true
        best_move = find_best_move_simple_optimized(X, D, V, R, B, M, usables_i, values_matrix, risk_vec, targets_lower, targets_upper, β, strategy)
        best_move === nothing && break
        apply_move_simple_optimized!(X, values_matrix, risk_vec, best_move, V, R)
        Weight += best_move.weight_diff
    end
    Weight = dot(X, D)
    return Solution(instance, X, Y, Weight, solution.Time)
end
function find_best_move_simple_optimized(X::Matrix{Int}, D::Matrix{Int}, V::Vector{Vector{Int}}, R::Vector{Int}, B::Int, M::Int, usables_i::Vector{Int}, values_matrix::Matrix{Float32}, risk_vec::Vector{Float32}, targets_lower::SVector{3,Float32}, targets_upper::SVector{3,Float32}, β::Int64, strategy::Symbol)
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
function apply_move_simple_optimized!(X::Matrix{Int}, values_matrix::Matrix{Float32}, risk_vec::Vector{Float32}, move, V::Vector{Vector{Int}}, R::Vector{Int})
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

function mainLocal(; path="solucion_grasp_16_625_feas.jld2")
    solution = read_solution(path)
    newSolution = localSearch(solution)
    println(isFactible(newSolution))
    println(newSolution.Weight)
    write_solution(newSolution, "sol_ls_1250_11.jld2")
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

    original_weight = 10000000000000
    if factible
        println("Factible")
        original_weight = solution.Weight
        #println(original_weight)
        factible_after_repair = true
    else
        println("Reparando")
        repaired_1 = repair_solution1(oldSol, constraints, targets_lower, targets_upper, remove, add)
        fac_repaired_1, cons = isFactible(repaired_1, false)
        #println(fac_repaired_1)
        if !fac_repaired_1
            repaired_2 = repair_solution2(oldSol, constraints, targets_lower, targets_upper, remove, add)
            fac_repaired_2, cons = isFactible(repaired_2, false)
            if fac_repaired_2
                factible_after_repair = true
                repaired = repaired_2
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
        println(@benchmark deactivate_center_improve($oldSol, $targets_lower, $targets_upper))
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
        improvement = any(improvements)
    end
    println(oldSol.Weight)
    return oldSol
end

if abspath(PROGRAM_FILE) == @__FILE__
    mainLocal(; path=ARGS[1])
end