using Types, BenchmarkTools, DelimitedFiles
using TimerOutputs
using LinearAlgebra
using LoopVectorization
using Dates

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
                    println("μ: ", round(Int, Y[i] * μ[m][i] * (1 - T[m])))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= round(Int, Y[i] * μ[m][i] * (1 + T[m])))
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", round(Int, Y[i] * μ[m][i] * (1 + T[m])))
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

function get_magnitude_contraints(usables_i, B, M, V, R, X, values_matrix, risk_vec, targets_lower, targets_upper, β)
    value = 0
    for i in usables_i
        for m in 1:M
            if values_matrix[i, m] > targets_upper[m]
                value += abs(targets_upper[m] - values_matrix[i, m])
            elseif values_matrix[i, m] < targets_lower[m]
                value += abs(targets_lower[m] - values_matrix[i, m])
            end
        end
        if risk_vec[i] > β
            value += abs(β - risk_vec[i])
        end
    end
    return value
end

function calculate_targets(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    targets_lower = Vector{Int64}(undef, M)
    targets_upper = Vector{Int64}(undef, M)
    for m in 1:M
        targets_lower[m] = round(Int, (1 * μ[m][1] * (1 - T[m])))
        targets_upper[m] = round(Int, (1 * μ[m][1] * (1 + T[m])))
    end
    return targets_lower, targets_upper
end

function start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
    for i in 1:S
        for m in 1:M
            values_matrix[i, m] = sum(X[i, j] * V[m][j] for j in 1:B)
        end
        risk_vec[i] = sum(X[i, j] * R[j] for j in 1:B)
    end
    return values_matrix, risk_vec
end

function start_constraints_blazing_fast(X::Matrix{Int}, V::Array{Vector{Int},1}, R::Vector{Int})
    values_matrix = X * hcat(V...)
    risk_vec = vec(sum(X .* R', dims=2))
    return values_matrix, risk_vec
end

function start_constraints_optimized(S, B, M, V, R, X, values_matrix, risk_vec)
    @inbounds for i in 1:S
        risk_vec[i] = 0
        for m in 1:M
            values_matrix[i, m] = 0
        end
        @inbounds for j in 1:B
            x_val = X[i, j]
            risk_vec[i] += x_val * R[j]
            for m in 1:M
                values_matrix[i, m] += x_val * V[m][j]
            end
        end
    end
    return values_matrix, risk_vec
end

function start_constraints_optimized_v3(S, B, M, V, R, X, values_matrix, risk_vec)
    @inbounds for i in 1:S
        risk_vec[i] = 0
        for m in 1:M
            tmp = 0
            @simd for j in 1:B
                tmp += X[i, j] * V[m][j]
            end
            values_matrix[i, m] = tmp
        end

        @simd for j in 1:B
            risk_vec[i] += X[i, j] * R[j]
        end
    end
    return values_matrix, risk_vec
end

function start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    for m in eachindex(V)
        mul!(view(values_matrix, :, m), X, V[m])
    end
    mul!(risk_vec, X, R)
    return values_matrix, risk_vec
end

@inline function prepare_values_risk_vec(S, M, values_matrix, risk_vec)
    @inbounds for i in 1:S
        risk_vec[i] = 0
        for m in 1:M
            values_matrix[i, m] = 0
        end
    end
    return values_matrix, risk_vec
end

function start_constraints_optimized_v4(S, B, M, V, R, X, values_matrix, risk_vec)
    values_matrix, risk_vec = prepare_values_risk_vec(S, M, values_matrix, risk_vec)
    @turbo for i in 1:S
        for j in 1:B
            x_val = X[i, j]
            risk_vec[i] += x_val * R[j]
            for m in 1:M
                values_matrix[i, m] += x_val * V[m][j]
            end
        end
    end
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

function maximums3(M, n)
    v = vec(M)
    l = length(v)
    ix = [1:l;]
    partialsortperm!(ix, v, (l-n+1):l, initialized=true)
    indices = CartesianIndices(M)[ix[(l-n+1):l]]
    return reverse!(indices)
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
    values_matrix = Matrix{Int64}(undef, S, M)
    risk_vec = Vector{Int64}(undef, S)

    values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
    current_constraints_value = get_magnitude_contraints(usables_i, B, M, V, R, X, values_matrix, risk_vec, targets_lower, targets_upper, β)
    cant_fix = false
    for ĩ in remove
        if cant_fix
            break
        end
        values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
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
    #current_constraints_value = get_magnitude_contraints(S, B, M, V, R, X, values_matrix, risk_vec, targets_lower, targets_upper, β)
    #println(current_constraints_value)
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
    current_constraints_value = get_magnitude_contraints(usables_i, B, M, V, R, X, values_matrix, risk_vec, targets_lower, targets_upper, β)
    newSol = Solution(instance, X, Y, Weight, solution.Time)
    return newSol
end

# no funciona por los indices en julia, no se recorre el arreglo bien

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
    values_matrix = Matrix{Int64}(undef, S, M)
    risk_vec = Vector{Int64}(undef, S)
    remove_vec = collect(remove) |> sort!
    remove_vec_fixed = collect(remove) |> sort!
    add_vec = collect(add) |> sort!
    add_vec_fixed = collect(add) |> sort!
    remove_factible = [false for i in remove_vec]
    add_factible = [false for i in add_vec]

    values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
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

function interchange_bu_improve(solution, targets_lower, targets_upper, strategy)
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
    Weight = solution.Weight
    values_matrix = Matrix{Int64}(undef, S, M)
    risk_vec = Vector{Int64}(undef, S)

    values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
    improvement = true
    while improvement
        improvement = false
        best_ĩ = 0
        best_i✶ = 0
        best_j₁ = 0
        best_j₂ = 0
        best_improvement = 0
        for j₁ in 1:B
            ĩ = findfirst(==(1), @views X[:, j₁])
            for j₂ in j₁+1:B
                i✶ = findfirst(==(1), @views X[:, j₂])
                weight_diff = -D[ĩ, j₁] - D[i✶, j₂] + D[ĩ, j₂] + D[i✶, j₁]
                if weight_diff < 0 # si el peso mejora, revisa el estado del movimiento
                    can_do_move = true
                    value_ĩ_1 = values_matrix[ĩ, 1] - V[1][j₁] + V[1][j₂]
                    value_ĩ_2 = values_matrix[ĩ, 2] - V[2][j₁] + V[2][j₂]
                    value_ĩ_3 = values_matrix[ĩ, 3] - V[3][j₁] + V[3][j₂]
                    value_i✶_1 = values_matrix[i✶, 1] + V[1][j₁] - V[1][j₂]
                    value_i✶_2 = values_matrix[i✶, 2] + V[2][j₁] - V[2][j₂]
                    value_i✶_3 = values_matrix[i✶, 3] + V[3][j₁] - V[3][j₂]
                    if value_ĩ_1 > targets_upper[1]
                        can_do_move = false
                    end
                    if value_ĩ_2 > targets_upper[2]
                        can_do_move = false
                    end
                    if value_ĩ_3 > targets_upper[3]
                        can_do_move = false
                    end
                    if value_i✶_1 > targets_upper[1]
                        can_do_move = false
                    end
                    if value_i✶_2 > targets_upper[2]
                        can_do_move = false
                    end
                    if value_i✶_3 > targets_upper[3]
                        can_do_move = false
                    end
                    if value_ĩ_1 < targets_lower[1]
                        can_do_move = false
                    end
                    if value_ĩ_2 < targets_lower[2]
                        can_do_move = false
                    end
                    if value_ĩ_3 < targets_lower[3]
                        can_do_move = false
                    end
                    if value_i✶_1 < targets_lower[1]
                        can_do_move = false
                    end
                    if value_i✶_2 < targets_lower[2]
                        can_do_move = false
                    end
                    if value_i✶_3 < targets_lower[3]
                        can_do_move = false
                    end
                    risk_ĩ = risk_vec[ĩ] - R[j₁] + R[j₂]
                    risk_i✶ = risk_vec[i✶] + R[j₁] - R[j₂]
                    if risk_ĩ > β
                        # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                        can_do_move = false
                    end
                    if risk_i✶ > β
                        can_do_move = false
                    end
                    if can_do_move
                        newWeight = Weight + weight_diff
                        abs_improvement = -weight_diff #inviertele el signo para hacerlo positivo
                        if newWeight < Weight
                            if strategy == :ff
                                #println("Intercambiando $ĩ por $i✶ la j $j₁ por la j $j₂")
                                X[ĩ, j₁] = 0
                                X[i✶, j₂] = 0
                                X[ĩ, j₂] = 1
                                X[i✶, j₁] = 1
                                improvement = true
                                Weight = newWeight
                                for m in 1:M
                                    values_matrix[ĩ, m] -= V[m][j₁] # restale a ĩ, previo
                                    values_matrix[i✶, m] -= V[m][j₂] # restale a i✶, previo
                                    values_matrix[ĩ, m] += V[m][j₂] # sumale a ĩ, nuevo
                                    values_matrix[i✶, m] += V[m][j₁] # sumale a i✶, nuevo
                                end
                                risk_vec[ĩ] -= R[j₁] # restale a ĩ el viejo
                                risk_vec[i✶] -= R[j₂] # restale a i✶ el viejo
                                risk_vec[ĩ] += R[j₂] # sumale a i✶ el nuevo
                                risk_vec[i✶] += R[j₁] # sumale a ĩ el nuevo
                                aux = ĩ
                                ĩ = i✶
                                i✶ = aux
                                #newSol = Solution(instance, X, Y, Weight, solution.Time)
                                #println(isFactible(newSol, true))
                            else
                                if abs_improvement > best_improvement
                                    best_ĩ = ĩ
                                    best_i✶ = i✶
                                    best_j₁ = j₁
                                    best_j₂ = j₂
                                    best_improvement = abs_improvement
                                end
                            end
                        end
                    end
                end
            end
        end
        if strategy == :bf
            if best_ĩ != 0
                X[best_ĩ, best_j₁] = 0
                X[best_i✶, best_j₂] = 0
                X[best_ĩ, best_j₂] = 1
                X[best_i✶, best_j₁] = 1
                indices = findall(x -> x == 1, X)
                newWeight = 0
                for indice in indices
                    newWeight += instance.D[indice]
                end
                Weight = newWeight
                improvement = true
                for m in 1:M
                    values_matrix[best_ĩ, m] -= V[m][best_j₁] # restale a ĩ, previo
                    values_matrix[best_i✶, m] -= V[m][best_j₂] # restale a i✶, previo
                    values_matrix[best_ĩ, m] += V[m][best_j₂] # sumale a ĩ, nuevo
                    values_matrix[best_i✶, m] += V[m][best_j₁] # sumale a i✶, nuevo
                end
                risk_vec[best_ĩ] -= R[best_j₁] # restale a ĩ el viejo
                risk_vec[best_i✶] -= R[best_j₂] # restale a i✶ el viejo
                risk_vec[best_ĩ] += R[best_j₂] # sumale a i✶ el nuevo
                risk_vec[best_i✶] += R[best_j₁] # sumale a ĩ el nuevo
            end
        end
    end
    indices = findall(x -> x == 1, X)
    Weight = 0
    for indice in indices
        Weight += instance.D[indice]
    end
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

function get_best_clients_for_centers(D, N)
    num_centers, num_clients = size(D)
    best_clients = Dict{Int,Vector{Int64}}()
    for i in 1:num_centers
        # Use a temporary array to store the client opportunity costs for this center
        costs = Tuple{Int64,Int}[]
        for j in 1:num_clients
            push!(costs, (D[i, j], j))
        end
        # Sort the costs
        sort!(costs)
        # Store the top N clients for this center
        best_clients[i] = [cost[2] for cost in costs[1:N]] # Extract only the client index
    end
    return best_clients
end

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

function deactivate_center_improve(solution, targets_lower, targets_upper, strategy=:bf)
    to = TimerOutput()
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
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    P = instance.P
    Weight = solution.Weight
    not_usables_i = Set(findall(==(0), Y))
    usables_i = Set(findall(==(1), Y))
    values_matrix = Matrix{Int64}(undef, S, M)
    risk_vec = Vector{Int64}(undef, S)
    values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
    @timeit to "Get best assignments clients" best_assignments_clients = get_best_assignments(D, P)
    @timeit to "Get best clients for centers" best_clients_for_centers = get_best_clients_for_centers(D, B)
    count = count_k(usables_i, Sk)
    improvement = true
    total_time_spent = 0
    @timeit to "Main while" while improvement
        improvement = false
        for ĩ in usables_i
            for i✶ in not_usables_i
                start = now()
                #@timeit to "Start constraints" values_matrix, risk_vec = start_constraints_optimized_v3(S, B, M, V, R, X, values_matrix, risk_vec)
                @timeit to "Start constraints" values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
                #@timeit to "Start contraints" values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
                finish = now()
                delta = finish - start
                delta_millis = round(delta, Millisecond)
                delta_val = delta_millis.value
                total_time_spent += delta_val
                #@timeit to "Calc start constraints" values_matrix, risk_vec = start_constraints_blazing_fast(X, V, R)
                #@timeit to "Calc start constraints" values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
                modified_X = Dict{Tuple{Int,Int},Int}()
                useful = true
                ĩₖ = node_type(ĩ, Sk)
                i✶ₖ = node_type(i✶, Sk)
                count_ĩ = count[ĩₖ] - 1
                count_i✶ = count[i✶ₖ] + 1
                if count_ĩ <= Uk[ĩₖ] && count_ĩ >= Lk[ĩₖ] && count_i✶ <= Uk[i✶ₖ] && count_i✶ >= Lk[i✶ₖ]
                    js_assigned = findall(==(1), @views X[ĩ, :])
                    for j in js_assigned
                        if !haskey(modified_X, (ĩ, j))
                            modified_X[(ĩ, j)] = X[ĩ, j]
                        end
                    end
                    X[ĩ, :] .= 0
                    js_assigned_set = Set(js_assigned)
                    weight_old_branch = sum(D[ĩ, js_assigned]) # total weight of old branch/center
                    weight_new_branch = 0
                    factible_yet = false
                    fulls_m = zeros(Int, M)
                    @timeit to "Asignar clientes a i✶" for client in best_clients_for_centers[i✶]
                        potential_assignment_valid = true
                        previous_i_client = findfirst(==(1), @views X[:, client])
                        if previous_i_client !== nothing # si es nothing entonces estaba asignado a ĩ
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
                                    values_matrix[previous_i_client, m] -= V[m][client]
                                    values_matrix[i✶, m] += V[m][client]
                                    if values_matrix[i✶, m] > targets_lower[m]
                                        fulls_m[m] = 1
                                    end
                                end
                                risk_vec[previous_i_client] -= R[client]
                                risk_vec[i✶] += R[client]
                            end
                        else
                            for m in 1:M
                                values_matrix[i✶, m] += V[m][client]
                                if values_matrix[i✶, m] > targets_lower[m]
                                    fulls_m[m] = 1
                                end
                            end
                            risk_vec[i✶] += R[client]
                        end
                        if potential_assignment_valid
                            if !haskey(modified_X, (i✶, client))
                                modified_X[(i✶, client)] = X[i✶, client]
                            end
                            if previous_i_client !== nothing
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
                                factible_yet = true
                                # ya llenamos el centro i✶
                                break
                            end
                        end
                    end
                    useful = true
                    if !factible_yet
                        useful = false
                    end
                    @timeit to "Asignar clientes huerfanos" for orphaned_client in js_assigned_set
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
                                        values_matrix[center, m] += V[m][orphaned_client]
                                    end
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
                    else
                        for item in modified_X
                            i, j = item[1]
                            val = item[2]
                            X[i, j] = val
                        end
                    end
                end
            end
        end
    end
    indices = findall(x -> x == 1, X)
    Weight = 0
    for indice in indices
        Weight += instance.D[indice]
    end
    show(to)
    #@show total_time_spent
    return Solution(instance, X, Y, Weight, solution.Time)
end

function deactivate_center_improve2(solution, targets_lower, targets_upper, strategy=:bf)
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
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    P = instance.P
    Weight = solution.Weight
    n = P
    not_usables_i = Set(findall(==(0), Y))
    usables_i = Set(findall(==(1), Y))
    values_matrix = Matrix{Int64}(undef, S, M)
    risk_vec = Vector{Int64}(undef, S)

    for i in not_usables_i
        D[i, :] .= 10000000000000
    end

    values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
    count = count_k(usables_i, Sk)
    improvement = true
    # hacerlo solo FF

    while improvement
        improvement = false
        for ĩ in usables_i
            for i✶ in not_usables_i
                ij_prev_assignments = []
                useful = true
                ĩₖ = node_type(ĩ, Sk)
                i✶ₖ = node_type(i✶, Sk)
                count_ĩ = count[ĩₖ] - 1
                count_i✶ = count[i✶ₖ] + 1
                if count_ĩ <= Uk[ĩₖ] && count_ĩ >= Lk[ĩₖ] && count_i✶ <= Uk[i✶ₖ] && count_i✶ >= Lk[i✶ₖ]
                    # tenemos que mantener cuales son las asignaciones de ĩ antes de apagarlas
                    # para esto las guardo en un arreglo mejor
                    js_assigned = findall(==(1), X[ĩ, :])
                    X[ĩ, :] .= 0
                    js_assigned_set = Set(js_assigned)
                    weight_old_branch = sum(D[ĩ, js_assigned]) # peso total representado por la rama
                    weight_new_branch = 0
                    D[i✶, :] .= instance.D[i✶, :]
                    D[ĩ, :] .= 10000000000000
                    factible_yet = false
                    candidates_BUs = minimums2(D[i✶, :], B)
                    for m in 1:M
                        values_matrix[ĩ, m] = 0
                    end
                    risk_vec[ĩ] = 0
                    while !factible_yet
                        if length(candidates_BUs) == 0
                            useful = false
                            break
                        end
                        j = popfirst!(candidates_BUs)[1]
                        i_old = findfirst(==(1), X[:, j]) # asignacion previa de esta j
                        if i_old === nothing
                            i_old = ĩ
                        end
                        factible_yet = true
                        can_do_move = true
                        for m in 1:M
                            values_matrix[i_old, m] -= V[m][j] # restale a ĩ, previo
                            values_matrix[i✶, m] += V[m][j] # sumale a i, nuevo
                            if values_matrix[i✶, m] > targets_upper[m]
                                # no deberia de pasar porque entonces la infactibilidad cambia de razon
                                factible_yet = false
                                can_do_move = false
                            end
                            if values_matrix[i✶, m] < targets_lower[m]
                                factible_yet = false
                            end
                            if i_old != ĩ
                                if values_matrix[i_old, m] < targets_lower[m]
                                    factible_yet = false
                                    can_do_move = false
                                    # no deberia de pasar porque entonces i_old es infactible ahora
                                end
                            end
                        end
                        risk_vec[i_old] -= R[j] # restale a ĩ el viejo
                        risk_vec[i✶] += R[j] # sumale a i el nuevo
                        if risk_vec[i✶] > β
                            # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                            can_do_move = false
                            factible_yet = false
                        end
                        if can_do_move
                            X[i✶, j] = 1
                            X[i_old, j] = 0
                            if j ∉ js_assigned
                                push!(ij_prev_assignments, (i_old, j))
                            else
                                delete!(js_assigned_set, j) # ya la asignamos
                            end
                            weight_new_branch += D[i✶, j]
                        else
                            for m in 1:M
                                values_matrix[i_old, m] += V[m][j] # corrige el valor de la NO asignacion
                                values_matrix[i✶, m] -= V[m][j]
                            end
                            risk_vec[i_old] += R[j]
                            risk_vec[i✶] -= R[j]
                        end
                    end


                    for j in js_assigned_set # para las js que aun queden sin asignar
                        candidates_is = minimums2(D[:, j], n)
                        j_assigned = false
                        if !useful
                            break
                        end
                        while !j_assigned
                            i = 0
                            if length(candidates_is) == 0 # si no podemos asignar esta j, entonces el movimiento se invalida
                                useful = false
                                break
                            else
                                i = popfirst!(candidates_is)
                                i = i[1]
                            end
                            can_do_move = true
                            #println("probando cambio en $i por $ĩ")
                            for m in 1:M
                                values_matrix[i, m] += V[m][j] # sumale a i, nuevo
                                if values_matrix[i, m] > targets_upper[m]
                                    #println("no es factible por que se pasa el otro upper: ", values_matrix[i, m])
                                    can_do_move = false
                                end
                            end
                            risk_vec[i] += R[j] # sumale a i el nuevo
                            if risk_vec[i] > β
                                # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                                can_do_move = false
                            end
                            if can_do_move
                                X[i, j] = 1
                                j_assigned = true
                                weight_new_branch += D[i, j]
                            else
                                # si no puedo hacer el movimiento, restaura el valor de la ev parcial
                                for m in 1:M
                                    values_matrix[i, m] -= V[m][j]
                                end
                                risk_vec[i] -= R[j]
                            end
                        end
                    end

                    if weight_new_branch < weight_old_branch && useful
                        # arreglar esto
                    else
                        useful = false
                    end

                    if useful
                        Y[i✶] = 1
                        Y[ĩ] = 0
                        count[ĩₖ] -= 1
                        count[i✶ₖ] += 1
                        delete!(not_usables_i, i✶)
                        delete!(usables_i, ĩ)
                        push!(usables_i, i✶)
                        push!(not_usables_i, ĩ)
                        weight_old_branch = weight_new_branch
                        indices = findall(x -> x == 1, X)
                        Weight = 0
                        for indice in indices
                            Weight += instance.D[indice]
                        end
                        aux = ĩ
                        ĩ = i✶
                        i✶ = aux
                        improvement = true
                    else
                        # restaurar el valor de values_matrix y risk_vec
                        # devolver el valor de las demas afectadas
                        for m in 1:M
                            values_matrix[ĩ, m] = sum(V[m][js_assigned])
                        end
                        risk_vec[ĩ] = sum(R[js_assigned])
                        for j in js_assigned
                            i = findfirst(==(1), X[:, j])
                            if i === nothing
                                i = ĩ
                            end
                            for m in 1:M
                                values_matrix[i, m] -= V[m][j]
                            end
                            risk_vec[i] -= R[j]
                            X[:, j] .= 0 # deshacemos asignaciones de los clientes
                        end

                        X[i✶, :] .= 0 #deshaz cualquier asignacion al nuevo centro

                        for (i_old, j_old) in ij_prev_assignments
                            X[i_old, j_old] = 1 # pon la asignacion previa de vuelta
                        end

                        for m in 1:M
                            values_matrix[i✶, m] = 0
                        end
                        risk_vec[i✶] = 0
                        X[ĩ, js_assigned] .= 1 # haz la asignacion vieja de nuevo
                        D[i✶, :] .= 10000000000000
                        D[ĩ, :] .= instance.D[ĩ, :]
                    end
                end
            end
        end
    end

    indices = findall(x -> x == 1, X)
    Weight = 0
    for indice in indices
        Weight += instance.D[indice]
    end
    return Solution(instance, X, Y, Weight, solution.Time)
end

function simple_bu_improve(solution, targets_lower, targets_upper, strategy)
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
    Weight = solution.Weight
    n = round(Int, (P / 2))
    not_usables_i = Set(findall(==(0), Y))
    usables_i = Set(findall(==(1), Y))
    values_matrix = Matrix{Int64}(undef, S, M)
    risk_vec = Vector{Int64}(undef, S)

    for j in not_usables_i
        D[j, :] .= 10000000000000
    end

    values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
    improvement = true

    while improvement
        improvement = false
        best_new_i = 0
        best_old_i = 0
        best_j = 0
        best_improvement = 0
        for j in 1:B
            ĩ = findfirst(==(1), @views X[:, j]) # asignacion previa ĩ, j = 1
            min_index = argmin(D[:, j]) # minimo de j
            if ĩ == min_index
                #println("saltandome mejor asignacion")
                continue # no tiene caso intentar mover una BU que ya tiene asignacion optima
            end
            usables_i_without_ĩ = setdiff(usables_i, ĩ)
            for i in usables_i_without_ĩ
                weight_diff = -D[ĩ, j] + D[i, j]
                if weight_diff < 0 # si el peso mejora, revisa el estado del movimiento
                    can_do_move = true
                    value_ĩ_1 = values_matrix[ĩ, 1] - V[1][j]
                    value_ĩ_2 = values_matrix[ĩ, 2] - V[2][j]
                    value_ĩ_3 = values_matrix[ĩ, 3] - V[3][j]
                    value_i_1 = values_matrix[i, 1] + V[1][j]
                    value_i_2 = values_matrix[i, 2] + V[2][j]
                    value_i_3 = values_matrix[i, 3] + V[3][j]
                    if value_ĩ_1 > targets_upper[1]
                        can_do_move = false
                    end
                    if value_ĩ_2 > targets_upper[2]
                        can_do_move = false
                    end
                    if value_ĩ_3 > targets_upper[3]
                        can_do_move = false
                    end
                    if value_i_1 > targets_upper[1]
                        can_do_move = false
                    end
                    if value_i_2 > targets_upper[2]
                        can_do_move = false
                    end
                    if value_i_3 > targets_upper[3]
                        can_do_move = false
                    end
                    if value_ĩ_1 < targets_lower[1]
                        can_do_move = false
                    end
                    if value_ĩ_2 < targets_lower[2]
                        can_do_move = false
                    end
                    if value_ĩ_3 < targets_lower[3]
                        can_do_move = false
                    end
                    if value_i_1 < targets_lower[1]
                        can_do_move = false
                    end
                    if value_i_2 < targets_lower[2]
                        can_do_move = false
                    end
                    if value_i_3 < targets_lower[3]
                        can_do_move = false
                    end
                    risk_ĩ = risk_vec[ĩ] - R[j]
                    risk_i = risk_vec[i] + R[j]
                    if risk_ĩ > β
                        # si yo le agrego, no deberia de pasar esto porque cambia infactibilidad
                        can_do_move = false
                    end
                    if risk_i > β
                        can_do_move = false
                    end
                    if can_do_move
                        newWeight = Weight + weight_diff
                        abs_improvement = -weight_diff
                        if newWeight < Weight
                            if strategy == :ff # first found
                                X[ĩ, j] = 0
                                X[i, j] = 1
                                improvement = true
                                for m in 1:M
                                    values_matrix[ĩ, m] -= V[m][j] # restale a ĩ, previo
                                    values_matrix[i, m] += V[m][j] # restale a i✶, previo
                                end
                                risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
                                risk_vec[i] += R[j] # restale a i✶ el viejo
                                Weight = newWeight
                                ĩ = i # duh
                            elseif strategy == :bf
                                if abs_improvement > best_improvement
                                    best_new_i = i
                                    best_old_i = ĩ
                                    best_j = j
                                    best_improvement = abs_improvement
                                end
                            end
                        end
                    end
                end
            end
        end
        if strategy == :bf
            if best_new_i != 0
                X[best_new_i, best_j] = 1
                X[best_old_i, best_j] = 0
                for m in 1:M
                    values_matrix[best_old_i, m] -= V[m][best_j] # restale a ĩ, previo
                    values_matrix[best_new_i, m] += V[m][best_j] # restale a i✶, previo
                end
                risk_vec[best_new_i] += R[best_j]
                risk_vec[best_old_i] -= R[best_j]
                improvement = true
            else
                improvement = false
            end
        end
    end
    indices = findall(x -> x == 1, X)
    Weight = 0
    for indice in indices
        Weight += instance.D[indice]
    end
    return Solution(instance, X, Y, Weight, solution.Time)
end

function mainLocal(; path="sol_1_625_78_32_relax_opp.jld2")
    solution = read_solution(path)
    newSolution = @time localSearch(solution)
    write_solution(newSolution, "sol_ls_1_312.jld2")
    plot_solution(newSolution, "plot_sol_2_1250_viejo_relax_ls.png")
    return newSolution
end

function localSearch(solution)
    oldSol = solution
    oldTime = oldSol.Time
    instance = solution.Instance
    factible_after_repair = false
    targets_lower, targets_upper = calculate_targets(instance)
    factible, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper)
    repaired = oldSol
    println(remove)
    println(add)
    original_weight = 10000000000000
    if factible
        println("Factible")
        original_weight = solution.Weight
        println(original_weight)
        factible_after_repair = true
    else
        println("Reparando")
        repaired_1 = repair_solution1(oldSol, constraints, targets_lower, targets_upper, remove, add)
        fac_repaired_1, cons = isFactible(repaired_1, false)
        println(fac_repaired_1)
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
            println(repaired.Weight)
            println("Reparada ")
        end
    end
    if !factible_after_repair
        @error "INFACTIBLE"
        return 0
    end
    oldSol = repaired
    improvement = true
    loop = 0
    while improvement
        improvement = false  # Reset the flag at the start of each loop iteration
        prev_weight = oldSol.Weight

        # Array to keep track of individual improvements
        improvements = Bool[]
        """
        # First improvement function
        sol_moved_bu = simple_bu_improve(oldSol, targets_lower, targets_upper, :bf)
        new_weight_moved = sol_moved_bu.Weight
        push!(improvements, new_weight_moved < prev_weight)
        if improvements[end]
            prev_weight = new_weight_moved
            #println("En el loop loop el movimiento simple mejora con un new_weight_moved")
            oldSol = sol_moved_bu  # Update oldSol if there was an improvement
        end
        #println(isFactible(sol_moved_bu, true))
        """
        """
        # Second improvement function
        sol_interchanged_bu = interchange_bu_improve(oldSol, targets_lower, targets_upper, :bf)
        new_weight_moved = sol_interchanged_bu.Weight
        push!(improvements, new_weight_moved < prev_weight)
        if improvements[end]
            prev_weight = new_weight_moved
            #println("En el loop loop el movimiento intercambio mejora con un new_weight_moved")
            oldSol = sol_interchanged_bu  # Update oldSol if there was an improvement
        end
        """
        #println(isFactible(sol_interchanged_bu, true))

        # Third improvement function
        sol_deactivated_center = deactivate_center_improve(oldSol, targets_lower, targets_upper)
        new_weight_moved = sol_deactivated_center.Weight
        push!(improvements, new_weight_moved < prev_weight)
        if improvements[end]
            prev_weight = new_weight_moved
            #println("En el loop $loop el movimiento desactivar mejora con un $new_weight_moved")
            oldSol = sol_deactivated_center  # Update oldSol if there was an improvement
        end
        #println(isFactible(sol_deactivated_center, true))

        # Check for any improvements

        improvement = any(improvements)
        #println(isFactible(oldSol, true))
    end

    println(isFactible(oldSol, true))
    println(oldSol.Weight)
    println("\a")
    return oldSol
end

if abspath(PROGRAM_FILE) == @__FILE__
    mainLocal(; path=ARGS[1])
else
    mainLocal(; path="sol_1_625_78_32_pdisp_opp.jld2")
end