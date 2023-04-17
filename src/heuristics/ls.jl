using Types

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
            values_matrix[i,m] = sum(X[i,j] * V[m][j] for j in 1:B)
        end
        risk_vec[i] = sum(X[i,j] * R[j] for j in 1:B)
    end
    return values_matrix, risk_vec
end

function update_constraints(original_cons, M, V, R, ĩ, i, j, values_matrix, risk_vec, targets_lower, targets_upper, β)
    constraints = original_cons
    for m in 1:M
        values_matrix[ĩ,m] -= V[m][j] # restale a ĩ, previo
        values_matrix[i,m] += V[m][j] # sumale a i, nuevo
        if values_matrix[i,m] > targets_upper[m]
            constraints += 1
        end
        if values_matrix[i,m] < targets_lower[m]
            constraints += 1
        end
        if values_matrix[ĩ,m] < targets_lower[m]
            constraints += 1
            # naturalmente si al hacer el cambio, el ĩ incumple con lower, no lo queremos hacer
        end
    end
    risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
    risk_vec[i] += R[j] # sumale a i el nuevo
    if risk_vec[i] > β
        constraints += 1
    end
    return constraints
end

function checkFactibilityCenter(M, V, R, ĩ, i, j, values_matrix, risk_vec, targets_lower, targets_upper, β)
    factible = false
    for m in 1:M
        values_matrix[ĩ,m] -= V[m][j] # restale a ĩ, previo
        values_matrix[i,m] += V[m][j] # sumale a i, nuevo
        if values_matrix[i,m] < targets_upper[m]
            factible = true
        end
        if values_matrix[i,m] > targets_lower[m]
            factible = true
        end
        if values_matrix[ĩ,m] > targets_lower[m]
            factible = true
            # naturalmente si al hacer el cambio, el ĩ incumple con lower, no lo queremos hacer
        end
    end
    risk_vec[ĩ] -= R[j] # restale a ĩ el viejo
    risk_vec[i] += R[j] # sumale a i el nuevo
    if risk_vec[i] < β
        factible = true
    end
    if risk_vec[ĩ] < β
        factible = true
    end
    return factible
end

function isFactible4(solution::Types.Solution, targets_lower, targets_upper, verbose=true)
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
    P = instance.S
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

    if sum(Y) > P
        if verbose
            println("Violando número de centros asignados ", sum(Y))
        end
        number_constraints_violated += 1
    end

    for j in 1:B
        if sum(X[i, j] for i in usables_i) ≠ 1
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

struct move_bu
    ĩ::Int64
    i::Int64
    j::Int64
    cons::Int64
end

function get_cons(var::move_bu)
    return var.cons
end

# seguir agregando o quitando hasta que sea factible+

function minimums(v, n)
    l = length(v)
    ix = [1:l;]
    partialsortperm!(ix, v, (l-n+1):l, rev=true, initialized=true)
    vals = reverse!(v[ix[(l-n+1):l]])
    indices = reverse!(CartesianIndices(v)[ix[(l-n+1):l]])
    return indices
end

function maximums(vec, n)::Tuple{Vector{Int64}, Vector{CartesianIndex{1}}}
    type = eltype(vec)
    vals = zeros(type, n)
    indices = Array{Int64}(undef,n)
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

function move_bu_repair(solution, cons, targets_lower, targets_upper, remove, add)
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
    n = round(Int, (P/2))
    not_usables_i = Set(findall(==(0), Y))
    values_matrix = Matrix{Int64}(undef, S, M)
    risk_vec = Vector{Int64}(undef, S)

    values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
    for i in add
        prev_assigned_bus = findall(==(1), X[i,:]) #indices de las bus asignadas previamente
        Dᵢ = D[i,:]
        for prev in prev_assigned_bus
            Dᵢ[prev] = 1e9
        end
        # las previamente asignadas las ponemos muy grandes, buscamos los minimoss
        candidates_bus = minimums(Dᵢ, n) #indices de las bus mas cercanas al centro i
        factible_yet = false
        while !factible_yet
            j = popfirst!(candidates_bus)
            ĩ = findfirst(==(1), X[:,j]) # obten la asignacion previa de cada j
            factible_yet = true
            can_do_move = true
            for m in 1:M
                values_matrix[ĩ,m] -= V[m][j] # restale a ĩ, previo
                values_matrix[i,m] += V[m][j] # sumale a i, nuevo
                if values_matrix[i,m] > targets_upper[m]
                    # no deberia de pasar porque entonces la infactibilidad cambia de razon
                    factible_yet = false
                    can_do_move = false
                end
                if values_matrix[i,m] < targets_lower[m]
                    factible_yet = false
                end
                if values_matrix[ĩ,m] < targets_lower[m]
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
                    values_matrix[ĩ,m] += V[m][j] # corrige el valor de la NO asignacion
                    values_matrix[i,m] -= V[m][j]
                end
                risk_vec[ĩ] += R[j]
                risk_vec[i] -= R[j]
            end
        end
    end
    
    values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
    for ĩ in remove
        values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
        # aqui el chiste es quitar las BUs ASIGNADAS MAS lejanas
        # para asegurar que solo son las asignadas, multiplicamos la D por X
        all_asigneds = findall(==(1), X[ĩ, :])
        Dᵢ = copy(D[ĩ,:]) .* X[ĩ,:] # al hacer 0s los no asignados, no los agarra maximums
        vals, candidates_bus = maximums(Dᵢ, length(all_asigneds))
        factible_yet = false
        while !factible_yet
            factible_yet = true
            j = popfirst!(candidates_bus)[1]
            # cuales is podemos agarrar? las n mas cercanas? probamos todas? 
            # agarramos las n iₛ mas cercanas a j
            distance_to_bu = copy(D[:,j])
            for not_usable_i in not_usables_i
                distance_to_bu[not_usable_i] = 1e9 # a las is en Yi = 0, hacemos la distancia super grande para no fomentar su uso
            end
            j_assigned = false
            candidates_is = minimums(distance_to_bu, n)
            while !j_assigned
                i = popfirst!(candidates_is)
                i = i[1]
                can_do_move = true
                if i == ĩ
                    can_do_move = false
                end
                for m in 1:M
                    values_matrix[ĩ,m] -= V[m][j] # restale a ĩ, previo
                    values_matrix[i,m] += V[m][j] # sumale a i, nuevo
                    if values_matrix[i,m] > targets_upper[m]
                        factible_yet = false
                    end
                    if values_matrix[i,m] < targets_lower[m]
                        factible_yet = false
                        can_do_move = false
                        # no deberia de pasar porque la infactibilidad cambia de razon entonces
                    end
                    if values_matrix[ĩ,m] < targets_lower[m]
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
                if risk_vec[ĩ] > B
                    factible_yet = false
                end
                if can_do_move
                    X[ĩ, j] = 0
                    X[i, j] = 1
                    j_assigned = true
                else
                    # si no puedo hacer el movimiento, restaura el valor de la ev parcial
                    for m in 1:M
                        values_matrix[ĩ,m] += V[m][j] # corrige el valor de la NO asignacion
                        values_matrix[i,m] -= V[m][j]
                    end
                    risk_vec[ĩ] += R[j]
                    risk_vec[i] -= R[j]
                end
            end
        end
    end
    Weight = 0
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += D[indice]
    end
    newSol = Solution(instance, X, Y, Weight, solution.Time)
    println(isFactible(solution, false))
    println(isFactible(newSol, false))
    return newSol
end

function move_bu_improve(solution, strategy=:bf)
end

function localSearch(solPath, plotPath, solution::Types.Solution)
    oldSol = solution
    oldTime = oldSol.Time
    instance = solution.Instance
    targets_lower, targets_upper = calculate_targets(instance)
    factible, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper)
    original_weight = 1e9
    if factible
        println("Factible")
        original_weight = solution.Weight
    else
        println("Infactible")
    end

    
end