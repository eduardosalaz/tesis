using Types
using Dates

function isFactible(solution::Solution, verbose=true)
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
            if !(trunc(Int, Y[i] * μ[m][i] * (1 - T[m])) <= sum(X[i, j] * V[m][j] for j in 1:B))
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

function assigntounused(solution::Types.Solution) # en esta funcion agarramos las Y[i] = 1 pero que no tienen asignacion y hacemos una asignacion simple
    instance = solution.Instance
    Y = solution.Y
    X = solution.X
    X2 = copy(X)
    D = instance.D
    for y in eachindex(Y)
        if Y[y] == 1
            assignments_y = X2[y, :]
            if !any(x -> x == 1, assignments_y)
                nearest_bu = argmin(D[y, :])
                X2[:, nearest_bu] .= 0
                X2[y, nearest_bu] = 1 # asignamos a la bu mas cercana para no dejar sola la branch
            end
        end
    end
    weight = evalWeight(X2, D)
    newSol = Types.Solution(instance, X2, Y, weight)
    return newSol
end

function localSearch(solPath, plotPath, solution::Solution)
    oldSol = solution
    oldTime = oldSol.Time
    factible, cons_v = isFactible(solution, false)
    weight_before = 1e9
    if factible
        println("Solucion Factible")
        weight_before = solution.Weight
    else
        println("Solucion no factible")
    end
    before = now()
    solution, cons_v = simple_move_bu(solution, cons_v)
    println("Simple terminado")
    solution, cons_v = interchange_bus(solution, cons_v)
    println("Intercambio terminado")
    solution, cons_v = deactivateBranch(solution, cons_v)
    after = now()
    delta = after - before
    secs = round(delta, Second)
    time = secs.value
    println("Desactivar y Activar terminado")
    println("Fin de Busqueda Local")
    factible, cons_v = isFactible(solution, false)
    if factible
        println("Solucion factible despues de busqueda local")
    else
        println("Bad luck")
    end
    newSol = solution
    weight_now = newSol.Weight
    X_now = newSol.X
    Y_now = newSol.Y
    instance = newSol.instance
    newTime = oldTime + time

    lsSol = Solution(instance, X_now, Y_now, weight_now, newTime)
    @show weight_now
    
    Types.write_solution(lsSol, solPath)
    Types.plot_solution(lsSol, plotPath)
    return newSol
end


struct interchangeMove
    newConstraints::Int64
    newSolution::Types.Solution
    bu_1::Int64
    bu_2::Int64
    i::Int64
    j::Int64
    weight::Int64
end

function interchange_bus(solution::Types.Solution, cons_v;) # en este movimiento se intercambian dos BUs 1 y 2 de territorios i y j
    original_solution = solution
    original_X = original_solution.X
    original_Y = original_solution.Y
    original_w = original_solution.Weight
    original_cons = cons_v
    instance = original_solution.Instance
    factible = false
    if original_cons == 0
        factible = true
    end
    S = instance.S
    B = instance.B
    D = instance.D
    max_iters = 300
    iter = 0
    while !factible
        if iter > max_iters
            return solution, original_cons
        end
        Xs = []
        for j in 1:B
            moves = interchangeMove[]
            i_original = findall(x -> x == 1, original_X[:, j])
            if length(i_original) > 0
                i_original = i_original[1]
            else
                continue
            end
            for j2 in 1:B
                iter += 1
                if iter > max_iters
                    return solution, original_cons
                end
                if j ≠ j2
                    j_original = findall(x -> x == 1, original_X[:, j2])
                    if length(j_original) > 0
                        j_original = j_original[1]
                    else
                        continue
                    end
                    X_copy = copy(original_X)
                    X_copy[i_original, j] = 0
                    X_copy[i_original, j2] = 1
                    X_copy[j_original, j2] = 0
                    X_copy[j_original, j] = 1
                    if X_copy ∉ Xs # sin duplicados
                        push!(Xs, X_copy)
                        Y_copy = copy(original_Y)
                        for i in 1:S
                            if any(entry -> entry == 1, X_copy[i, :]) # revisando cada fila si hay un 1, si hay entonces Y[i] se usa
                                Y_copy[i] = 1
                            else # si no hay ni un solo 1, entonces Y[i] = 0
                                Y_copy[i] = 0
                            end
                        end
                        Weight = evalWeight(X_copy, D)
                        newSol = Types.Solution(instance, X_copy, Y_copy, Weight, original_solution.Time)
                        factible, cons_v = isFactible(newSol, false)
                        interchangedMove = interchangeMove(cons_v, newSol, j, j2, i_original, j_original, Weight)
                        push!(moves, interchangedMove)
                    end
                end
            end
            sort!(moves, by=move -> move.newConstraints) # ordena en base a las menores constraints violadas
            if length(moves) == 0
                return solution, original_cons
            else
                bestMove = moves[1] # obten el mejor movimiento
                if bestMove.newConstraints < original_cons
                    original_X = bestMove.newSolution.X # actualizamos la X
                    original_cons = bestMove.newConstraints
                    solution = bestMove.newSolution
                    println(bestMove.newConstraints)
                    if bestMove.newConstraints == 0 # empezamos a minimizar entonces el valor de la funcion objetivo
                        @info "Factible solution"
                        factible = true
                    end
                end
            end
        end
    end
    original_w = solution.Weight
    iter = 1
    for iter in 1:max_iters
        Xs = []
        for j in 1:B
            moves = interchangeMove[]
            i_original = findall(x -> x == 1, original_X[:, j])
            if length(i_original) > 0
                i_original = i_original[1]
            else
                continue
            end
            for j2 in 1:B
                iter += 1
                if iter > max_iters
                    return solution, original_cons
                end
                if j ≠ j2
                    j_original = findall(x -> x == 1, original_X[:, j2])
                    if length(j_original) > 0
                        j_original = j_original[1]
                    else
                        continue
                    end
                    X_copy = copy(original_X)
                    X_copy[i_original, j] = 0
                    X_copy[i_original, j2] = 1
                    X_copy[j_original, j2] = 0
                    X_copy[j_original, j] = 1
                    if X_copy ∉ Xs # sin duplicados
                        push!(Xs, X_copy)
                        Y_copy = copy(original_Y)
                        for i in 1:S
                            if any(entry -> entry == 1, X_copy[i, :]) # revisando cada fila si hay un 1, si hay entonces Y[i] se usa
                                Y_copy[i] = 1
                            else # si no hay ni un solo 1, entonces Y[i] = 0
                                Y_copy[i] = 0
                            end
                        end
                        Weight = evalWeight(X_copy, D)
                        newSol = Types.Solution(instance, X_copy, Y_copy, Weight, original_solution.Time)
                        factible, cons_v = isFactible(newSol, false)
                        interchangedMove = interchangeMove(cons_v, newSol, j, j2, i_original, j_original, Weight)
                        push!(moves, interchangedMove)
                    end
                end
            end
            sort!(moves, by=move -> move.weight) # ordena en base a las menores constraints violadas
            if length(moves) == 0
                return solution, original_cons
            else
                canMove = false
                while !canMove
                    bestMove = moves[1] # obten el mejor movimiento
                    canMove, _ = isFactible(bestMove.newSolution, false)
                    if bestMove.weight < original_w && canMove
                        println(bestMove.weight)
                        println(original_w)
                        w = bestMove.weight
                        @info "Pasando a una mejor solucion con peso:" w
                        original_X = bestMove.newSolution.X # actualizamos la X
                        original_w = bestMove.weight
                        solution = bestMove.newSolution
                    else
                        popfirst!(moves) # si no es factible el movimiento, quitalo y busca el siguiente
                    end
                end  
            end
        end
    end
    return solution, original_cons
end

struct deactivatedMove
    newConstraints::Int64
    newSolution::Types.Solution
    deactivated_branch::Int64
    activated_branch::Int64
    weight::Int64
end

# en este movimiento hay que apagar una branch completa y prender otra
function deactivateBranch(solution::Types.Solution, cons_v)
    original_solution = solution
    original_X = original_solution.X
    original_Y = original_solution.Y
    Y = copy(original_Y)
    X = copy(original_X)
    original_w = original_solution.Weight
    original_cons = cons_v
    instance = original_solution.Instance
    factible = false
    S = instance.S
    B = instance.B
    D = instance.D
    max_iters = 300
    iter = 0
    while !factible
        if iter > max_iters
            break
        end
        I = findall(y -> y == 1, Y)
        for i in I
            moves = []
            J = findall(y -> y == 0, Y)
            bus = findall(x -> x == 1, X[i, :])
            for j in J
                iter += 1
                Y2 = copy(Y)
                Y2[i] = 0
                Y2[j] = 1
                X2 = copy(X)
                X2[i, :] .= 0
                for bu in bus
                    # asignar a la nueva branch
                    X2[j, bu] = 1
                end
                Weight = evalWeight(X2, D)
                newSol = Types.Solution(instance, X2, Y2, Weight, original_solution.Time)
                factible, cons_v = isFactible(newSol, false)
                deactivatedmove = deactivatedMove(cons_v, newSol, i, j, Weight)
                push!(moves, deactivatedmove)
            end
            sort!(moves, by=move -> move.newConstraints) # ordena en base a las menores constraints violadas
            if length(moves) == 0
                return solution, original_cons
            else
                bestMove = moves[1] # obten el mejor movimiento
                if bestMove.newConstraints < original_cons
                    X = bestMove.newSolution.X # actualizamos la X
                    Y = bestMove.newSolution.Y
                    original_cons = bestMove.newConstraints
                    solution = bestMove.newSolution
                    println(bestMove.newConstraints)
                    if bestMove.newConstraints == 0 # empezamos a minimizar entonces el valor de la funcion objetivo
                        @info "Factible solution"
                        factible = true
                    end
                end
            end
        end
        if !factible
            return solution, cons_v
        end
        iter = 1
        for iter in 1:max_iters
            I = findall(y -> y == 1, Y)
            for i in I
                moves = []
                J = findall(y -> y == 0, Y)
                bus = findall(x -> x == 1, X[i, :])
                for j in J
                    iter += 1
                    Y2 = copy(Y)
                    Y2[i] = 0
                    Y2[j] = 1
                    X2 = copy(X)
                    X2[i, :] .= 0
                    for bu in bus
                        # asignar a la nueva branch
                        X2[j, bu] = 1
                    end
                    Weight = evalWeight(X2, D)
                    newSol = Types.Solution(instance, X2, Y2, Weight, original_solution.Time)
                    factible, cons_v = isFactible(newSol, false)
                    deactivatedmove = deactivatedMove(cons_v, newSol, i, j, Weight)
                    push!(moves, deactivatedmove)
                end
                sort!(moves, by=move -> move.weight) # ordena en base al peso
                if length(moves) == 0
                    return solution, original_cons
                else
                    canMove = false
                    while !canMove
                        bestMove = moves[1] # obten el mejor movimiento
                        canMove, _ = isFactible(bestMove.newSolution, false)
                        if bestMove.weight < original_w && canMove
                            w = bestMove.weight
                            @info "Pasando a una mejor solucion con peso:" w
                            original_X = bestMove.newSolution.X # actualizamos la X
                            original_w = bestMove.weight
                            solution = bestMove.newSolution
                        else
                            popfirst!(moves) # si no es factible el movimiento, quitalo y busca el siguiente
                        end
                    end  
                end
            end
        end
    end
    return solution, original_cons
end

struct simpleMove
    newConstraints::Int64
    newSolution::Types.Solution
    original_branch::Int64
    new_branch::Int64
    bu::Int64
    weight::Int64
end


function simple_move_bu(solution::Types.Solution, cons_v)
    # este movimiento se define de la siguiente manera:
    # Muevo una unidad básica ψ de un territorio i a un territorio τ
    # La unidad básica originalmente se define como ψ = X[i,j] = 1
    # El movimiento consiste en iterar para cada j dentro de B:
    #   X[τ, j] = 1, donde τ es distinto de i
    # Se evalua cada solucion buscando la minimización del número de constraints violadas
    original_solution = solution
    original_X = original_solution.X
    original_w = original_solution.Weight
    original_cons = cons_v
    instance = original_solution.Instance
    factible = false
    if original_cons == 0
        factible = true
    end
    S = instance.S
    B = instance.B
    max_iters = 250
    iter = 0
    while !factible
        println("No es factible")
        if iter > max_iters
            break
        end
        for j in 1:B
            println(iter)            
            if iter > max_iters
                break
            end
            moves = simpleMove[]
            i_original = findall(x -> x == 1, original_X[:, j])[1] # el nodo asociado inicialmente a la BU
            for i in 1:S
                if i ≠ i_original
                    new_cons_v, new_sol = move_bu(solution, j, i_original, i)
                    move = simpleMove(new_cons_v, new_sol, i_original, i, j, new_sol.Weight)
                    push!(moves, move)
                end
            end
            sort!(moves, by=move -> move.newConstraints) # ordena en base a las menores constraints violadas
            bestMove = moves[1] # obten el mejor movimiento
            if bestMove.newConstraints < original_cons
                original_X = bestMove.newSolution.X # actualizamos la X
                original_cons = bestMove.newConstraints
                solution = bestMove.newSolution
                println(bestMove.newConstraints)
                if bestMove.newConstraints == 0 # empezamos a minimizar entonces el valor de la funcion objetivo
                    @info "Factible solution"
                    factible = true
                end
            end
        end
    end
    if !factible
        return solution, original_cons
    end
    iter = 1
    for iter in 1:max_iters
        for j in 1:B
            moves = simpleMove[]
            i_original = findall(x -> x == 1, original_X[:, j])[1] # el nodo asociado inicialmente a la BU
            for i in 1:S
                if i ≠ i_original
                    iter += 1
                    new_cons_v, new_sol = move_bu(solution, j, i_original, i)
                    move = simpleMove(new_cons_v, new_sol, i_original, i, j, new_sol.Weight)
                    push!(moves, move)
                end
            end
            sort!(moves, by=move -> move.weight) # ordena en base al peso
            canMove = false
            if length(moves) == 0
                canMove = true
            end
            while !canMove
                if length(moves) == 0
                    canMove = true
                    return solution, original_cons
                end
                bestMove = moves[1] # obten el mejor movimiento
                canMove, _ = isFactible(bestMove.newSolution, false)
                if bestMove.weight < original_w && canMove
                    w = bestMove.weight
                    @info "Pasando a una mejor solucion con peso:" w
                    original_X = bestMove.newSolution.X # actualizamos la X
                    original_w = bestMove.weight
                    solution = bestMove.newSolution
                else
                    popfirst!(moves) # si no es factible el movimiento, quitalo y busca el siguiente
                end
            end
        end
    end
    return solution, original_cons
end

function evalWeight(X, D)
    Weight = 0
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += D[indice]
    end
    return Weight
end

function move_bu(solution::Types.Solution, j, previous, new)
    # movimiento simple de basic unit
    X = solution.X
    Y = solution.Y
    instance = solution.Instance
    S, B = size(X)
    D = instance.D
    X_copy = copy(X)
    Y_copy = copy(Y)
    X_copy[previous, j] = 0
    X_copy[new, j] = 1
    Weight = evalWeight(X_copy, D)
    # ahora hay que revisar que Y no haya cambiado
    for i in 1:S
        if any(entry -> entry == 1, X_copy[i, :]) # revisando cada fila si hay un 1, si hay entonces Y[i] se usa
            Y_copy[i] = 1
        else # si no hay ni un solo 1, entonces Y[i] = 0
            Y_copy[i] = 0
        end
    end
    newSol = Types.Solution(instance, X_copy, Y_copy, Weight, solution.Time)
    factible, cons_v = isFactible(newSol, false)
    if factible
        return 0, newSol
    else
        return cons_v, newSol
    end
end


function main_local(; path="inst", read_file=true, sol_obj=nothing)
    cons_solution = 1 # scope
    newSolPath = ""
    newPlotPath = ""
    if read_file
        cons_solution = read_solution(path)
    else
        cons_solution = sol_obj
    end
    newSolPath = replace(path, ".jld2" => "_ls.jld2")
    newPlotPath = replace(newSolPath, ".jld2" => ".png")
    newSol = localSearch(newSolPath, newPlotPath, cons_solution)
    return newSol
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_local(; path=ARGS[1])
end