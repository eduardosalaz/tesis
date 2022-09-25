using Types

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
            if !(Y[i] * μ[m][i] * (1 - T[m]) <= sum(X[i, j] * V[m][j] for j in 1:B))
                if verbose
                    println("violando V inferior en i: $i y m: $m")
                    println("μ: ", Y[i] * μ[m][i] * (1 - T[m]))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= Y[i] * μ[m][i] * (1 + T[m]))    
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", Y[i] * μ[m][i] * (1 + T[m]))
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

function assigntounused(solution::Solution) # en esta funcion agarramos las Y[i] = 1 pero que no tienen asignacion y hacemos una asignacion simple
    number_constraints_violated = 0
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
    newSol = Solution(instance, X2, Y, weight)
    return newSol
end

function localSearch(solution::Solution)
    factible, cons_v = isFactible(solution, false)
    println(cons_v)
    println("______________________________________________________")
    solution2 = assigntounused(solution)
    factibleAfter, cons_v_after = isFactible(solution2, false)
    move = :simple_move_bu
    if !factibleAfter
        if move == :simple_move_bu
            newSol, cons_v_after = simple_move_bu(solution2, cons_v_after)
            println("simple done ", cons_v_after)
            newSol, cons_v_after = interchange_bus(newSol, cons_v_after)
            println("interchange done ", cons_v_after)
            newSol, cons_v_after = deactivateBranch(newSol, cons_v_after)
            println("daectivate done ", cons_v_after)
            newSol, cons_v_after = simple_move_bu(newSol, cons_v_after)
            println("simple done 2 ", cons_v_after)
            newSol, cons_v_after = interchange_bus(newSol, cons_v_after)
            println("interchange done 2 ", cons_v_after)
            newSol, cons_v_after = deactivateBranch(newSol, cons_v_after)
            println("deactivate done 2 ", cons_v_after)
            isFactible(newSol, true)
            return newSol
        end
    end
    
end


struct interchangeMove
    newConstraints::Int64
    newSolution::Solution
    bu_1::Int64
    bu_2::Int64
    i::Int64
    j::Int64
end

function interchange_bus(solution::Solution, cons_v) # en este movimiento se intercambian dos BUs 1 y 2 de territorios i y j

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
            if iter > max_iters
                return solution, original_cons
            end
            moves = interchangeMove[]
            i_original = findall(x->x ==1, original_X[:,j])
            if length(i_original)> 0
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
                    j_original = findall(x->x ==1, original_X[:,j2])
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
                        newSol = Solution(instance, X_copy, Y_copy, Weight)
                        factible, cons_v = isFactible(newSol, false)
                        interchangedMove = interchangeMove(cons_v, newSol, j, j2, i_original, j_original)
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
                    if bestMove.newConstraints == 0 # empezamos a minimizar entonces el valor de la funcion objetivo
                        @info "Factible solution"
                        factible = true
                    end
                end
            end
        end
    end
    return solution, original_cons
end

struct deactivatedMove
    newConstraints::Int64
    newSolution::Solution
    deactivated_branch::Int64
    activated_branch::Int64
end


# en este movimiento hay que apagar una branch completa y prender otra
function deactivateBranch(solution::Solution, cons_v)
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
        I = findall(y->y==1, Y)
        for i in I
            moves = []
            J = findall(y-> y == 0, Y)
            bus = findall(x->x==1, X[i,:])
            for j in J
                iter += 1
                Y2 = copy(Y)
                Y2[i] = 0
                Y2[j] = 1
                X2 = copy(X)
                X2[i,:] .= 0
                for bu in bus
                    # asignar a la nueva branch
                    X2[j, bu] = 1
                end
                #=
                
                =#
                Weight = evalWeight(X2, D)
                newSol = Solution(instance, X2, Y2, Weight)
                factible, cons_v = isFactible(newSol, false)
                deactivatedmove = deactivatedMove(cons_v, newSol, i, j)
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
                    if bestMove.newConstraints == 0 # empezamos a minimizar entonces el valor de la funcion objetivo
                        @info "Factible solution"
                        factible = true
                    end
                end
            end
        end
    end
    return solution, original_cons
end

struct simpleMove
    newConstraints::Int64
    newSolution::Solution
    original_branch::Int64
    new_branch::Int64
    bu::Int64
end


function simple_move_bu(solution::Solution, cons_v)
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
    S = instance.S
    B = instance.B
    max_iters = 300
    iter = 0
    while !factible
        if iter > max_iters
            break
        end
        for j in 1:B
            moves = simpleMove[]
            i_original = findall(x -> x == 1, original_X[:, j])[1] # el nodo asociado inicialmente a la BU
            for i in 1:S
                if i ≠ i_original
                    iter += 1
                    new_cons_v, new_sol = move_bu(solution, j, i_original, i)
                    move = simpleMove(new_cons_v, new_sol, i_original, i, j)
                    push!(moves, move)
                end
            end
            sort!(moves, by=move -> move.newConstraints) # ordena en base a las menores constraints violadas
            bestMove = moves[1] # obten el mejor movimiento
            if bestMove.newConstraints < original_cons
                original_X = bestMove.newSolution.X # actualizamos la X
                original_cons = bestMove.newConstraints
                solution = bestMove.newSolution
                if bestMove.newConstraints == 0 # empezamos a minimizar entonces el valor de la funcion objetivo
                    @info "Factible solution"
                    factible = true
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

#=
for j in 1:B
           i_original = findall(x->x==1, X[:,j])[1]
           for i in 1:S
               if i ≠ i_original
                   X_copy = copy(X)
                   X_copy[i_original, j] = 0
                   X_copy[i, j] = 1
                   count += 1
                   if count < max_count
                       println(X_copy)
                       println()
                       println(X)
                   end
               end
           end
       end
=#

#=
contador = 0
for j in 1:B
    i_original = findall(x->x==1, X[:,j])[1]
    @show j, i_original
    for i in 1:S
        if i ≠ i_original
            X_copy = copy(X)
            X_copy[i_original, j] = 0
            X_copy[i, j] = 1
            Y_copy = copy(Y)
            contador += 1
            for aux in 1:S
                if any(entry->entry==1, X_copy[aux,:]) # revisando cada fila si hay un 1, si hay entonces Y[i] se usa
                    Y_copy[aux] = 1
                else # si no hay ni un solo 1, entonces Y[i] = 0
                    Y_copy[aux] = 0
                end
            end
        end
    end
end
=#


function move_bu(solution::Solution, j, previous, new)
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
    newSol = Solution(instance, X_copy, Y_copy, Weight)
    factible, cons_v = isFactible(newSol, false)
    if factible
        return 0, newSol
    else
        return cons_v, newSol
    end
end

function interchange()
    
end


function main()
    path = ARGS[1]
    #path = "sol_34_pdisp.jld2"
    solution = read_solution(path)
    newSol = localSearch(solution)
    write_solution(newSol, "sol6_1_200_ls_pdisp2.jld2")
    plot_solution(newSol, "plot6_sol_1_200_ls_pdisp2.png")
end

main()