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
            if !any(x->x==1, assignments_y)
                if verbose
                    println("Violando asignación de Y en: $y")
                end
                number_constraints_violated +=1
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
                    println("μ: ", Y[i] * μ[m][i] * (1 - T[m]))
                    println("V: ", sum(X[i,j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= Y[i] * μ[m][i] * (1 + T[m]))
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", Y[i] * μ[m][i] * (1 + T[m]))
                    println("V: ", sum(X[i,j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
        end
    end

    for i in 1:S
        if sum(X[i, j] * R[j] for j in 1:B) > β[i]
            if verbose
                println("violando riesgo en $i")
                println("β: ",  β[i])
                println("R: ", sum(X[i, j] * R[j] for j in 1:B))
            end
            number_constraints_violated += 1
        end
    end
    println("Restricciones violadas: $number_constraints_violated")
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

function localSearch(solution::Solution)
    factible, cons_v = isFactible(solution, false)
    move = :simple_move_bu
    if !factible
        if move == :simple_move_bu
            newSol = simple_move_bu(solution, cons_v)
        end
    end
    return newSol
end

function simple_move_bu(solution::Solution, cons_v)
    # este movimiento se define de la siguiente manera:
    # Muevo una unidad básica ψ de un territorio i a un territorio τ
    # La unidad básica originalmente se define como ψ = X[i,j] = 1
    # El movimiento consiste en iterar para cada j dentro de B:
    #   X[τ, j] = 1, donde τ es distinto de i
    # Se evalua cada solucion buscando la minimización del número de constraints violadas
    original_cons = cons_v
    original_sol = solution
    instance = solution.Instance
    best_ψiτ = (1,1,1)
    B = instance.B
    S = instance.S
    max_steps = 10
    step = 1
    while original_cons ≠ 0
        X = solution.X
        for j in 1:B
            ψ = j
            original_b = findall(x->x==1, X[:,j])[1] # asignacion inicial de la sucursal
            for i in 1:S
                if i ≠ original_b
                    new_cons, newSol = move_bu(solution, original_b, i, j)
                    a = 1
                    if new_cons < original_cons
                        if new_cons == 0
                            println("Solucion factible, no hay constraints violadas")
                        end
                        original_cons = new_cons
                        τ = i
                        best_ψiτ = (ψ, original_b, τ)
                        solution = newSol
                        println("new_cons: $new_cons")
                        println("pichula")
                        println("best tuple: $best_ψiτ")
                    else
                        break
                    end
                end
                step += 1
            end
        end
        if step > max_steps
            break
        end
        println(step)
    end
    return solution
end


function evalWeight(X, D)
    Weight = 0
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += D[indice]
    end
    return Weight
end


function move_bu(solution::Solution, original_b, i, j)
    # movimiento simple de basic unit
    X = solution.X
    Y = solution.Y
    instance = solution.Instance
    D = instance.D
    diff_X = X
    diff_X[original_b, j] = 0
    diff_X[i, j] = 1
    Weight = evalWeight(diff_X, D)
    newSol = Solution(instance, diff_X, Y, Weight)
    factible, cons_v = isFactible(newSol, false)
    if factible
        return 0, newSol
    else
        return cons_v, newSol
    end

end


function main()
    # path = ARGS[1]
    path = "sol_34_pdisp.jld2"
    solution = read_solution(path)
    newSol = localSearch(solution)
    write_solution(newSol, "sol_34_ls3_pdisp.jld2")
    plot_solution(newSol, "plot_sol_34_ls3_pdisp.png")
end

main()