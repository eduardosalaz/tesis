using Types

function isFactible(solution::Solution)
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
            println("Violando Lk en $k")
            number_constraints_violated += 1
        end
        if counts_k[k] > Uk[k]
            println("Violando Uk en $k")
            number_constraints_violated += 1
        end
    end

    for i in 1:S
        for m in 1:M
            if Y[i] * μ[m][i] * (1 - T[m]) > sum(X[i,j] * V[m][j] for j in 1:B)
                println("violando V inferior en i: $i y m: $m")
                println("μ: ", Y[i] * μ[m][i] * (1 - T[m]))
                println("V: ", sum(X[i,j] * V[m][j] for j in 1:B))
                number_constraints_violated += 1
            end
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= Y[i] * μ[m][i] * (1 + T[m]))
                println("violando V superior en i: $i y m: $m")
                println("μ: ", Y[i] * μ[m][i] * (1 + T[m]))
                println("V: ", sum(X[i,j] * V[m][j] for j in 1:B))
                number_constraints_violated += 1
            end
        end
    end

    for i in 1:S
        if sum(X[i, j] * R[j] for j in 1:B) > β[i]
            println("violando riesgo en $i")
            number_constraints_violated += 1
        end
    end
    println("Restricciones violadas: $number_constraints_violated")
        # esto esta mal btw
        # intercambiar nodo i por nodo j 
        # mover nodo de territorio i a territorio j
        # apagar un branch y prender otro branch
        # checar ILS ??
    if number_constraints_violated ≠ 0
        return false
    else
        return true
    end
end


function main()
    path = ARGS[1]
    solution = read_solution(path)
    factible = isFactible(solution)
    println(factible)
end

main()