include("constructive.jl")
include("local.jl")
using Types
using Dates
using DelimitedFiles
using .Threads

function grasp(α, iters, instance)
    start = now()
    bestSol = nothing
    bestWeight = 1e12
    instancia = instance
    B = instancia.B
    S = instancia.S
    P = instancia.P
    X = Matrix{Int64}
    Y = Vector{Int64}
    D = instancia.D
    Weight = 0
    iter_value = Matrix{Int64}(undef, iters, 2)

    Y, time_init = relax_init(instance)
    println(time_init)
    facilities_on = findall(x -> x == 1, Y)
    facilities_off = findall(x -> x == 0, Y)
    modifications = round(Int, α * length(Y))
    lockVar = ReentrantLock()
    Threads.@threads for i in 1:iters
        println("iter: ", i)
        before_assignment = now()
        println("Entrando a opp cost assignment desde: $(threadid())")
        X = oppCostAssignmentGrasp(Y, instance, α, facilities_on, facilities_off, modifications)
        after_assignment = now()
        delta = after_assignment - before_assignment
        time_secs = round(delta, Second)
        println(delta)
        peso = evalWeight(X, D)
        total_time = time_secs.value + time_init
        sol = Solution(instance, X, Y, peso, total_time)
        factible, cons_v_og = isFactible(sol, false)
        weight_before = 1e12
        if factible
            weight_before = peso
        end
        println(cons_v_og)
        before_ls = Dates.now()
        solution, cons_v = simple_move_bu(sol, cons_v_og)
        #println("Constraints despues de simple: ", cons_v)
        mejora_cons_sim = cons_v_og - cons_v
        #println("Mejora constraints simple: ", mejora_cons_sim)
        mejora_peso_sim = 0
        cons_v_og = cons_v
        if cons_v_og == 0
            valor_previo = solution.Weight
            if factible
                mejora_peso_sim = weight_before - valor_previo
            end
        end
        println("Mejora valor simple: ", mejora_peso_sim, " thread: $(threadid())")
        #println("Simple terminado con tiempo: ", delta_simple_milli)
        solution, cons_v = interchange_bus(solution, cons_v)
        #println("Constraints despues de intercambio: ", cons_v)
        mejora_cons_int = cons_v_og - cons_v
        #println("Mejora constraints intercambio: ", mejora_cons_int)
        mejora_peso_int = 0
        cons_v_og = cons_v
        if cons_v_og == 0
            mejora_peso_int = valor_previo - solution.Weight
            valor_previo = solution.Weight
        end
        println("Mejora valor intercambio: ", mejora_peso_int, " thread: $(threadid())")
        after_ls = Dates.now()
        delta_ls = after_ls - before_ls
        time_secs_ls = round(delta_ls, Second)
        peso = solution.Weight
        println("Bloqueando ")
        lock(lockVar)
        try
            if peso < bestWeight
                println("Mejor valor que el incumbente: ", peso)
                bestWeight = peso
                total_time = time_secs.value + time_init + time_secs_ls.value
                bestSol = Solution(instance, X, Y, peso, total_time)
            end
            iter_value[i, 1] = i
            iter_value[i, 2] = bestWeight
        finally
            println("Desbloqueando ")
            unlock(lockVar)
        end
        println("peso de: ", peso)
        println("tiempo de: ", total_time)
    end

    open("delim_iters_1200.txt", "w") do io
        writedlm(io, iter_value)
    end

    finish = now()

    println(finish - start)

    return bestSol, bestWeight
end

function oppCostAssignmentGrasp(Y, instance::Types.Instance, α, facilities_on, facilities_off, modifications)
    D = copy(instance.D)
    P = instance.P
    S = instance.S
    B = instance.B
    X = zeros(Int64, S, B)
    modified_Y = copy(Y)
    to_off = Int64[]
    to_on = Int64[]
    while length(to_off) < modifications
        off = facilities_on[rand(1:end)]
        if off ∉ to_off
            push!(to_off, off)
        end
    end
    while length(to_on) < modifications
        on = facilities_off[rand(1:end)]
        if on ∉ to_on
            push!(to_on, on)
        end
    end
    for idx in to_on
        modified_Y[idx] = 1
    end
    for idx in to_off
        modified_Y[idx] = 0
    end
    count = 0
    not_assigned_y = findall(y -> y == 0, Y)
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
    n = trunc(Int, P - 10)
    while !todos
        values, indices::Vector{CartesianIndex{2}} = maximums(diff, n)
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
            bestVal = values[1]
            cutoff = trunc(Int, (bestVal + (bestVal * α))) # todas las posibles asignaciones que esten dentro del rango de la mejor posible
            rcl = [x for x in values if x <= cutoff]
            idx_rcl = [x for x in eachindex(values) if values[x] <= cutoff]
            idx_random = idx_rcl[rand(1:end)] # escoge uno al azar de la rcl
            indice = indices[idx_random]::CartesianIndex{2} 
            col = indice[2]::Int64
            row = findfirst(x -> x == 0, diff[:, col])
            picked = CartesianIndex(row, col)
        else
            if all(x -> x ≠ 0, constraints) # si todas violan constraints
                original_cons, idx = findmin(constraints) # agarra el que viole menos
                cutoff_cons = trunc(Int, (original_cons + (original_cons * α)))
                rcl_inner = [x for x in constraints if x <= cutoff_cons]
                idx_rcl_inner = [x for x in eachindex(constraints) if constraints[x] <= cutoff_cons]
                idx_random = idx_rcl_inner[rand(1:end)] # escoge uno al azar de la rcl
                indice = indices[idx_random] # de nuevo escoge una al azar
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

function evalWeight(X, D)
    Weight = 0
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += D[indice]
    end
    return Weight
end

function main_grasp()
    file_name = ARGS[1]
    instance = read_instance(file_name)
    α = 0.2
    iters = 10
    solucion, peso = grasp(α, iters, instance)
    write_solution(solucion, "solucion_grasp_multithread_1_310.jld2")
end

main_grasp()