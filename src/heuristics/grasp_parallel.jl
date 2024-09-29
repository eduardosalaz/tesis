include("constructive.jl")
include("ls.jl")
using Types
using Dates
using DelimitedFiles
using TimerOutputs
using .Threads

function grasp(αₗ, αₐ, max_iters, instance)
    start = now()
    bestSol = nothing
    bestWeight = 100000000
    instancia = instance
    B = instancia.B
    S = instancia.S
    P = instancia.P
    D = instancia.D
    targets_lower, targets_upper = calculate_targets(instance)
    count_repair_1 = 0
    count_repair_2 = 0
    lockVar = ReentrantLock()
    Threads.@threads for _ in 1:max_iters
        X = Matrix{Int64}(undef, S, B)
        Y = Vector{Int64}(undef, S)
        Weight = 0
        start_iter = now()
        #println(iter)
        Y, time_loc = pdisp_grasp(instance, αₗ)
        X, time_alloc = oppCostQueueGRASP(Y, instance, αₐ)
        Weight = 0
        indices = findall(x -> x == 1, X)
        for indice in indices
            Weight += D[indice]
        end
        #@show Weight
        oldSol = Types.Solution(instance, X, Y, Weight, time_loc + time_alloc)
        #println(isFactible(oldSol))
        repair_delta = 0
        factible_after_repair = false
        factible, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper, false)
        if !factible
            repaired_1 = repair_solution1(oldSol, constraints, targets_lower, targets_upper, remove, add)
            fac_repaired_1, cons = isFactible(repaired_1, false)
            if !fac_repaired_1
                repaired_2 = repair_solution2(oldSol, constraints, targets_lower, targets_upper, remove, add)
                fac_repaired_2, cons = isFactible(repaired_2, false)
                if fac_repaired_2
                    repair_algorithm = 2
                    count_repair_2 += 1
                    factible_after_repair = true
                    repaired = repaired_2
                end
            else
                count_repair_1 += 1
                repaired = repaired_1
                factible_after_repair = true
            end
            if factible_after_repair
                original_weight = repaired.Weight
                weight_before = repaired.Weight
                #println("Reparada")
            end
        end
        #println(original_weight)
        #println(iter, " ", factible_after_repair)
        if factible_after_repair
            oldSol = repaired
            #println(isFactible(oldSol, true))
        end
        improvement = true
        while improvement && factible_after_repair
            improvement = false  # Reset the flag at the start of each loop iteration
            prev_weight = oldSol.Weight

            # Array to keep track of individual improvements
            improvements = Bool[]

            # First improvement function
            sol_moved_bu = simple_bu_improve(oldSol, targets_lower, targets_upper, :ff)
            new_weight_moved = sol_moved_bu.Weight
            push!(improvements, new_weight_moved < prev_weight)
            if improvements[end]
                prev_weight = new_weight_moved
                #println("En el loop $loop el movimiento simple mejora con un $new_weight_moved")
                oldSol = sol_moved_bu  # Update oldSol if there was an improvement
            end
            #println(isFactible(sol_moved_bu, true))
            # Second improvement function
            sol_interchanged_bu = interchange_bu_improve(oldSol, targets_lower, targets_upper, :ff)
            new_weight_moved = sol_interchanged_bu.Weight
            push!(improvements, new_weight_moved < prev_weight)
            if improvements[end]
                prev_weight = new_weight_moved
                #println("En el loop loop el movimiento intercambio mejora con un new_weight_moved")
                oldSol = sol_interchanged_bu  # Update oldSol if there was an improvement
            end
            #println(isFactible(sol_interchanged_bu, true))
            # Third improvement function
            
            sol_deactivated_center = deactivate_center_improve(oldSol, targets_lower, targets_upper)
            new_weight_moved = sol_deactivated_center.Weight
            push!(improvements, new_weight_moved < prev_weight)
            if improvements[end]
                prev_weight = new_weight_moved
                #println("En el loop loop el movimiento desactivar mejora con un new_weight_moved")
                oldSol = sol_deactivated_center  # Update oldSol if there was an improvement
            end
            
            #println(isFactible(sol_deactivated_center, true))
            # Check for any improvements

            improvement = any(improvements)
            #println(isFactible(oldSol, true))
        end
        end_iter = now()
        #@show end_iter
        delta_total = end_iter - start
        delta_total_millis = round(delta_total, Millisecond)
        if delta_total_millis.value >= 1800000 # 30 minutos en milisegundos
            break
        end
        #delta_iter = end_iter - start_iter
        #delta_iter_millis = round(delta_iter, Millisecond)
        #delta_iter_value = delta_iter_millis.value

        lock(lockVar)
        try
            current_timestamp = Dates.datetime2unix(now(Dates.UTC))
            if oldSol !== nothing
                #=
                if !factible_after_repair
                    results[iter] = (iter=iter, weight=0, time=delta_iter_value, factible=false, improving=false, curr_epoch=current_timestamp)
                else
                    =#
                    improving = false
                    if oldSol.Weight < bestWeight
                        bestSol = oldSol
                        bestWeight = oldSol.Weight
                        improving = true
                    end
                    #results[iter] = (iter=iter, weight=oldSol.Weight, time=delta_iter_value, factible=true, improving=improving, curr_epoch=current_timestamp)
                #end
            #else
                #if !factible_after_repair
                    #results[iter] = (iter=iter, weight=0, time=delta_iter_value, factible=false, improving=false, curr_epoch=current_timestamp)
                #end
            end
        finally
            unlock(lockVar)
        end
        #@show oldSol.Weight

        #gap_repaired = (1 - (weight_exac / weight_before)) * 100
        #gap_improved = (1 - (weight_exac / weight_after)) * 100
        #abs_improved = weight_before - weight_after
        #rel_improved = ((abs_improved) / weight_after) * 100
    end
    finish = now()
    # @show count_repair_1
    #@show count_repair_2
    delta = finish - start
    delta_millis = round(delta, Millisecond)
    # println(delta_millis.value)
    return bestSol, delta_millis.value # , results
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


function calculate_targets_upper(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    targets = Vector{Float32}(undef, M)
    for m in 1:M
        targets[m] = (1 * μ[m][1] * (1 + T[m])) # corregir esto
    end
    return targets
end

function calculate_targets_lower(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    targets = Vector{Float32}(undef, M)
    for m in 1:M
        targets[m] = (1 * μ[m][1] * (1 - T[m])) # corregir esto
    end
    return targets
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

function start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    for m in eachindex(V)
        mul!(view(values_matrix, :, m), X, V[m])
    end
    mul!(risk_vec, X, R)
    return values_matrix, risk_vec
end

function update_constraints(M, V, R, i, j, values_matrix, risk_vec, targets, β)
    constraints = 0
    for m in 1:M
        values_matrix[i, m] += V[m][j]
        if values_matrix[i, m] > targets[m]
            constraints += 1
        end
    end
    risk_vec[i] += R[j]
    if risk_vec[i] > β
        constraints += 1
    end
    return constraints
end

function compute_assignments_and_opportunity_costs(D::Matrix{Int64}, Y::Vector{Int}, N::Int)
    _, num_clients = size(D)
    max_facilities = length(Y)
    best_assignments = Dict{Int,Vector{Int64}}()
    pq = PriorityQueue{Int,Int64}(Base.Order.Reverse)

    for j in 1:num_clients
        # Use a temporary array to store the facility opportunity costs for this client
        costs = Tuple{Int64,Int}[]
        for (i, yi) in enumerate(Y)
            if yi == 1
                push!(costs, (D[i, j], i))
            end
        end
        # Sort the costs
        sort!(costs)
        # Store the top N assignments for this client
        best_assignments[j] = best_assignments[j] = [cost[2] for cost in costs[1:N]] # extrae el indice nada mas
        # Calculate the opportunity cost as the largest difference among all possible assignments
        opp_cost = costs[end][1] - costs[1][1]  # Difference between largest and smallest costs
        # Use the smallest facility index for this client to form the unique key
        # unique_key = encode_key(j, costs[1][2], max_facilities)
        enqueue!(pq, j, opp_cost)
    end
    return best_assignments, pq
end

function pick_center_from_rcl(rcl, full_centers)
    feasible_rcl = setdiff(rcl, full_centers)
    if isempty(feasible_rcl)
        return nothing  # Or some sentinel value indicating no feasible center found
    end
    return rand(feasible_rcl)
end

function oppCostQueueGRASP(Y, instance::Types.Instance, α)
    before_alloc = Dates.now()
    D = copy(instance.D)
    X = zeros(Int, size(D))
    targets = calculate_targets_lower(instance)
    β = instance.β[1]
    S = instance.S
    B = instance.B
    M = instance.M
    V = instance.V
    P = instance.P
    R = instance.R
    values_matrix = Matrix{Float32}(undef, S, M)
    risk_vec = Vector{Float32}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    #values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, zeros(S, M), zeros(S))
    N = instance.P # podemos hacerlo en proporcion a P, no necesariamente tiene que ser P
    best_assignments, pq = compute_assignments_and_opportunity_costs(D, Y, N)
    full_centers = Set()
    assigned_bus = Set()
    unassigned_bus = Set(collect(1:B))
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        if bu ∉ assigned_bus
            best_value = first(best_assignments[bu])[1]
            threshold = best_value * (1 + α)
            RCL = [facility for facility in best_assignments[bu] if facility[1] <= threshold]
            center = pick_center_from_rcl(RCL, full_centers)
            if center !== nothing
                full = false
                fulls_m = zeros(Int, M)
                if center ∉ full_centers
                    for m in 1:M
                        values_matrix[center, m] += V[m][bu]
                        if values_matrix[center, m] > targets[m]
                            fulls_m[m] = 1
                        end
                    end
                    if all(x -> x == 1, fulls_m)
                        push!(full_centers, center)
                    end
                    risk_vec[center] += R[bu]
                    if risk_vec[center] > β
                        push!(full_centers, center)
                    end
                    push!(assigned_bus, bu)
                    pop!(unassigned_bus, bu)
                    X[center, bu] = 1
                    break
                end
            end
        end
    end
    todos = true
    count = 0
    for col in eachcol(X)
        if all(x -> x == 0, col)
            count += 1
            todos = false
        end
    end
    X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    after_alloc = Dates.now()
    delta_alloc = after_alloc - before_alloc
    delta_alloc_milli = round(delta_alloc, Millisecond)
    return X, delta_alloc_milli.value
end

function handle_unassigned_clients2(X, instance, best_assignments, unassigned_clients, values_matrix, risk_vec)
    targets_upper = calculate_targets_upper(instance)
    S = instance.S
    M = instance.M
    V = instance.V
    β = instance.β[1]
    R = instance.R

    # Only loop over unassigned clients
    for client in unassigned_clients
        assigned = false

        # Sort facilities for this client based on their remaining risk capacity
        sorted_best_assignments = sort(best_assignments[client], by=f -> β - risk_vec[f], rev=true)

        for facility in sorted_best_assignments
            potential_assignment_valid = true
            for m in 1:M
                if values_matrix[facility, m] + V[m][client] > targets_upper[m]
                    potential_assignment_valid = false
                    break
                end
            end

            if risk_vec[facility] + R[client] > β
                potential_assignment_valid = false
            end

            if potential_assignment_valid
                X[facility, client] = 1
                assigned = true
                delete!(unassigned_clients, client)
                for m in 1:M
                    values_matrix[facility, m] += V[m][client]
                end
                risk_vec[facility] += R[client]
                break  # Assign to the first valid facility and then exit the inner loop
            end
        end

        if !assigned
            @error "Client $client could not be assigned."
        end
    end

    # println("Number of unassigned clients: $(length(unassigned_clients))")
    return X
end

function update_center_capacity(center, bu, values_matrix, risk_vec, targets, β, M, V, R)
    for m in 1:M
        values_matrix[center, m] += V[m][bu]
        if values_matrix[center, m] > targets[m]
            return true
        end
    end
    risk_vec[center] += R[center]
    if risk_vec[center] > β
        return true
    end
    return falses
end

function pdisp_simple_grasp(d, p, N, α)
    maxdist = 0
    bestpair = (0, 1)
    for i in 1:N
        for j in i+1:N
            if d[i, j] > maxdist
                maxdist = d[i, j]
                bestpair = (i, j)
            end
        end
    end
    P = Set([])
    push!(P, bestpair[1])
    push!(P, bestpair[2])

    while length(P) < p
        maxdist = 0
        vbest = 0
        minimal = Inf
        for v in 1:N
            if v in P
                continue
            end
            mindist = Inf
            for vprime in P
                if d[v, vprime] < mindist
                    mindist = d[v, vprime]
                end
            end
            if mindist > maxdist
                maxdist = mindist
                minimal = mindist
            end
        end
        maxdist = minimal
        vbest = 0
        rcl = []
        for v in 1:N
            if v in P
                continue
            end
            mindist = Inf
            for vprime in P
                if d[v, vprime] < mindist
                    mindist = d[v, vprime]
                end
            end
            if mindist >= round(Int, (maxdist - maxdist * α))
                push!(rcl, v)
            end
        end
        selected = false
        while !selected
            v_candidate = rand(rcl)
            if v_candidate != 0 && !(v_candidate in P)
                push!(P, v_candidate)
                selected = true
            end
        end
    end
    collection = collect(P)
    return collection
end

function pdisp_grasp(instance, αₗ)
    before_init = Dates.now()
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
    k = 5
    s_coords = instance.S_coords
    metric = Distances.Euclidean()
    d = trunc.(Int, Distances.pairwise(metric, s_coords, dims=1))
    coords_S1 = s_coords[Sk[1], :]
    coords_S2 = s_coords[Sk[2], :]
    coords_S3 = s_coords[Sk[3], :]
    coords_S4 = s_coords[Sk[4], :]
    coords_S5 = s_coords[Sk[5], :]

    d1 = trunc.(Int, Distances.pairwise(metric, coords_S1, dims=1))
    d2 = trunc.(Int, Distances.pairwise(metric, coords_S2, dims=1))
    d3 = trunc.(Int, Distances.pairwise(metric, coords_S3, dims=1))
    d4 = trunc.(Int, Distances.pairwise(metric, coords_S4, dims=1))
    d5 = trunc.(Int, Distances.pairwise(metric, coords_S5, dims=1))
    N1 = length(Sk[1])
    N2 = length(Sk[2])
    N3 = length(Sk[3])
    N4 = length(Sk[4])
    N5 = length(Sk[5])
    p1 = Lk[1]
    p2 = Lk[2]
    p3 = Lk[3]
    p4 = Lk[4]
    p5 = Lk[5]

    pdisp1 = pdisp_simple_grasp(d1, p1, N1, αₗ)
    pdisp2 = pdisp_simple_grasp(d2, p2, N2, αₗ)
    pdisp3 = pdisp_simple_grasp(d3, p3, N3, αₗ)
    pdisp4 = pdisp_simple_grasp(d4, p4, N4, αₗ)
    pdisp5 = pdisp_simple_grasp(d5, p5, N5, αₗ)

    pdisp1_fixed = Sk[1][pdisp1]
    pdisp2_fixed = Sk[2][pdisp2]
    pdisp3_fixed = Sk[3][pdisp3]
    pdisp4_fixed = Sk[4][pdisp4]
    pdisp5_fixed = Sk[5][pdisp5]

    N = S
    pdisp_ok = Set(vcat([pdisp1_fixed, pdisp2_fixed, pdisp3_fixed, pdisp4_fixed, pdisp5_fixed]...))
    if length(pdisp_ok) != P
        count = count_k(pdisp_ok, Sk)
        while length(pdisp_ok) < P
            # Find the node v that maximizes the distance to its closest neighbor in P
            maxdist = 0
            vbest = 0
            best_dist = 0
            for v in 1:N
                if v in pdisp_ok
                    continue
                end
                k = node_type(v, Sk)
                if count[k] >= Uk[k]
                    continue
                end
                dist = minimum([d[v, vprime] for vprime in pdisp_ok])
                if dist > maxdist
                    maxdist = dist
                    best_dist = dist
                    #vbest = v
                end
            end
            maxdist = 0
            vbest = 0
            rcl = []
            #println(best_dist)
            #println(round(Int, best_dist - best_dist * αₗ))
            for v in 1:N
                if v in pdisp_ok
                    continue
                end
                k = node_type(v, Sk)
                if count[k] >= Uk[k]
                    continue
                end
                dist = minimum([d[v, vprime] for vprime in pdisp_ok])
                if dist >= round(Int, best_dist - best_dist * αₗ)
                    push!(rcl, v)
                    #vbest = v
                end
            end
            vbest = rand(rcl)
            #println("Escogiendo $vbest")
            # If no such node exists, stop the algorithm
            if vbest == 0
                @error "P DISP ERROR"
                break
            end
            # Add the node vbest to the set P and update the counts
            k = node_type(vbest, Sk)
            count[k] += 1
            push!(pdisp_ok, vbest)
        end
    end
    collection = collect(pdisp_ok)
    after_init = Dates.now()
    delta_init = after_init - before_init
    delta_init_milli = round(delta_init, Millisecond)
    Y_bool = zeros(Int, instance.S)
    for idx in collection
        Y_bool[idx] = 1
    end
    return Y_bool, delta_init_milli.value
end


function main_grasp(;path="solucion_grasp_16_625_feas.jld2", iters=10)
    #file_name = "instances\\625_78_32\\inst_1_625_78_32.jld2"
    instance = read_instance(path)
    αₗ = 0.3
    αₐ = 0.3
    iters = parse(Int, ARGS[2])
    bestSolution, totalTime = grasp(αₗ, αₐ, iters, instance)
    println(totalTime)
    println(bestSolution.Weight)
    #println(isFactible(bestSolution))
    #=
    sorted_results_desc = sort(results, by=p -> p.curr_epoch) # al ordenar por epoch podeoms comparar contra el inmediato anterior?
    for var in sorted_results_desc
        println(var)
    end
    =#
    write_solution(bestSolution, "solucion_grasp_1_1250_changes7.jld2")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_grasp(;path=ARGS[1], iters=parse(Int, ARGS[2]))
else
    main_grasp()
end