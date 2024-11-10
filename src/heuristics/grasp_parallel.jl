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
    targets_lower_op, targets_upper_op = calculate_targets_optimized(instance)
    count_repair_1 = 0
    count_repair_2 = 0
    lockVar = ReentrantLock()
    Threads.@threads for _ in 1:max_iters
        X = Matrix{Int64}(undef, S, B)
        Y = Vector{Int64}(undef, S)
        Weight = 0
        start_iter = now()
        #println(iter)
        Y, time_loc = pdisp_grasp_new(instance, αₗ)
        X, time_alloc = oppCostQueueGRASPnew(Y, instance, αₐ)
        Weight = 0
        indices = findall(x -> x == 1, X)
        for indice in indices
            Weight += D[indice]
        end
        oldSol = Types.Solution(instance, X, Y, Weight, time_loc + time_alloc)
        repair_delta = 0
        factible_after_repair = false
        println("fac 39 $factible_after_repair")
        factible, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper, false)
        if !factible
            println("nao nao")
            repaired_1 = repair_solution1(oldSol, constraints, targets_lower, targets_upper, remove, add)
            fac_repaired_1, cons = isFactible(repaired_1, false)
            println("cons1 $cons")
            println("fac 46 $factible_after_repair")
            if !fac_repaired_1
                repaired_2 = repair_solution2(oldSol, constraints, targets_lower, targets_upper, remove, add)
                fac_repaired_2, cons = isFactible(repaired_2, false)
                println("cons2 $cons")
                println("fac 51 $factible_after_repair")
                if fac_repaired_2
                    repair_algorithm = 2
                    count_repair_2 += 1
                    factible_after_repair = true
                    repaired = repaired_2
                    println("fac 57 $factible_after_repair")
                end
            else
                count_repair_1 += 1
                repaired = repaired_1
                factible_after_repair = true
                println("fac 62 $factible_after_repair")
            end
            if factible_after_repair
                println("yey")
                original_weight = repaired.Weight
                weight_before = repaired.Weight
            end
        else
            repaired = oldSol
            factible_after_repair = true
            println("fac 73 $factible_after_repair")
        end
        println("fac 75 $factible_after_repair")
        if factible_after_repair
            oldSol = repaired
        end
        improvement = true
        last_weights = Dict(:simple_bu => Inf, :interchange_bu => Inf, :deactivate_center => Inf)
        println(improvement && factible_after_repair)
        println(bestSol)
        while (improvement && factible_after_repair)
            
            println("fac 83 $factible_after_repair")
            println("entre al loop")
            improvement = false
            prev_weight = oldSol.Weight
            
            # Simple BU move
            if oldSol.Weight != last_weights[:simple_bu]
                sol_moved_bu = simple_bu_improve_optimized(oldSol, targets_lower_op, targets_upper_op, :ff)
                if sol_moved_bu.Weight < prev_weight
                    oldSol = sol_moved_bu
                    prev_weight = sol_moved_bu.Weight
                    improvement = true
                end
                last_weights[:simple_bu] = oldSol.Weight
            end
            
            # Interchange BU move
            if oldSol.Weight != last_weights[:interchange_bu]
                sol_interchanged_bu = interchange_bu_improve_optimized(oldSol, targets_lower_op, targets_upper_op, :ff)
                if sol_interchanged_bu.Weight < prev_weight
                    oldSol = sol_interchanged_bu
                    prev_weight = sol_interchanged_bu.Weight
                    improvement = true
                end
                last_weights[:interchange_bu] = oldSol.Weight
            end
            
            # Deactivate center move
            if oldSol.Weight != last_weights[:deactivate_center]
                sol_deactivated_center = deactivate_center_improve(oldSol, targets_lower, targets_upper)
                if sol_deactivated_center.Weight < prev_weight
                    oldSol = sol_deactivated_center
                    prev_weight = sol_deactivated_center.Weight
                    improvement = true
                end
                last_weights[:deactivate_center] = oldSol.Weight
            end
        end

        end_iter = now()
        delta_total = end_iter - start
        delta_total_millis = round(delta_total, Millisecond)
        if delta_total_millis.value >= 1800000 # 30 minutos en milisegundos
            break
        end
        if factible_after_repair
            lock(lockVar)
            try
                if oldSol !== nothing
                        improving = false
                        if oldSol.Weight < bestWeight
                            println("updateando")
                            bestSol = oldSol
                            bestWeight = oldSol.Weight
                            improving = true
                        end
                end
            finally
                unlock(lockVar)
            end
        end
        
    end
    println(bestSol)
    finish = now()
    delta = finish - start
    delta_millis = round(delta, Millisecond)
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

    # O(B)
    for j in 1:num_clients
        # Use a temporary array to store the facility opportunity costs for this client
        # O(Y)
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


function oppCostQueueGRASPnew(Y, instance::Types.Instance, α)
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
            # Get the sorted list of assignments for this bus
            assignments = best_assignments[bu]
            
            # Calculate ϕ values for all non-full centers
            ϕ = OrderedDict()
            for center in best_assignments[bu]
                if center ∉ full_centers
                    ϕ[center] = D[center, bu]  # Direct distance calculation
                end
            end
            
            # If all centers are full, continue to the next bus
            if isempty(ϕ)
                continue
            end
            
            # ϕ_min and ϕ_max are now the first and last values in the ordered dictionary
            ϕ_values = collect(values(ϕ))
            ϕ_min = ϕ_values[1]
            ϕ_max = ϕ_values[end]
            
            # Define RCL based on actual distances for minimization
            RCL = [center for (center, dist) in ϕ if dist <= ϕ_min + α * (ϕ_max - ϕ_min)]
            
            # Select a random center from RCL
            if !isempty(RCL)
                center = rand(RCL)
                
                # Update values_matrix and check targets
                fulls_m = zeros(Int, M)
                for m in 1:M
                    values_matrix[center, m] += V[m][bu]
                    if values_matrix[center, m] > targets[m]
                        fulls_m[m] = 1
                    end
                end

                # Check if center is full for all target types
                if all(x -> x == 1, fulls_m)
                    push!(full_centers, center)
                end

                # Update risk and check risk threshold
                risk_vec[center] += R[bu]
                if risk_vec[center] > β
                    push!(full_centers, center)
                end

                # Assign bus to center
                push!(assigned_bus, bu)
                pop!(unassigned_bus, bu)
                X[center, bu] = 1
            end
        end
    end
    # Check if all buses are assigned
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
    # @show length(unassigned_clients)

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
            if mindist >= maxdist - α * (maxdist - mindist) # kms
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

function pdisp_simple_grasp_new(d, p, N, α)
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
    
    P = Set([bestpair[1], bestpair[2]])
    
    while length(P) < p
        ϕ = Dict()
        ϕ_max = -Inf
        ϕ_min = Inf
        
        for v in 1:N
            if v in P
                continue
            end
            ϕ[v] = minimum(d[v, vprime] for vprime in P)
            ϕ_max = max(ϕ_max, ϕ[v])
            ϕ_min = min(ϕ_min, ϕ[v])
        end
        
        rcl = [v for v in keys(ϕ) if ϕ[v] >= ϕ_max - α * (ϕ_max - ϕ_min)]
        #println("ϕ_max: $(ϕ_max)")
        if !isempty(rcl)
            v_candidate = rand(rcl)
            push!(P, v_candidate)
            #println("pushing to P $v_candidate")
        else
            #println("empty")
            # If RCL is empty, choose the point with the maximum ϕ value
            v_candidate = argmax(ϕ)
            push!(P, v_candidate)
        end
    end
    
    return collect(P)
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

function pdisp_grasp_new(instance, αₗ)
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

    pdisp1 = pdisp_simple_grasp_new(d1, p1, N1, αₗ)
    pdisp2 = pdisp_simple_grasp_new(d2, p2, N2, αₗ)
    pdisp3 = pdisp_simple_grasp_new(d3, p3, N3, αₗ)
    pdisp4 = pdisp_simple_grasp_new(d4, p4, N4, αₗ)
    pdisp5 = pdisp_simple_grasp_new(d5, p5, N5, αₗ)

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
            ϕ = Dict()
            ϕ_max = -Inf
            ϕ_min = Inf
            
            # Calculate ϕ(v) for all eligible nodes
            for v in 1:N
                if v in pdisp_ok
                    continue
                end
                k = node_type(v, Sk)
                if count[k] >= Uk[k]
                    continue
                end
                ϕ[v] = minimum(d[v, vprime] for vprime in pdisp_ok)
                ϕ_max = max(ϕ_max, ϕ[v])
                ϕ_min = min(ϕ_min, ϕ[v])
            end
            
            # Create the RCL based on the new GRASP logic
            rcl = [v for v in keys(ϕ) if ϕ[v] >= ϕ_max - αₗ * (ϕ_max - ϕ_min)]
            
            # If RCL is empty, choose the point with the maximum ϕ value
            if isempty(rcl)
                vbest = argmax(ϕ)
            else
                vbest = rand(rcl)
            end
            
            # If no eligible node exists, stop the algorithm
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
    pattern = Regex("[t][_]\\d{1,3}")
    index = findfirst(pattern, path)
    almost_number = path[index]
    _, number = split(almost_number, "_")
    αₗ = 0.0
    αₐ = 0.0    
    iters = parse(Int, ARGS[2])
    bestSolution, totalTime = grasp(αₗ, αₐ, iters, instance)
    println(totalTime)
    println(bestSolution.Weight)
    println(isFactible(bestSolution))
    name = "solucion_grasp_newinstance_$number" * "_$αₗ" * "_$αₐ" *"_$iters" * "_$(nthreads())" * ".jld2"
    full_out_path = "out\\solutions\\1250_155_62\\grasp\\" * name
    write_solution(bestSolution, full_out_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_grasp(;path=ARGS[1], iters=parse(Int, ARGS[2]))
else
    #main_grasp()
end