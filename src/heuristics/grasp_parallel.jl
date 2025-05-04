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
        
        X, time_alloc = grasp_allocation_with_queue2(Y, instance, αₐ)
        Weight = 0
        indices = findall(x -> x == 1, X)
        for indice in indices
            Weight += D[indice]
        end
        oldSol = Types.Solution(instance, X, Y, Weight, time_loc + time_alloc)
        #println(isFactible(oldSol))
        repair_delta = 0
        factible_after_repair = false
        #println(isFactible(oldSol, false))
        factible, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper, false)
        if !factible
            repaired_1 = repair_solution1(oldSol, constraints, targets_lower, targets_upper, remove, add)
            fac_repaired_1, cons = isFactible(repaired_1, false)
            if !fac_repaired_1
                #repaired_2 = repair_solution2(oldSol, constraints, targets_lower, targets_upper, remove, add)
                #fac_repaired_2, cons = isFactible(repaired_2, false)
                #if fac_repaired_2
                #    repair_algorithm = 2
                #    count_repair_2 += 1
                #    factible_after_repair = true
                #    repaired = repaired_2
                #end
            else
                count_repair_1 += 1
                repaired = repaired_1
                factible_after_repair = true
            end
            if factible_after_repair
                original_weight = repaired.Weight
                weight_before = repaired.Weight
            end
        else
            repaired = oldSol
            factible_after_repair = true
        end
        if factible_after_repair
            oldSol = repaired
        end
        if factible_after_repair
            println("yes")
            println("old weight: $(threadid()) ,", oldSol.Weight)
            #println(isFactible(oldSol))
            oldSol = vnd_local_search(oldSol, targets_lower_op, targets_upper_op, targets_lower, targets_upper)
            println("new weight: $(threadid()) ,", oldSol.Weight)
            #println(isFactible(oldSol))
        end
        

        end_iter = now()
        delta_total = end_iter - start
        delta_total_millis = round(delta_total, Millisecond)
        if delta_total_millis.value >= 1800000 # 30 minutos en milisegundos
            break
        end
        if factible_after_repair
            #println("FACTIBLE")
            lock(lockVar)
            try
                if oldSol !== nothing
                        improving = false
                        if oldSol.Weight < bestWeight
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
    finish = now()
    delta = finish - start
    delta_millis = round(delta, Millisecond)
    return bestSol, delta_millis.value # , results
end

function calculate_targets(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    S = instance.S
    targets_lower = Matrix{Int}(undef, S, M)
    targets_upper = Matrix{Int}(undef, S, M)
    for s in 1:S
        for m in 1:M
            targets_lower[s, m] = ceil(Int, (1 * μ[m][s] * (1 - T[m])))
            targets_upper[s, m] = floor(Int, (1 * μ[m][s] * (1 + T[m])))
        end
    end
    return targets_lower, targets_upper
end


function calculate_targets_upper(instance)
    M = instance.M
    S = instance.S
    μ = instance.μ
    T = instance.T
    # Now we need targets for each center and each activity type
    targets = Matrix{Int}(undef, S, M)
    for s in 1:S
        for m in 1:M
            targets[s, m] = floor(Int, μ[m][s] * (1 + T[m]))
        end
    end
    return targets
end

function calculate_targets_lower(instance)
    M = instance.M
    S = instance.S
    μ = instance.μ
    T = instance.T
    targets = Matrix{Int}(undef, S, M)
    for s in 1:S
        for m in 1:M
            targets[s, m] = ceil(Int, μ[m][s] * (1 - T[m]))
        end
    end
    return targets
end

function calculate_targets_optimized(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    S = instance.S
    
    # Create matrices to store lower and upper bounds for each center and activity
    targets_lower = Matrix{Int}(undef, S, M)
    targets_upper = Matrix{Int}(undef, S, M)
    
    # Calculate bounds for each center i and activity m
    for i in 1:S
        for m in 1:M
            targets_lower[i,m] = ceil(Int, (μ[m][i] * (1 - T[m])))
            targets_upper[i,m] = floor(Int, (μ[m][i] * (1 + T[m])))
        end
    end
    
    # Convert to static matrices
    targets_lower = SMatrix{S,M,Int}(targets_lower)
    targets_upper = SMatrix{S,M,Int}(targets_upper)
    
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

# no
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
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
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

# no
function graspOppCostQueue(Y, instance::Types.Instance, alpha::Float64)
    D = copy(instance.D)
    X = zeros(Int, size(D))
    targets = calculate_targets_lower(instance)
    targets_upper = calculate_targets_upper(instance)
    β = instance.β[1]
    S = instance.S
    B = instance.B
    M = instance.M
    V = instance.V
    P = instance.P
    R = instance.R
    
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    
    N = instance.P
    best_assignments, pq = compute_assignments_and_opportunity_costs(D, Y, N)
    
    # Track progress towards targets for each facility
    target_progress = zeros(Float32, S, M)
    for i in 1:S
        for m in 1:M
            target_progress[i,m] = values_matrix[i,m] / targets[m]
        end
    end
    
    full_centers = Set()
    assigned_bus = Set()
    unassigned_bus = Set(collect(1:B))
    
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        if bu ∉ assigned_bus
            # Get all possible centers for this BU
            possible_centers = best_assignments[bu]
            
            # Create RCL based on distance
            center_distances = [(center, D[center, bu]) for center in possible_centers if center ∉ full_centers]
            if isempty(center_distances)
                continue
            end
            
            # Sort by distance
            sort!(center_distances, by = x -> x[2])
            
            # Calculate RCL bounds
            min_dist = center_distances[1][2]
            max_dist = center_distances[end][2]
            threshold = min_dist + alpha * (max_dist - min_dist)
            
            # Create RCL
            rcl = [center for (center, dist) in center_distances if dist ≤ threshold]
            
            if isempty(rcl)
                continue
            end
            
            # Evaluate centers in RCL based on improvement in target progress
            center_improvements = Float64[]
            valid_centers = Int[]
            
            for center in rcl
                # Calculate improvement in target progress
                current_min_progress = minimum(target_progress[center,:])
                
                # Simulate adding this BU
                temp_progress = copy(target_progress[center,:])
                for m in 1:M
                    temp_progress[m] = (values_matrix[center,m] + V[m][bu]) / targets[m]
                end
                
                new_min_progress = minimum(temp_progress)
                improvement = new_min_progress - current_min_progress
                
                # Check if this assignment would exceed upper bounds
                exceeds_upper = false
                for m in 1:M
                    if values_matrix[center,m] + V[m][bu] > targets_upper[m]
                        exceeds_upper = true
                        break
                    end
                end
                
                # Check risk constraint
                if !exceeds_upper && risk_vec[center] + R[bu] <= β
                    push!(center_improvements, improvement)
                    push!(valid_centers, center)
                end
            end
            
            if !isempty(valid_centers)
                # Select randomly from the best improvements
                max_improvement = maximum(center_improvements)
                min_improvement = minimum(center_improvements)
                improvement_threshold = max_improvement - alpha * (max_improvement - min_improvement)
                
                candidates = [(c, imp) for (c, imp) in zip(valid_centers, center_improvements) 
                            if imp ≥ improvement_threshold]
                
                # Randomly select from candidates
                selected_idx = rand(1:length(candidates))
                best_center = candidates[selected_idx][1]
                
                # Make the assignment
                for m in 1:M
                    values_matrix[best_center,m] += V[m][bu]
                    target_progress[best_center,m] = values_matrix[best_center,m] / targets[m]
                end
                risk_vec[best_center] += R[bu]
                
                # Check if center should be marked as full
                min_progress = minimum(target_progress[best_center,:])
                if min_progress >= 0.95  # Within 5% of targets
                    push!(full_centers, best_center)
                end
                
                push!(assigned_bus, bu)
                X[best_center, bu] = 1
                pop!(unassigned_bus, bu)
            end
        end
    end
    
    unassigned_count = length(unassigned_bus)
    if unassigned_count > 0
        X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    end
    
    return X, unassigned_count
end

# no
function graspOppCostQueue2(Y, instance::Types.Instance, alpha::Float64)
    D = copy(instance.D)
    X = zeros(Int, size(D))
    targets = calculate_targets_lower(instance)
    targets_upper = calculate_targets_upper(instance)
    β = instance.β[1]
    S = instance.S
    B = instance.B
    M = instance.M
    V = instance.V
    P = instance.P
    R = instance.R
    
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    
    N = instance.P
    best_assignments, pq = compute_assignments_and_opportunity_costs(D, Y, N)
    
    # Track progress towards targets for each facility
    target_progress = zeros(Float32, S, M)
    for i in 1:S
        for m in 1:M
            target_progress[i,m] = values_matrix[i,m] / targets[m]
        end
    end
    
    full_centers = Set()
    assigned_bus = Set()
    unassigned_bus = Set(collect(1:B))
    
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        if bu ∉ assigned_bus
            # Get distances to all possible centers
            center_distances = [(center, D[center, bu]) for center in best_assignments[bu] 
                              if center ∉ full_centers]
            
            if !isempty(center_distances)
                # Create RCL based purely on distance
                sort!(center_distances, by = x -> x[2])  # Sort by distance
                min_dist = center_distances[1][2]
                max_dist = center_distances[end][2]
                threshold = min_dist + alpha * (max_dist - min_dist)
                
                # RCL contains centers within distance threshold
                rcl = [center for (center, dist) in center_distances if dist ≤ threshold]
                
                # Find best improvement among RCL centers
                best_center = nothing
                best_improvement = -Inf
                
                for center in rcl
                    # Calculate improvement in target progress
                    current_min_progress = minimum(target_progress[center,:])
                    temp_progress = copy(target_progress[center,:])
                    
                    for m in 1:M
                        temp_progress[m] = (values_matrix[center,m] + V[m][bu]) / targets[m]
                    end
                    
                    new_min_progress = minimum(temp_progress)
                    improvement = new_min_progress - current_min_progress
                    
                    # Check constraints
                    exceeds_upper = false
                    for m in 1:M
                        if values_matrix[center,m] + V[m][bu] > targets_upper[m]
                            exceeds_upper = true
                            break
                        end
                    end
                    
                    # Update best if feasible and better
                    if !exceeds_upper && risk_vec[center] + R[bu] <= β && improvement > best_improvement
                        best_improvement = improvement
                        best_center = center
                    end
                end
                
                # Make assignment if we found a valid center
                if best_center !== nothing
                    for m in 1:M
                        values_matrix[best_center,m] += V[m][bu]
                        target_progress[best_center,m] = values_matrix[best_center,m] / targets[m]
                    end
                    risk_vec[best_center] += R[bu]
                    
                    # Check if center should be marked as full
                    min_progress = minimum(target_progress[best_center,:])
                    if min_progress >= 0.95  # Within 5% of targets
                        push!(full_centers, best_center)
                    end
                    
                    push!(assigned_bus, bu)
                    X[best_center, bu] = 1
                    pop!(unassigned_bus, bu)
                end
            end
        end
    end
    
    unassigned_count = length(unassigned_bus)
    if unassigned_count > 0
        X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    end
    
    return X, unassigned_count
end
# no
function graspOppCostQueueCombined(Y, instance::Types.Instance, alpha::Float64)
    D = copy(instance.D)
    X = zeros(Int, size(D))
    targets = calculate_targets_lower(instance)
    targets_upper = calculate_targets_upper(instance)
    β = instance.β[1]
    S = instance.S
    B = instance.B
    M = instance.M
    V = instance.V
    P = instance.P
    R = instance.R
    
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    
    N = instance.P
    best_assignments, pq = compute_assignments_and_opportunity_costs(D, Y, N)
    
    # Track progress towards targets for each facility
    target_progress = zeros(Float32, S, M)
    for i in 1:S
        for m in 1:M
            target_progress[i,m] = values_matrix[i,m] / targets[m]
        end
    end
    
    full_centers = Set()
    assigned_bus = Set()
    unassigned_bus = Set(collect(1:B))
    
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        if bu ∉ assigned_bus
            # Get all valid centers with their distances and improvements
            candidates = Tuple{Int, Float64, Float64}[]  # (center, distance, improvement)
            
            for center in best_assignments[bu]
                if center ∉ full_centers
                    # Get distance
                    dist = D[center, bu]
                    
                    # Calculate improvement in target progress
                    current_min_progress = minimum(target_progress[center,:])
                    temp_progress = copy(target_progress[center,:])
                    
                    for m in 1:M
                        temp_progress[m] = (values_matrix[center,m] + V[m][bu]) / targets[m]
                    end
                    
                    new_min_progress = minimum(temp_progress)
                    improvement = new_min_progress - current_min_progress
                    
                    # Check constraints
                    exceeds_upper = false
                    for m in 1:M
                        if values_matrix[center,m] + V[m][bu] > targets_upper[m]
                            exceeds_upper = true
                            break
                        end
                    end
                    
                    # Add to candidates if feasible
                    if !exceeds_upper && risk_vec[center] + R[bu] <= β
                        push!(candidates, (center, dist, improvement))
                    end
                end
            end
            
            if !isempty(candidates)
                # Find bounds for both metrics
                min_dist = minimum(x[2] for x in candidates)
                max_dist = maximum(x[2] for x in candidates)
                min_imp = minimum(x[3] for x in candidates)
                max_imp = maximum(x[3] for x in candidates)
                
                # Create RCL with candidates that are good in either metric
                rcl = Int[]
                for (center, dist, imp) in candidates
                    dist_threshold = min_dist + alpha * (max_dist - min_dist)
                    imp_threshold = max_imp - alpha * (max_imp - min_imp)
                    
                    if dist <= dist_threshold || imp >= imp_threshold
                        push!(rcl, center)
                    end
                end
                
                if !isempty(rcl)
                    # Random selection from RCL
                    selected_center = rcl[rand(1:length(rcl))]
                    
                    # Make the assignment
                    for m in 1:M
                        values_matrix[selected_center,m] += V[m][bu]
                        target_progress[selected_center,m] = values_matrix[selected_center,m] / targets[m]
                    end
                    risk_vec[selected_center] += R[bu]
                    
                    # Check if center should be marked as full
                    min_progress = minimum(target_progress[selected_center,:])
                    if min_progress >= 0.95  # Within 5% of targets
                        push!(full_centers, selected_center)
                    end
                    
                    push!(assigned_bus, bu)
                    X[selected_center, bu] = 1
                    pop!(unassigned_bus, bu)
                end
            end
        end
    end
    
    unassigned_count = length(unassigned_bus)
    if unassigned_count > 0
        X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    end
    
    return X, unassigned_count
end

# si ?
function grasp_allocation_with_queue(Y, instance::Types.Instance, alpha::Float64)
    before_alloc = Dates.now()
    D = copy(instance.D)
    X = zeros(Int, size(D))
    targets = calculate_targets_lower(instance)
    targets_upper = calculate_targets_upper(instance)
    β = instance.β
    S = instance.S
    B = instance.B
    M = instance.M
    V = instance.V
    P = instance.P
    R = instance.R
    
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    
    N = instance.P
    best_assignments, pq = compute_assignments_and_opportunity_costs(D, Y, N)
    
    # Track progress towards targets for each facility
    target_progress = zeros(Float32, S, M)
    for i in 1:S
        for m in 1:M
            target_progress[i,m] = values_matrix[i,m] / targets[i,m]
        end
    end
    
    full_centers = Set()
    assigned_bus = Set()
    
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        if bu ∉ assigned_bus
            # Create RCL based on feasible assignments
            rcl = Vector{Tuple{Int, Float64}}()
            
            

            for center in best_assignments[bu]
                if center ∉ full_centers
                    # Calculate feasibility score
                    current_min_progress = minimum(target_progress[center,:])
                    
                    temp_progress = copy(target_progress[center,:])
                    for m in 1:M
                        temp_progress[m] = (values_matrix[center,m] + V[m][bu]) / targets[center,m]
                    end
                    
                    new_min_progress = minimum(temp_progress)
                    improvement = new_min_progress - current_min_progress
                    
                    # Check feasibility conditions
                    feasible = true
                    
                    # Check upper bounds
                    for m in 1:M
                        if values_matrix[center,m] + V[m][bu] > targets_upper[center,m]
                            feasible = false
                            break
                        end
                    end
                    
                    # Check risk constraint
                    if risk_vec[center] + R[bu] > β[center]
                        feasible = false
                    end
                    
                    # If assignment is feasible, add to RCL candidates
                    if feasible
                        push!(rcl, (center, improvement))
                    end
                end
            end
            
            # If we found feasible assignments, create RCL and select from it
            if !isempty(rcl)
                # Sort by improvement score
                sort!(rcl, by = x -> x[2], rev = true)
                
                # Calculate cutoff for RCL
                best_score = rcl[1][2]
                worst_score = rcl[end][2]
                threshold = best_score - alpha * (best_score - worst_score)
                
                # Create RCL with elements above threshold
                final_rcl = [center for (center, score) in rcl if score >= threshold]

                if length(final_rcl) != 1

                #println("final rcl para $bu")
                #println(final_rcl)
                for (center, score) in rcl
                    #println("Center: $center, Score: $score")
                end
                #println("Best score: $best_score")
                #println("Threshold: $threshold")
            end
                
                
                if !isempty(final_rcl)
                    # Randomly select center from RCL
                    center = rand(final_rcl)
                    
                    # Make the assignment
                    for m in 1:M
                        values_matrix[center,m] += V[m][bu]
                        target_progress[center,m] = values_matrix[center,m] / targets[center,m]
                    end
                    risk_vec[center] += R[bu]
                    
                    # Check if center should be marked as full
                    min_progress = minimum(target_progress[center,:])
                    if min_progress >= 1.0
                        push!(full_centers, center)
                    end
                    
                    push!(assigned_bus, bu)
                    X[center, bu] = 1
                end
            end
        end
    end
    
    # Handle any remaining unassigned BUs
    unassigned_bus = setdiff(Set(1:B), assigned_bus)
    unassigned_count = length(unassigned_bus)
    
    if unassigned_count > 0
        X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    end
    after_alloc = Dates.now()
    delta_alloc = after_alloc - before_alloc
    delta_alloc_milli = round(delta_alloc, Millisecond)
    
    return X, delta_alloc_milli.value
end

function grasp_allocation_with_queue2(Y, instance::Types.Instance, alpha::Float64)
    before_alloc = Dates.now()
    D = copy(instance.D)
    X = zeros(Int, size(D))
    targets = calculate_targets_lower(instance)
    targets_upper = calculate_targets_upper(instance)
    β = instance.β
    S = instance.S
    B = instance.B
    M = instance.M
    V = instance.V
    P = instance.P
    R = instance.R
    
    values_matrix = Matrix{Int}(undef, S, M)
    risk_vec = Vector{Int}(undef, S)
    values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
    
    N = instance.P
    best_assignments, pq = compute_assignments_and_opportunity_costs(D, Y, N)
    
    # Track progress towards targets for each facility
    target_progress = zeros(Float32, S, M)
    for i in 1:S
        for m in 1:M
            target_progress[i,m] = values_matrix[i,m] / targets[i,m]
        end
    end
    
    full_centers = Set()
    assigned_bus = Set()
    unassigned_bus = Set(collect(1:B))  # Added to track unassigned BUs
    
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        if bu ∉ assigned_bus
            # Create RCL based on feasible assignments
            rcl = Vector{Tuple{Int, Float64}}()
            
            for center in best_assignments[bu]
                if center ∉ full_centers
                    # Calculate feasibility score
                    current_min_progress = minimum(target_progress[center,:])
                    
                    temp_progress = copy(target_progress[center,:])
                    for m in 1:M
                        temp_progress[m] = (values_matrix[center,m] + V[m][bu]) / targets[center,m]
                    end
                    
                    new_min_progress = minimum(temp_progress)
                    improvement = new_min_progress - current_min_progress
                    
                    # Check feasibility conditions
                    feasible = true
                    
                    # Check upper bounds
                    for m in 1:M
                        if values_matrix[center,m] + V[m][bu] > targets_upper[center,m]
                            feasible = false
                            break
                        end
                    end
                    
                    # Check risk constraint
                    if risk_vec[center] + R[bu] > β[center]
                        feasible = false
                    end
                    
                    # If assignment is feasible, add to RCL candidates
                    if feasible
                        push!(rcl, (center, improvement))
                    end
                end
            end
            
            # If we found feasible assignments, create RCL and select from it
            if !isempty(rcl)
                # Sort by improvement score
                sort!(rcl, by = x -> x[2], rev = true)
                
                # Calculate cutoff for RCL
                best_score = rcl[1][2]
                worst_score = rcl[end][2]
                threshold = best_score - alpha * (best_score - worst_score)
                
                # Create RCL with elements above threshold
                final_rcl = [center for (center, score) in rcl if score >= threshold]

                if length(final_rcl) != 1
                    for (center, score) in rcl
                    end
                end
                
                if !isempty(final_rcl)
                    # Randomly select center from RCL
                    center = rand(final_rcl)
                    
                    # Make the assignment
                    for m in 1:M
                        values_matrix[center,m] += V[m][bu]
                        target_progress[center,m] = values_matrix[center,m] / targets[center,m]
                    end
                    risk_vec[center] += R[bu]
                    
                    # NEW: Check if center should be marked as full by checking remaining unassigned BUs
                    can_accept_more = false
                    for remaining_bu in unassigned_bus
                        if remaining_bu ∉ assigned_bus && remaining_bu != bu  # Exclude current BU
                            # Check if adding this BU would violate any constraint
                            would_exceed_upper = false
                            for m in 1:M
                                if values_matrix[center,m] + V[m][remaining_bu] > targets_upper[center,m]
                                    would_exceed_upper = true
                                    break
                                end
                            end
                            
                            # If at least one remaining BU can be added without violations,
                            # the center is not full
                            if !would_exceed_upper && risk_vec[center] + R[remaining_bu] <= β[center]
                                can_accept_more = true
                                break
                            end
                        end
                    end
                    
                    # Only mark as full if NO remaining BU can be added
                    if !can_accept_more
                        push!(full_centers, center)
                    end
                    
                    push!(assigned_bus, bu)
                    pop!(unassigned_bus, bu)  # Remove from unassigned set
                    X[center, bu] = 1
                end
            end
        end
    end
    
    # Handle any remaining unassigned BUs
    unassigned_count = length(unassigned_bus)
    
    if unassigned_count > 0
        X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    end
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
    β = instance.β
    R = instance.R
    # Only loop over unassigned clients
    for client in unassigned_clients
        assigned = false

        # Sort facilities for this client based on their remaining risk capacity
        sorted_best_assignments = sort(best_assignments[client], by=f -> β[f] - risk_vec[f], rev=true)

        for facility in sorted_best_assignments
            potential_assignment_valid = true
            for m in 1:M
                if values_matrix[facility, m] + V[m][client] > targets_upper[m]
                    potential_assignment_valid = false
                    break
                end
            end

            if risk_vec[facility] + R[client] > β[facility]
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
            #@error "Client $client could not be assigned."
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
        #println("rcl para individual: ")
       # println(rcl)
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
    k = 4
    s_coords = instance.S_coords
    metric = Distances.Euclidean()
    d = trunc.(Int, Distances.pairwise(metric, s_coords, dims=1))
    coords_S1 = s_coords[Sk[1], :]
    coords_S2 = s_coords[Sk[2], :]
    coords_S3 = s_coords[Sk[3], :]
    coords_S4 = s_coords[Sk[4], :]
    #coords_S5 = s_coords[Sk[5], :]

    d1 = trunc.(Int, Distances.pairwise(metric, coords_S1, dims=1))
    d2 = trunc.(Int, Distances.pairwise(metric, coords_S2, dims=1))
    d3 = trunc.(Int, Distances.pairwise(metric, coords_S3, dims=1))
    d4 = trunc.(Int, Distances.pairwise(metric, coords_S4, dims=1))
    #d5 = trunc.(Int, Distances.pairwise(metric, coords_S5, dims=1))
    N1 = length(Sk[1])
    N2 = length(Sk[2])
    N3 = length(Sk[3])
    N4 = length(Sk[4])
    #N5 = length(Sk[5])
    p1 = Lk[1]
    p2 = Lk[2]
    p3 = Lk[3]
    p4 = Lk[4]
    #p5 = Lk[5]

    pdisp1 = pdisp_simple_grasp_new(d1, p1, N1, αₗ)
    pdisp2 = pdisp_simple_grasp_new(d2, p2, N2, αₗ)
    pdisp3 = pdisp_simple_grasp_new(d3, p3, N3, αₗ)
    pdisp4 = pdisp_simple_grasp_new(d4, p4, N4, αₗ)
    #pdisp5 = pdisp_simple_grasp_new(d5, p5, N5, αₗ)

    pdisp1_fixed = Sk[1][pdisp1]
    pdisp2_fixed = Sk[2][pdisp2]
    pdisp3_fixed = Sk[3][pdisp3]
    pdisp4_fixed = Sk[4][pdisp4]
    #pdisp5_fixed = Sk[5][pdisp5]

    N = S
    pdisp_ok = Set(vcat([pdisp1_fixed, pdisp2_fixed, pdisp3_fixed, pdisp4_fixed]...))
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

            #println("RCL para completo: ")
            #println(rcl)
            
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

function vnd_local_search(initial_sol, targets_lower_op, targets_upper_op, targets_lower, targets_upper)
    # Current best solution
    current_sol = deepcopy(initial_sol)
    
    # Define neighborhood structures and their corresponding functions
    neighborhoods = [
        (:simple_bu, (sol) -> simple_bu_improve_optimized(sol, targets_lower_op, targets_upper_op, :bf)),
        (:interchange_bu, (sol) -> interchange_bu_improve_optimized(sol, targets_lower_op, targets_upper_op, :bf)),
        (:deactivate_center, (sol) -> deactivate_center_improve(sol, targets_lower, targets_upper))
    ]
    
    k = 1  # Start with first neighborhood
    
    while k <= length(neighborhoods)
        # Get current neighborhood and its move function
        neighborhood_name, move_function = neighborhoods[k]
        
        # Apply the move
        new_sol = move_function(current_sol)
        
        # If improvement found
        if new_sol.Weight < current_sol.Weight
            current_sol = new_sol
            k = 1  # Reset to first neighborhood
        else
            k += 1  # Move to next neighborhood
        end
    end
    
    return current_sol
end


function main_grasp(;path="solucion_grasp_16_625_feas.jld2", iters=10)
    #file_name = "instances\\625_78_32\\inst_1_625_78_32.jld2"
    instance = read_instance(path)
    pattern = Regex("[t][_]\\d{1,3}")
    index = findfirst(pattern, path)
    almost_number = path[index]
    _, number = split(almost_number, "_")
    αₗ = 0.1
    αₐ = 0.1 
    iters = parse(Int, ARGS[2])
    bestSolution, totalTime = grasp(αₗ, αₐ, iters, instance)
    println(totalTime)
    println(bestSolution.Weight)
    println(isFactible(bestSolution))
    name = "solucion_grasp_newinstance_newrangeinteger_$number" * "_$αₗ" * "_$αₐ" *"_$iters" * "_$(nthreads())" * ".jld2"
    #full_out_path = "out\\solutions\\1250_155_62\\grasp\\" * name
    write_solution(bestSolution, name)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_grasp(;path=ARGS[1], iters=parse(Int, ARGS[2]))
else
    #main_grasp()
end