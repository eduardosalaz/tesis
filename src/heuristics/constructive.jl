using Types
using JuMP
using Distances
using DelimitedFiles
using Gurobi
using Dates
using TimerOutputs
using MathOptInterface
using Dates
using Random
using DataStructures
using LinearAlgebra
using Statistics

function remove!(V, item)
    deleteat!(V, findall(x -> x == item, V))
end

function constructive(instance, id, init_method, assign_method; withdir=false, dir="")
    instancia = instance
    B = instancia.B
    S = instancia.S
    P = instancia.P
    X = Matrix{Int64}[]
    Y = Vector{Int64}[]
    D = instancia.D
    Weight = 0
    println("hola 4")

    before1 = Dates.now()
    time_init1 = 0
    if init_method == "relax"
        Y, time_init1 = localize_facs(instancia, init_method)
        before1 = Dates.now()
    else
        Y, time = localize_facs(instancia, init_method)
    end

    if init_method ≠ "relax"
        Y_bool = zeros(Int, instancia.S)
        for idx in Y
            Y_bool[idx] = 1
        end
    else
        Y_bool = Y
    end
    after_init = Dates.now()
    delta_init1 = after_init - before1
    if init_method == "relax"
        delta_init1 = time_init1
    end
    println("Y done with time: ", delta_init1)
    if assign_method == "naive"
        #println("NAIVE")
        X = naive_assign_bu(instancia, Y_bool)
    elseif assign_method == "opp"
        println("OPP COST")
        X = oppCostAssignment(Y_bool, instancia)
    elseif assign_method == "strat"
        println("strat")
        X, count = strategic_allocation(Y_bool, instancia)
    elseif assign_method == "queue"
        println("queue")
        X, count = oppCostQueue2(Y_bool, instancia)
    end
    after1 = Dates.now()
    delta_assign = after1 - after_init
    println("X done with time: ", delta_assign)
    delta1 = after1 - before1
    secs1 = round(delta1, Second)
    time1 = secs1.value
    if init_method == "relax"
        time1 += time_init1
    end

    before = Dates.now()
    time_init = 0
    if init_method == "relax"
        Y, time_init = localize_facs(instancia, init_method)
        before = Dates.now()
    else
        Y, time = localize_facs(instancia, init_method)
    end
    if init_method ≠ "relax"
        Y_bool = zeros(Int, instancia.S)
        for idx in Y
            Y_bool[idx] = 1
        end
    else
        Y_bool = Y
    end
    after_init2 = Dates.now()
    delta_init2 = after_init2 - before
    if init_method == "relax"
        delta_init2 = time_init
    end
    println("Y done with time: ", delta_init2)
    if assign_method == "naive"
        #println("NAIVE")
        X = naive_assign_bu(instancia, Y_bool)
    elseif assign_method == "opp"
        println("OPP COST")
        X = oppCostAssignment(Y_bool, instancia)
    elseif assign_method == "strat"
        println("strat")
        X, count = strategic_allocation(Y_bool, instancia)
    elseif assign_method == "queue"
        println("queue")
        X, count = oppCostQueue2(Y_bool, instancia)
    end
    after = Dates.now()
    delta_assign2 = after - after_init2
    println("X done with time: ", delta_assign2)
    delta = after - before
    secs = round(delta, Second)
    time = secs.value
    if init_method == "relax"
        time += time_init
    end
    println("X done")
    println(size(X))
    indices = findall(x -> x == 1, X)
    Weight = dot(X, D)
    if time > time1
        time = time1
    end
    println(time)

    str_path = "sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * init_method * "_" * assign_method
    plot_str_path = str_path * ".png"
    solution_str_path = str_path * ".jld2"

    solution = Types.Solution(instancia, X, Y_bool, Weight, time)
    Types.plot_solution(solution, plot_str_path)
    #Types.write_solution(solution, solution_str_path)
    println(isFactible(solution, true))
    return solution
end


function pdisp_simple(d, p, N)
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
                vbest = v
            end
        end
        if vbest != 0 && !(vbest in P)

            push!(P, vbest)
        end
    end
    collection = collect(P)
    return collection
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

function localize_facs(instance, method)
    # k_type = findall(x->x==1, idx_candidate .∈ Sk) # get k type of fac
    if method == "pdisp"
        println("P-DISP")
        return pdisp_2(instance)
    elseif method == "random"
        println("RANDOM")
        return random_init(instance)
    elseif method == "relax"
        println("RELAXATION")
        return relax_init(instance)
    elseif method == "multi"
        println("MULTI TYPE PDISP")
        return multi_type_pdp_heuristic_instance(instance, instance.D)
    end
end

function pdisp_2(instance)
    before_init = Dates.now()
    S = instance.S
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    P = instance.P

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

    pdisp1 = pdisp_simple(d1, p1, N1)
    pdisp2 = pdisp_simple(d2, p2, N2)
    pdisp3 = pdisp_simple(d3, p3, N3)
    pdisp4 = pdisp_simple(d4, p4, N4)
    #pdisp5 = pdisp_simple(d5, p5, N5)

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
            # Find the node v that maximizes the distance to its closest neighbor in P
            maxdist = 0
            vbest = 0
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
                    vbest = v
                end
            end
            # If no such node exists, stop the algorithm
            if vbest == 0
                @error "PDISP FAILED"
                println("*******************************************************************************************")
                break
            end
            # Add the node vbest to the set P and update the counts
            k = node_type(vbest, Sk)
            count[k] += 1
            push!(pdisp_ok, vbest)
        end
    end
    after_init = Dates.now()
    delta_init = after_init - before_init
    delta_init_milli = round(delta_init, Millisecond)
    collection = collect(pdisp_ok)
    return collection, delta_init_milli.value
end

function multi_type_pdp_heuristic_instance(instance, d; objective="minmin", local_search=true)
    # Extract data from instance object
    n = instance.S               # Total number of points
    Sk = instance.Sk             # Sets of points by type (Sk[k] contains indices of all points of type k)
    Lk = instance.Lk             # Lower bounds for each type (minimum number to select)
    Uk = instance.Uk             # Upper bounds for each type (maximum number to select)
    p = instance.P               # Total number of points to select
    k = length(Sk)               # Number of different types
    
    # Initialize solution structures
    selected = Int[]             # Will store indices of selected points
    selected_by_type = zeros(Int, k)  # Counter for how many points of each type are selected
    
    # Phase 1: Greedy construction with type constraints - ensuring lower bounds are met
    remaining = p                # Track how many points we still need to select
    
    # First, ensure we meet lower bounds for each type
    for k_idx in 1:k             # Iterate through each type
        # Create a list of candidate points of this type that haven't been selected yet
        candidates = copy(Sk[k_idx])
        filter!(i -> !(i in selected), candidates)
        
        # Determine how many points of this type we need to select to meet lower bound
        needed = Lk[k_idx]
        
        # Keep selecting points of this type until lower bound is met
        while needed > 0 && !isempty(candidates)
            best_idx = -1        # Will store index of best candidate
            best_value = -Inf    # Will store value of best candidate
            
            # Evaluate each candidate point
            for idx in candidates
                # For the first point being selected, use distance from center as criterion
                if isempty(selected)
                    # Calculate center of all points
                    center_x = sum(instance.S_coords[:, 1]) / n
                    center_y = sum(instance.S_coords[:, 2]) / n
                    # Calculate distance from center
                    value = sqrt((instance.S_coords[idx, 1] - center_x)^2 + (instance.S_coords[idx, 2] - center_y)^2)
                else
                    # For subsequent points, evaluate based on chosen objective
                    if objective == "sumsum"
                        # Sum of distances to all already selected points
                        value = sum(d[idx, j] for j in selected)
                    else # minmin
                        # Minimum distance to any already selected point
                        value = minimum(d[idx, j] for j in selected)
                    end
                end
                
                # Keep track of the best candidate
                if value > best_value
                    best_value = value
                    best_idx = idx
                end
            end
            
            # If we found a suitable candidate, add it to our solution
            if best_idx != -1
                push!(selected, best_idx)           # Add to selected set
                selected_by_type[k_idx] += 1        # Increment counter for this type
                needed -= 1                         # One less point needed for this type
                remaining -= 1                      # One less point needed overall
                filter!(x -> x != best_idx, candidates)  # Remove from candidates list
            else
                break  # No suitable candidates found, move to next type
            end
        end
    end
    
    # Phase 2: Fill remaining slots greedily while respecting upper bounds
    while remaining > 0          # Continue until we've selected enough points
        best_idx = -1            # Will store index of best candidate
        best_value = -Inf        # Will store value of best candidate
        
        # Consider points of all types, but respect upper bounds
        for k_idx in 1:k
            # Skip if we're already at upper bound for this type
            if selected_by_type[k_idx] >= Uk[k_idx]
                continue
            end
            
            # Evaluate all unselected points of this type
            for i in Sk[k_idx]
                if !(i in selected)
                    # For the first point being selected, use distance from center
                    if isempty(selected)
                        center_x = sum(instance.S_coords[:, 1]) / n
                        center_y = sum(instance.S_coords[:, 2]) / n
                        value = sqrt((instance.S_coords[i, 1] - center_x)^2 + (instance.S_coords[i, 2] - center_y)^2)
                    else
                        # For subsequent points, evaluate based on objective
                        if objective == "sumsum"
                            value = sum(d[i, j] for j in selected)
                        else # minmin
                            value = minimum(d[i, j] for j in selected)
                        end
                    end
                    
                    # Keep track of the best candidate
                    if value > best_value
                        best_value = value
                        best_idx = i
                    end
                end
            end
        end
        
        # If we couldn't find any valid point, stop the algorithm
        if best_idx == -1
            break
        end
        
        # Add the best point to our solution
        push!(selected, best_idx)
        
        # Update the type counter for the selected point
        for k_idx in 1:k
            if best_idx in Sk[k_idx]
                selected_by_type[k_idx] += 1
                break
            end
        end
        
        remaining -= 1  # One less point needed overall
    end
    
    # Phase 3: Local search improvement - try to improve the solution through swaps
    if local_search && length(selected) >= 2
        improved = true          # Flag to track if we made an improvement
        iterations = 0           # Counter for number of iterations
        max_iterations = 100     # Maximum number of iterations allowed
        
        # Continue as long as we're still finding improvements
        while improved && iterations < max_iterations
            improved = false
            iterations += 1
            
            # Calculate current objective value
            current_obj = 0
            if objective == "sumsum"
                # Sum of all pairwise distances between selected points
                current_obj = sum(d[s1, s2] for s1 in selected for s2 in selected if s1 < s2)
            else # minmin
                # Minimum distance between any pair of selected points
                # First create all pairs of selected points
                if length(selected) >= 2
                    # Calculate minimum distance between any pair of selected points
                    min_dist = Inf
                    for i in 1:length(selected)
                        for j in i+1:length(selected)
                            dist = d[selected[i], selected[j]]
                            if dist < min_dist
                                min_dist = dist
                            end
                        end
                    end
                    current_obj = min_dist
                end
            end
            
            # Try swapping each selected point with each non-selected point of same type
            for i in selected
                # Find type of this point
                i_type = 0
                for k_idx in 1:k
                    if i in Sk[k_idx]
                        i_type = k_idx
                        break
                    end
                end
                
                # Try all non-selected points of the same type
                for j in Sk[i_type]
                    if !(j in selected)
                        # Simulate swap by creating a modified selected set
                        new_selected = copy(selected)
                        replace!(new_selected, i => j)
                        
                        # Calculate objective value with the simulated swap
                        new_obj = 0
                        if objective == "sumsum"
                            new_obj = sum(d[s1, s2] for s1 in new_selected for s2 in new_selected if s1 < s2)
                        else # minmin
                            # Calculate minimum distance between any pair of points in new selection
                            if length(new_selected) >= 2
                                min_dist = Inf
                                for i2 in 1:length(new_selected)
                                    for j2 in i2+1:length(new_selected)
                                        dist = d[new_selected[i2], new_selected[j2]]
                                        if dist < min_dist
                                            min_dist = dist
                                        end
                                    end
                                end
                                new_obj = min_dist
                            end
                        end
                        
                        # If the swap improves the objective, make it permanent
                        if new_obj > current_obj
                            replace!(selected, i => j)
                            improved = true
                            break  # Move to next selected point
                        end
                    end
                end
                
                # If we improved the solution, start over with the new solution
                if improved
                    break
                end
            end
        end
    end
    
    # Calculate final objective value for the solution
    obj_value = 0
    if length(selected) >= 2
        if objective == "sumsum"
            # Sum of all pairwise distances
            obj_value = sum(d[s1, s2] for s1 in selected for s2 in selected if s1 < s2)
        else # minmin
            # Minimum distance between any pair of selected points
            min_dist = Inf
            for i in 1:length(selected)
                for j in i+1:length(selected)
                    dist = d[selected[i], selected[j]]
                    if dist < min_dist
                        min_dist = dist
                    end
                end
            end
            obj_value = min_dist
        end
    end
    
    return selected, selected_by_type, obj_value
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

function minimums(vec::Vector, n)::Tuple{Vector{Int64},Vector{Int64}}
    type = eltype(vec)
    vals = fill(10000000000000, n)
    arr = Array{Tuple{type,Int64}}(undef, n)
    indices = Array{Int64}(undef, n)
    @inbounds for i ∈ eachindex(vec)
        if vec[i] > 0
            biggest, index = findmax(vals)
            if vec[i] < biggest
                arr[index] = vec[i], i
                vals[index] = vec[i]
            end
        end
    end
    arr = sort(arr, by=x -> x[1])
    vals = [x[1] for x in arr]
    indices = [x[2] for x in arr]
    return vals, indices
end
#132407

function maximums3(M, n)
    v = vec(M)
    l = length(v)
    ix = [1:l;]
    partialsortperm!(ix, v, (l-n+1):l, initialized=true)
    indices = CartesianIndices(M)[ix[(l-n+1):l]]
    return reverse!(indices)
end


function calculate_target(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    targets = Vector{Float32}(undef, M)
    for m in 1:M
        targets[m] = 1 * μ[m][1] * (1 + T[m]) # corregir esto
    end
    return targets
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
            targets[s, m] = floor(Int, (μ[m][s] * (1 + T[m])))
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
            targets[s, m] = ceil(Int, (μ[m][s] * (1 - T[m])))
        end
    end
    return targets
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
        if values_matrix[i, m] > targets[i, m]  # Note the i,m indexing
            constraints += 1
        end
    end
    risk_vec[i] += R[j]
    if risk_vec[i] > β[i]  # Use center-specific β
        constraints += 1
    end
    return constraints
end


function compute_assignments_and_opportunity_costs(D::Matrix{Int64}, Y::Vector{Int}, N::Int)
    _, num_clients = size(D)
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
        enqueue!(pq, j, opp_cost)
    end
    return best_assignments, pq
end



"""
1. Calculate Opportunity Costs
    Initialize oppCostMatrix with zeros
    FOR each facility i
        FOR each client j
            Calculate opportunity cost for client j with facility i
            Insert into oppCostMatrix[i][j]

2. Compute Best N Assignments
    Initialize best_assignments as an empty dictionary
    FOR each client j
        sorted_facilities = sort facilities based on oppCostMatrix[:,j]
        best_assignments[j] = first N entries of sorted_facilities

3. Prepare the priority queue with clients (based on their highest opportunity costs)
    Initialize priority_queue with all clients, priority being the highest opp cost from oppCostMatrix

4. Assignment Loop
    Initialize assigned_clients as an empty set
    WHILE priority_queue is not empty AND length of assigned_clients < number of clients
        client = dequeue priority_queue
        IF client not in assigned_clients
            FOR each facility in best_assignments[client]
                IF can_assign(facility)
                    assign(client, facility)
                    ADD client to assigned_clients
                    UPDATE capacity of facility
                    IF facility is full
                        update_best_assignments_for_all_clients(facility, best_assignments)
                    BREAK out of loop, move to next client

    Assignment Loop Pseudocode
    assigned_clients = Set() # An empty set to keep track of clients who have been assigned
    assignment_dict = Dict() # A dictionary to store which facility each client is assigned to

    WHILE !isempty(priority_queue) AND length(assigned_clients) < number_of_clients
        client, _ = dequeue(priority_queue) # Get the client with the highest opportunity cost

        IF client not in assigned_clients
            FOR facility in best_assignments[client]
                IF can_assign(facility) # Check constraints or facility capacity
                    assignment_dict[client] = facility
                    ADD client to assigned_clients
                    DECREASE capacity of facility
                    
                    IF facility is full, full is defined as when a facility reaches the minimal values required in the constraints, so it promotes the use of other centers
                        update_best_assignments_for_all_clients(facility, best_assignments)
                    END IF

                    BREAK # Move to the next client in the priority queue
                END IF
            END FOR
        END IF
    END WHILE

    RETURN assignment_dict

5. Handle Unassigned Clients
    FOR each client in clients
        IF client not in assigned_clients
            assign client to any available facility (based on some fallback heuristic)

FUNCTION update_best_assignments_for_all_clients(facility, best_assignments)
    FOR each client in best_assignments.keys()
        IF facility in best_assignments[client]
            REMOVE facility from best_assignments[client]
            IF length of best_assignments[client] < N
                Find the next best facility not in best_assignments[client] for this client
                ADD this new facility to best_assignments[client]
"""
function oppCostQueue(Y, instance::Types.Instance)
    D = copy(instance.D)
    X = zeros(Int, size(D))
    targets_lower = calculate_targets_lower(instance)
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
    N = instance.P # podemos hacerlo en proporcion a P, no necesariamente tiene que ser P
    best_assignments, pq = compute_assignments_and_opportunity_costs(D, Y, N)
    full_centers = Set()
    assigned_bus = Set()
    unassigned_bus = Set(collect(1:B))
    #println(pq)
    #println(unassigned_bus)
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        # println("dequeueing $bu")
        if bu ∉ assigned_bus
            for center in best_assignments[bu]
                full = false
                fulls_m = zeros(Int, M)
                if center ∉ full_centers
                    for m in 1:M
                        values_matrix[center, m] += V[m][bu]
                        if values_matrix[center, m] > (targets_lower[m])
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
                    X[center, bu] = 1
                    #println("asigne $bu")
                    pop!(unassigned_bus, bu)
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
    unassigned_count = length(unassigned_bus)
    println(unassigned_count)
    X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)

    return X, unassigned_count
end

# Version 1 - Using real values (probably better as you noted)
function oppCostQueue2(Y, instance::Types.Instance)
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
            target_progress[i,m] = values_matrix[i,m] / targets[i,m] # Use center-specific targets
        end
    end
    
    full_centers = Set()
    assigned_bus = Set()
    unassigned_bus = Set(collect(1:B))
    
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        if bu ∉ assigned_bus
            # Calculate which center needs this BU most
            best_center = nothing
            best_improvement = -Inf
            
            for center in best_assignments[bu]
                if center ∉ full_centers
                    # Calculate improvement in target progress
                    current_min_progress = minimum(target_progress[center,:])
        
                    temp_progress = copy(target_progress[center,:])
                    for m in 1:M
                        temp_progress[m] = (values_matrix[center,m] + V[m][bu]) / targets[center,m]
                    end
                    
                    new_min_progress = minimum(temp_progress)
                    # Weight improvement by how far we are from minimum
                    improvement = (new_min_progress - current_min_progress)

                    
                    # Check if this assignment would exceed upper bounds
                    exceeds_upper = false
                    for m in 1:M
                        if values_matrix[center,m] + V[m][bu] > targets_upper[center,m]
                            exceeds_upper = true
                            break
                        end
                    end
                    
                    # Check risk constraint
                    if !exceeds_upper && risk_vec[center] + R[bu] <= β[center] && improvement > best_improvement
                        best_improvement = improvement
                        best_center = center
                    end
                end
            end
            
            if best_center !== nothing
                # Make the assignment
                center = best_center
                for m in 1:M
                    values_matrix[center,m] += V[m][bu]
                    target_progress[center,m] = values_matrix[center,m] / targets[center, m]
                end
                risk_vec[center] += R[bu]
                
                # Check if center should be marked as full
                min_progress = minimum(target_progress[center,:])
                if min_progress >= 0.95  # Within 100% of targets, can be adjusted
                    push!(full_centers, center)
                end
                
                push!(assigned_bus, bu)
                X[center, bu] = 1
                pop!(unassigned_bus, bu)
            end
        end
    end
    unassigned_count = length(unassigned_bus)
    println(unassigned_count)
    if unassigned_count > 0
        X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    end
    
    return X, unassigned_count
end

function oppCostQueue3(Y, instance::Types.Instance)
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
    unassigned_bus = Set(collect(1:B))
   
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        if bu ∉ assigned_bus
            # Calculate which center needs this BU most
            best_center = nothing
            best_improvement = -Inf
           
            for center in best_assignments[bu]
                if center ∉ full_centers
                    # Calculate improvement in target progress
                    current_min_progress = minimum(target_progress[center,:])
       
                    temp_progress = copy(target_progress[center,:])
                    for m in 1:M
                        temp_progress[m] = (values_matrix[center,m] + V[m][bu]) / targets[center,m]
                    end
                   
                    new_min_progress = minimum(temp_progress)
                    improvement = (new_min_progress - current_min_progress)
                   
                    # Check if this assignment would exceed upper bounds
                    exceeds_upper = false
                    for m in 1:M
                        if values_matrix[center,m] + V[m][bu] > targets_upper[center,m]
                            exceeds_upper = true
                            break
                        end
                    end
                   
                    # Check risk constraint
                    if !exceeds_upper && risk_vec[center] + R[bu] <= β[center] && improvement > best_improvement
                        best_improvement = improvement
                        best_center = center
                    end
                end
            end
           
            if best_center !== nothing
                # Make the assignment
                center = best_center
                for m in 1:M
                    values_matrix[center,m] += V[m][bu]
                    target_progress[center,m] = values_matrix[center,m] / targets[center,m]
                end
                risk_vec[center] += R[bu]
               
                # NEW: Check if center should be marked as full by checking if ANY remaining unassigned
                # business unit would violate constraints if added
                can_accept_more = false
                for remaining_bu in unassigned_bus
                    if remaining_bu ∉ assigned_bus
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
                X[center, bu] = 1
                pop!(unassigned_bus, bu)
            end
        end
    end
   
    unassigned_count = length(unassigned_bus)
    println(unassigned_count)
    if unassigned_count > 0
        X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    end
   
    return X, unassigned_count
end

function oppCostQueue4(Y, instance::Instance)
    D = copy(instance.D)
    X = zeros(Int, size(D))
    M = instance.M
    S = instance.S
    B = instance.B
    V = instance.V
    R = instance.R
    β = instance.β
    
    # Initialize tracking matrices
    values_matrix = zeros(Int, S, M)
    risk_vec = zeros(Int, S)
    
    targets_lower = calculate_targets_lower(instance)
    targets_upper = calculate_targets_upper(instance)
    
    # Track centers that have met their lower bounds
    center_status = Dict{Int, Symbol}() # :open, :lower_met, :full
    for i in 1:S
        center_status[i] = :open
    end
    
    # Calculate initial assignments and costs
    N = instance.P
    best_assignments, pq = compute_assignments_and_opportunity_costs(D, Y, N)
    
    assigned_bus = Set{Int}()
    unassigned_bus = Set(collect(1:B))
    
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        
        if bu ∉ assigned_bus
            for center in best_assignments[bu]
                if get(center_status, center, :full) == :full
                    continue
                end
                
                valid_assignment = true
                
                # Check risk constraint
                if risk_vec[center] + R[bu] > β[center]
                    valid_assignment = false
                    continue
                end
                
                # Check activity constraints
                new_values = [values_matrix[center, m] + V[m][bu] for m in 1:M]
                for m in 1:M
                    if new_values[m] > targets_upper[center, m]
                        valid_assignment = false
                        break
                    end
                end
                
                if valid_assignment
                    # Make assignment
                    X[center, bu] = 1
                    push!(assigned_bus, bu)
                    pop!(unassigned_bus, bu)
                    
                    # Update matrices
                    for m in 1:M
                        values_matrix[center, m] += V[m][bu]
                    end
                    risk_vec[center] += R[bu]
                    
                    # Check if center has met lower bounds
                    lower_bounds_met = true
                    for m in 1:M
                        if values_matrix[center, m] < targets_lower[center, m]
                            lower_bounds_met = false
                            break
                        end
                    end
                    
                    if lower_bounds_met && center_status[center] == :open
                        center_status[center] = :lower_met
                    end
                    
                    # Check if center is near full (close to upper bounds)
                    is_full = true
                    for m in 1:M
                        remaining_capacity = targets_upper[center, m] - values_matrix[center, m]
                        if remaining_capacity > maximum([V[m][j] for j in unassigned_bus])
                            is_full = false
                            break
                        end
                    end
                    
                    if is_full
                        center_status[center] = :full
                    end
                    
                    break
                end
            end
        end
    end
    
    unassigned_count = length(unassigned_bus)
    X = handle_unassigned_clients3(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    return X, unassigned_count
end


# Function to calculate how much an assignment violates constraints
function calculate_violation(facility, client, targets_upper, V, M, R, β, values_matrix, risk_vec)
    violation = 0
    for m in 1:M
        excess = values_matrix[facility, m] + V[m][client] - targets_upper[facility, m]
        if excess > 0
            violation += excess
        end
    end
    excess_risk = risk_vec[facility] + R[client] - β
    if excess_risk > 0
        violation += excess_risk
    end
    return violation
end


function handle_unassigned_clients2(X, instance, best_assignments, unassigned_clients, values_matrix, risk_vec)
    targets_upper = calculate_targets_upper(instance)
    S = instance.S
    M = instance.M
    V = instance.V
    β = instance.β
    R = instance.R
    counter = 0
    # Only loop over unassigned clients
    for client in unassigned_clients
        assigned = false
        # Sort facilities for this client based on their remaining risk capacity
        sorted_best_assignments = sort(best_assignments[client], by=f -> instance.β[f] - risk_vec[f], rev=true)
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
            min_violation = 10000000
            least_violating_facility = -1
            for facility in sorted_best_assignments
                violation = calculate_violation(facility, client, targets_upper, V, M, R, β[facility], values_matrix, risk_vec)
                #println(violation, facility)
                if violation < min_violation
                    min_violation = violation
                    least_violating_facility = facility
                end
            end
            # Assign the client to the least violating facility
            X[least_violating_facility, client] = 1
            for m in 1:M
                values_matrix[least_violating_facility, m] += V[m][client]
            end
            risk_vec[least_violating_facility] += R[client]
            counter += 1
            println("Client $client was assigned to facility $least_violating_facility with a violation of $min_violation.")
        end
    end
    println("Number of unassigned clients: $(length(unassigned_clients)), served $counter")
    return X
end

function handle_unassigned_clients3(X, instance, best_assignments, unassigned_clients, values_matrix, risk_vec)

    function make_assignment!(X, facility, client, values_matrix, risk_vec, V, M, R)
        # Update assignment matrix
        X[facility, client] = 1
        
        # Update activity values for all measures
        for m in 1:M
            values_matrix[facility, m] += V[m][client]
        end
        
        # Update risk vector
        risk_vec[facility] += R[client]
    end

    function is_valid_assignment(facility, client, targets_upper, targets_lower, V, M, R, facility_risk_cap, values_matrix, risk_vec)
        # Check risk capacity constraint
        if risk_vec[facility] + R[client] > facility_risk_cap
            return false
        end
        
        # Check all activity measures constraints
        for m in 1:M
            new_value = values_matrix[facility, m] + V[m][client]
            
            # Check against upper bound
            if new_value > targets_upper[m]
                return false
            end
        end
        
        return true
    end
    targets_upper = calculate_targets_upper(instance)
    targets_lower = calculate_targets_lower(instance)
    S = instance.S
    M = instance.M
    V = instance.V
    β = instance.β
    R = instance.R
    
    # Track violation counts per facility
    facility_violation_counts = zeros(Int, S)
    
    for client in unassigned_clients
        assigned = false
        # First try: Find valid assignment
        sorted_best_assignments = sort(best_assignments[client], 
                                     by=f -> (instance.β[f] - risk_vec[f], -facility_violation_counts[f]), 
                                     rev=true)
        
        # Try valid assignments first
        for facility in sorted_best_assignments
            if is_valid_assignment(facility, client, targets_upper, targets_lower, V, M, R, β[facility], values_matrix, risk_vec)
                make_assignment!(X, facility, client, values_matrix, risk_vec, V, M, R)
                assigned = true
                break
            end
        end
        
        # If no valid assignment found, find least violating facility
        if !assigned
            min_violation = typemax(Int)
            least_violating_facility = -1
            
            for facility in sorted_best_assignments
                # Consider both current facility violations and new violation
                total_violation = calculate_violation(facility, client, targets_upper, V, M, R, β[facility], values_matrix, risk_vec)
                weighted_violation = total_violation * (1 + facility_violation_counts[facility] * 0.5)
                
                if weighted_violation < min_violation
                    min_violation = weighted_violation
                    least_violating_facility = facility
                end
            end
            
            # Ensure we assign to somewhere, even with violations
            if least_violating_facility != -1
                make_assignment!(X, least_violating_facility, client, values_matrix, risk_vec, V, M, R)
                facility_violation_counts[least_violating_facility] += 1
                println("Client $client was assigned to facility $least_violating_facility with a violation of $min_violation")
            else
                # Absolute fallback: assign to closest facility regardless of violations
                closest_facility = argmin([instance.D[f, client] for f in 1:S])
                make_assignment!(X, closest_facility, client, values_matrix, risk_vec, V, M, R)
                println("FALLBACK: Client $client forced assignment to closest facility $closest_facility")
            end
        end
    end
    
    # Verify all clients are assigned
    for j in 1:instance.B
        if sum(X[:, j]) == 0
            println("WARNING: Client $j still unassigned after fallback!")
            # Force assign to closest facility as absolute last resort
            closest_facility = argmin(instance.D[:, j])
            make_assignment!(X, closest_facility, j, values_matrix, risk_vec, V, M, R)
            println("EMERGENCY: Forced assignment of client $j to facility $closest_facility")
        end
    end
    
    return X
end

#=
 \STATE \textbf{Input:} Distance matrix $d$, Decision Variable $Y$, $\lvert B \rvert$, $\lvert S \rvert$
        \STATE \textbf{Output:} Opportunity Cost Matrix
         \Procedure{OppCostMatrix}
            \FOR{$j$ \in $1$ \to $ \lvert B \rvert$}
                \STATE $minimal \gets$ \textbf{argmin}($d_[:, j]$)
                \FOR{$i$ \in $1$ \to $\lvert S \rvert$}
                    \STATE $oppMatrix[j, i] \gets D[j, i] - minimal$
                \ENDFOR
            \ENDFOR
            \STATE \textbf{return} $oppMatrix$
        \EndProcedure
=#

struct ActivityMetric
    avg::Float64
    min::Float64
    max::Float64
    range::Float64
end

mutable struct FacilityState
    id::Int
    current_values::Vector{Float64}
    min_thresholds::Vector{Float64}
    max_thresholds::Vector{Float64}
    normalized_gaps::Vector{Float64}
    current_risk::Float64
    risk_capacity::Float64
end

function get_activity_metrics(instance::Instance)
    activity_metrics = Vector{ActivityMetric}(undef, instance.M)
    
    for m in 1:instance.M
        values = Float64.(instance.V[m])  # Convert to Float64 for precision
        activity_metrics[m] = ActivityMetric(
            mean(values),
            minimum(values),
            maximum(values),
            maximum(values) - minimum(values)
        )
    end
    
    return activity_metrics
end

function initialize_facility_states(Y::Vector{Int}, instance::Instance)
    facility_states = Dict{Int, FacilityState}()
    
    for f in 1:instance.S
        if Y[f] == 1
            min_thresholds = [round(Int, instance.μ[m][f] * (1 - instance.T[m])) for m in 1:instance.M]
            max_thresholds = [round(Int, instance.μ[m][f] * (1 + instance.T[m])) for m in 1:instance.M]
            
            # Initial gaps are just the minimum thresholds since we start at zero
            normalized_gaps = min_thresholds  # Will be normalized in the calling function
            
            facility_states[f] = FacilityState(
                f,                      # id
                zeros(instance.M),      # current_values
                min_thresholds,         # min_thresholds
                max_thresholds,         # max_thresholds
                normalized_gaps,        # normalized_gaps
                0.0,                    # current_risk
                instance.β[f]           # risk_capacity
            )
        end
    end
    
    return facility_states
end

function get_most_critical_facility(facility_states)
    critical_facilities = []
    
    for (_, state) in facility_states
        max_gap = maximum(state.normalized_gaps)
        if max_gap > 0  # Has at least one activity below minimum
            push!(critical_facilities, (facility=state, priority=max_gap))
        end
    end
    
    isempty(critical_facilities) && return nothing
    
    # Return facility with highest priority (largest gap)
    return sort(critical_facilities, by = x -> x.priority, rev=true)[1].facility
end

function select_best_bu(facility_state, unassigned_bus, instance, activity_metrics)
    critical_activity = argmax(facility_state.normalized_gaps)
    
    best_bu = nothing
    best_contribution = -Inf
    
    for bu in unassigned_bus
        # Check risk constraint first
        if facility_state.current_risk + instance.R[bu] > facility_state.risk_capacity
            continue
        end
        
        # Evaluate contribution to critical activity
        contribution = instance.V[critical_activity][bu] / activity_metrics[critical_activity].avg
        
        # Check if this would exceed maximum thresholds
        would_exceed = false
        for m in 1:instance.M
            new_value = facility_state.current_values[m] + instance.V[m][bu]
            if new_value > facility_state.max_thresholds[m]
                would_exceed = true
                break
            end
        end
        
        if !would_exceed && contribution > best_contribution
            best_contribution = contribution
            best_bu = bu
        end
    end
    
    return best_bu
end

function make_assignment!(X, facility_id, bu, facility_states, instance)
    # Update assignment matrix
    X[facility_id, bu] = 1
    
    # Update facility state
    state = facility_states[facility_id]
    
    # Update values and risk
    for m in 1:instance.M
        state.current_values[m] += instance.V[m][bu]
    end
    state.current_risk += instance.R[bu]
    
    # Update normalized gaps
    for m in 1:instance.M
        gap = state.min_thresholds[m] - state.current_values[m]
        state.normalized_gaps[m] = max(0, gap)
    end
end

function handle_balanced_phase(facility_states, unassigned_bus, instance, activity_metrics)
    # When no critical gaps, allocate based on distance to max thresholds
    best_score = -Inf
    best_pair = nothing
    
    for (_, state) in facility_states
        # Skip if risk capacity exhausted
        remaining_risk = state.risk_capacity - state.current_risk
        remaining_risk <= 0 && continue
        
        for bu in unassigned_bus
            # Basic feasibility check
            instance.R[bu] > remaining_risk && continue
            
            # Check maximum thresholds
            would_exceed = false
            for m in 1:instance.M
                if state.current_values[m] + instance.V[m][bu] > state.max_thresholds[m]
                    would_exceed = true
                    break
                end
            end
            would_exceed && continue
            
            # Simple scoring: prefer assignments furthest from max thresholds
            score = minimum(
                (state.max_thresholds[m] - (state.current_values[m] + instance.V[m][bu])) / 
                activity_metrics[m].avg for m in 1:instance.M
            )
            
            if score > best_score
                best_score = score
                best_pair = (facility=state.id, bu=bu)
            end
        end
    end
    
    return best_pair
end

function strategic_allocation(Y, instance::Instance)
    println("Starting allocation for instance with $(instance.B) BUs and $(sum(Y)) open facilities")
    
    X = zeros(Int, instance.S, instance.B)
    activity_metrics = get_activity_metrics(instance)
    println("Activity metrics:")
    for (m, metric) in enumerate(activity_metrics)
        println("Activity $m: avg=$(round(metric.avg, digits=2)), range=[$(metric.min), $(metric.max)]")
    end
    
    facility_states = initialize_facility_states(Y, instance)
    unassigned_bus = Set(1:instance.B)
    
    iteration = 1
    while !isempty(unassigned_bus)
        println("\nIteration $iteration: $(length(unassigned_bus)) BUs remaining")
        
        critical_facility = get_most_critical_facility(facility_states)
        if critical_facility !== nothing
            println("Critical facility $(critical_facility.id) found")
            println("Normalized gaps: ", round.(critical_facility.normalized_gaps, digits=3))
            
            best_bu = select_best_bu(critical_facility, unassigned_bus, instance, activity_metrics)
            if best_bu !== nothing
                println("Selected BU $best_bu for critical facility")
                println("Before assignment - Values: ", round.(critical_facility.current_values, digits=2))
                make_assignment!(X, critical_facility.id, best_bu, facility_states, instance)
                println("After assignment - Values: ", round.(critical_facility.current_values, digits=2))
                delete!(unassigned_bus, best_bu)
            else
                println("No feasible BU found for critical facility")
                assignment = handle_balanced_phase(facility_states, unassigned_bus, instance, activity_metrics)
                if isnothing(assignment)
                    println("No feasible assignments possible")
                    break
                end
                println("Balanced assignment: facility $(assignment.facility), BU $(assignment.bu)")
                make_assignment!(X, assignment.facility, assignment.bu, facility_states, instance)
                delete!(unassigned_bus, assignment.bu)
            end
        else
            println("No critical facilities - entering balanced phase")
            assignment = handle_balanced_phase(facility_states, unassigned_bus, instance, activity_metrics)
            isnothing(assignment) && break
            println("Balanced assignment: facility $(assignment.facility), BU $(assignment.bu)")
            make_assignment!(X, assignment.facility, assignment.bu, facility_states, instance)
            delete!(unassigned_bus, assignment.bu)
        end
        
        if iteration % 10 == 0
            println("\nFacility states after iteration $iteration:")
            for (f, state) in facility_states
                println("Facility $f:")
                println("  Values: ", round.(state.current_values, digits=2))
                println("  Gaps: ", round.(state.normalized_gaps, digits=3))
                println("  Risk: $(round(state.current_risk, digits=2))/$(state.risk_capacity)")
            end
        end
        
        iteration += 1
    end
    println("\nAllocation completed. Unassigned BUs: $(length(unassigned_bus))")
    handle_unassigned_bus!(X, unassigned_bus, facility_states, instance, activity_metrics)
    println(size(X))
    return X, 0
end

function handle_unassigned_bus!(X, unassigned_bus, facility_states, instance, activity_metrics)
    println("\nHandling $(length(unassigned_bus)) unassigned BUs")
    
    function calculate_violation(facility_state, bu)
        total_violation = 0.0
        
        # Check activity violations
        for m in 1:instance.M
            new_value = facility_state.current_values[m] + instance.V[m][bu]
            
            if new_value < facility_state.min_thresholds[m]
                total_violation += (facility_state.min_thresholds[m] - new_value) / activity_metrics[m].avg
            elseif new_value > facility_state.max_thresholds[m]
                total_violation += (new_value - facility_state.max_thresholds[m]) / activity_metrics[m].avg
            end
        end
        
        # Check risk violation
        if facility_state.current_risk + instance.R[bu] > facility_state.risk_capacity
            total_violation += (facility_state.current_risk + instance.R[bu] - facility_state.risk_capacity) / 
                             facility_state.risk_capacity
        end
        
        return total_violation
    end
    
    for bu in unassigned_bus
        min_violation = Inf
        best_facility = nothing
        
        # Try each facility
        for (_, state) in facility_states
            violation = calculate_violation(state, bu)
            if violation < min_violation
                min_violation = violation
                best_facility = state
            end
        end
        
        println("BU $bu assigned to facility $(best_facility.id) with violation $min_violation")
        make_assignment!(X, best_facility.id, bu, facility_states, instance)
    end
end


function oppCostAssignment(Y, instance::Types.Instance)
    D = copy(instance.D)
    P = instance.P
    S = instance.S
    B = instance.B
    M = instance.M
    V = instance.V
    M = instance.M
    R = instance.R
    β = instance.β[1]
    X = zeros(Int64, S, B)
    count = 0
    not_assigned_y = findall(y -> y == 0, Y)
    #println(not_assigned_y)
    for j in not_assigned_y
        D[j, :] .= -1
    end
    oppMatrix = copy(D)
    for i in 1:B
        minimal = 0
        minimals, _ = minimums(D[:, i], 1) # cambiar por findmin en 1 solo
        minimal = minimals[1]
        oppMatrix[:, i] .= D[:, i] .- minimal
    end
    count = 0
    todos = false
    targets = calculate_target(instance)
    values_matrix = Matrix{Float32}(undef, S, M)
    risk_vec = Vector{Float32}(undef, S)
    n = round(Int, (P / 2)) # falta tweakear
    while !todos
        indices::Vector{CartesianIndex{2}} = maximums3(oppMatrix, n)
        values_matrix, risk_vec = start_constraints_optimized_v5(S, B, M, V, R, X, values_matrix, risk_vec)
        constraints = Int64[]
        # aqui calculamos inicialmente las constraints
        # vjm <= μim para cada i en 1:S, para cada m en 1:3, son S*3 constraints
        # Rj > βi para cada i en 1:S, S constraints
        for indice::CartesianIndex in indices
            #X_copy = Matrix{Int64}(undef, S, B)
            #unsafe_copyto!(X_copy, 1, X, 1, S * B)
            col = indice[2]
            row = findfirst(x -> x == 0, oppMatrix[:, col])
            #X[row, col] = 1
            constraints_v = update_constraints(M, V, R, row, col, values_matrix, risk_vec, targets, β)
            push!(constraints, constraints_v)
        end
        picked = CartesianIndex(1, 1)::CartesianIndex{2}
        if all(x -> x == 0, constraints) # si no se violan constraints, agarra el 0 de la columna del costo maximo
            indice = indices[1]::CartesianIndex{2}
            col = indice[2]::Int64
            row = findfirst(x -> x == 0, oppMatrix[:, col])
            picked = CartesianIndex(row, col)
        else
            if all(x -> x ≠ 0, constraints) # si todas violan constraints
                original_cons, idx = findmin(constraints) # agarra el que viole menos constraints
                indice = indices[idx]
                col = indice[2]
                original_row = findfirst(x -> x == 0, oppMatrix[:, col])::Int64
                # queremos buscar en esa columna el siguiente valor más cercano a 0 en oppMatrix
                # al hacerlo, nos acercamos al minimo valor de la matriz de distancias
                busqueda = round(Int, (P / 2)) # a lo mejor cambiarlo despues
                _, idxs_inner::Array{Int64} = minimums(oppMatrix[:, col], busqueda)
                # el primer minimo es el 0 de nuestra localizacion optima
                # lo desechamos para darle variedad a la busqueda
                idxs_inner2::Array{Int64} = idxs_inner[2:end]
                for row in idxs_inner2
                    #X_copy = Matrix{Int64}(undef, S, B)
                    #unsafe_copyto!(X_copy, 1, X, 1, S * B)
                    try
                        #X_copy[row, col] = 1
                        constraints_v2 = update_constraints(M, V, R, row, col, values_matrix, risk_vec, targets, β)
                        if constraints_v2 < original_cons
                            original_cons = constraints_v2
                            original_row = row
                            break
                        end
                    catch
                        # no se porque hay un try catch aqui
                    end
                end
                picked = CartesianIndex(original_row, col)
            else # si hay una que no viola constraints
                for idx in eachindex(constraints)
                    if constraints[idx] == 0 # agarrala
                        indice = indices[idx]
                        col = indice[2]
                        row = findfirst(x -> x == 0, oppMatrix[:, col])
                        picked = CartesianIndex(row, col)
                        break
                    end
                end
            end
        end
        X[picked] = 1
        column = picked[2]::Int64
        oppMatrix[:, column] .= -1 # "apagamos" esa columna
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

function pdisp(instance)
    P = instance.P
    s_coords = instance.S_coords
    S = instance.S
    metric = Distances.Euclidean()
    s_distances = trunc.(Int, Distances.pairwise(metric, s_coords, dims=1))
    idx_1, idx_2 = Tuple(argmax(s_distances)) # two furthest facs
    S_sol = [idx_1, idx_2] # S is our solution of our p-disp problem
    T = collect(1:S) # nodes not in our solution
    # se puede usar la suma de xⱼ a S_sol
    # o usar la facility mas cercana de xⱼ
    for s in S_sol
        remove!(T, s)
    end

    while length(S_sol) < P
        distances = Tuple{Int64,Int64}[]
        idx_max_dist = 0
        for i in T
            distance = 0
            for j in S
                distance += s_distances[i, j]
            end
            tupla = (distance, i)
            push!(distances, tupla)
        end
        distances_max = [distance[1] for distance in distances]::Array{Int64}
        indexes = [distance[2] for distance in distances]::Array{Int64}
        idx_max = argmax(distances_max)::Int64
        idx_max_dist = indexes[idx_max]
        remove!(T, idx_max_dist)
        push!(S_sol, idx_max_dist)
    end
    return S_sol
end

function random_init(instance)
    # tengo que agarrar los parametros para asignar una factible
    P = instance.P
    S = instance.S
    Y = Random.shuffle(collect(1:S))
    Y = Y[1:P]
    return Y
end

function relax_init(instance)
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
    m = instance.M
    k = instance.K

    model = JuMP.Model(Gurobi.Optimizer) # THIS IS WHERE THE FUN BEGINS

    if solver_name(model) == "Gurobi"
        MOI.set(model, MOI.RawOptimizerAttribute("SolutionLimit"), 1) # retorna la primera factible encontrada, solo funciona para gurobi    
    end

    JuMP.@variable(model, x[1:S, 1:B], lower_bound = 0, upper_bound = 1)
    # num suc and num bu, Xᵢⱼ
    JuMP.@variable(model, y[1:S], Bin)
    # Yᵢ

    JuMP.@objective(model, Min, sum(D .* x))
    # Xᵢⱼ * Dᵢⱼ

    JuMP.@constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)

    # ∑ᵢ∈S Xᵢⱼ = 1, ∀ j ∈ B

    JuMP.@constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y[i])

    # Xᵢⱼ ≤ Yᵢ , ∀ i ∈ S, j ∈ B

    JuMP.@constraint(model, cardinality, sum(y) == P)

    # ∑ i ∈ S Yᵢ = p

    JuMP.@constraint(model, risk[j in 1:B, i in 1:S], x[i, j] * R[j] <= β[i])

    # ∑ j ∈ B Xᵢⱼ Rⱼ ≤ βᵢ, ∀ i ∈ S

    JuMP.@constraint(
        model,
        tol_l[i in 1:S, M in 1:m],
        y[i] * μ[M][i] * (1 - T[M]) <= sum(x[i, j] * V[m][j] for j in 1:B),
    )

    JuMP.@constraint(
        model,
        tol_u[i in 1:S, M in 1:m],
        sum(x[i, j] * V[m][j] for j in 1:B) <= y[i] * μ[M][i] * (1 + T[M]),
    )

    # Yᵢμₘⁱ(1-tᵐ) ≤ ∑i∈S Xᵢⱼ vⱼᵐ ≤ Yᵢμₘʲ(1+tᵐ) ∀ j ∈ B, m = 1 … 3

    JuMP.@constraint(
        model,
        low_k[K in 1:k],
        Lk[K] <= sum(y[i] for i in Sk[K]),
    )

    JuMP.@constraint(
        model,
        upp_k[K in 1:k],
        sum(y[i] for i in Sk[K]) <= Uk[K],
    )

    JuMP.set_silent(model)
    JuMP.optimize!(model)
    tiempo = MOI.get(model, MOI.SolveTimeSec())
    time_int = round(Int, tiempo)
    Y = round.(Int, JuMP.value.(model[:y]))
    return Y, time_int
end

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
    usables_i = Set(findall(==(1), Y))
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
            println("Violando número de centros asignados ", sum(Y), " $P")
        end
        number_constraints_violated += 1
    end

    for j in 1:B
        if sum(X[i, j] for i in 1:S) ≠ 1
            if verbose
                println("Violando servicio a la BU: ", j)
                number_constraints_violated += 1
            end

        end
    end

    for j in 1:B
        assigned = findfirst(==(1), X[:, j])
        if assigned ∉ usables_i
            println("Violando $j asignada a una y no abierta $assigned")
            number_constraints_violated += 1
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
            if !(ceil(Int, (Y[i] * μ[m][i] * (1 - T[m]))) <= sum(X[i, j] * V[m][j] for j in 1:B))
                if verbose
                    println("violando V inferior en i: $i y m: $m")
                    println("μ: ", round(Int, (Y[i] * μ[m][i] * (1 - T[m]))))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= floor(Int, Y[i] * μ[m][i] * (1 + T[m])))
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ", round(Int, (Y[i] * μ[m][i] * (1 + T[m]))))
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

function main_constructive(init_method, assign_method; path="inst", read_file=true, instance_obj=nothing, id=0)
    println("hola 2")
    instance = 1 # para traerlo al scope
    if read_file
        pattern = Regex("[t][_]\\d{1,3}")
        index = findfirst(pattern, path)
        almost_number = path[index]
        _, number = split(almost_number, "_")
        instance = Types.read_instance(path)
    else
        instance = instance_obj
        number = id
    end
    B = instance.B
    P = instance.P
    S = instance.S
    println("hola 3")

    solution = constructive(instance, number, init_method, assign_method)
    println(isFactible(solution))
    sol_path = "sol_$number" * "_$B" * "_$S" * "_$P" * "_$init_method" * "_$assign_method" * "_int_010_010_010_newranges_queue2.jld2"
    write_solution(solution, sol_path)
    println("\a")
    return solution
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("hola")
    main_constructive(ARGS[2], ARGS[3]; path=ARGS[1])
end
