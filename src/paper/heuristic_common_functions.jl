# heuristic_common_functions.jl
# Shared functions between Y initialization and local search ablation experiments

# Common functions for both experiments


function calculate_dist_min(D)
    B = size(D, 2)
    dist_min = zeros(B)
    for j in 1:B
        dist_min[j] = minimum(D[:, j])
    end
    return dist_min
end

function calculate_p_max(D, P)
    B = size(D, 2)
    p_max = zeros(B)
    for j in 1:B
        sorted_dists = sort(D[:, j])
        p_max[j] = sorted_dists[min(P, end)]
    end
    return p_max
end

function solve_multi_type_p_median_strong(instance, time_limit, log_file)
    """
    Strong formulation of the multi-type p-median problem
    Returns only the Y vector (binary array indicating selected centers)
    """
    S = instance.S  # number of centers
    P = instance.P  # number to select
    B = instance.B  # number of branches
    D = instance.D  # distance matrix
    
    # Create model
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", time_limit)
    set_optimizer_attribute(model, "MIPGap", 1e-4)
    set_optimizer_attribute(model, "LogFile", log_file)
    #set_optimizer_attribute(model, "OutputFlag", 0)  # Suppress output
    
    # Variables
    @variable(model, y[1:S], Bin)  # 1 if center i is selected
    @variable(model, x[1:S, 1:B], Bin)  # 1 if center i serves branch j
    
    # Objective: minimize total distance
    @objective(model, Min, sum(D[i,j] * x[i,j] for i in 1:S, j in 1:B))
    
    # Each branch must be served by exactly one center
    @constraint(model, serve[j in 1:B], sum(x[i,j] for i in 1:S) == 1)
    
    # Can only assign to open centers
    @constraint(model, link[i in 1:S, j in 1:B], x[i,j] <= y[i])
    
    # Select exactly P centers
    @constraint(model, cardinality, sum(y) == P)
    
    # Type constraints
    k = instance.K
    @constraint(model, low_k[K in 1:k], 
        instance.Lk[K] <= sum(y[i] for i in instance.Sk[K]))
    @constraint(model, upp_k[K in 1:k], 
        sum(y[i] for i in instance.Sk[K]) <= instance.Uk[K])
    
    # Solve
    optimize!(model)
    
    status = termination_status(model)
    if status in [MOI.OPTIMAL, MOI.TIME_LIMIT]
        Y_solution = round.(Int, value.(y))
        return Y_solution
    else
        return nothing
    end
end

function solve_multi_type_p_median_benders(instance, time_limit, log_file)
    """
    Benders decomposition for the multi-type p-median problem
    Returns only the Y vector (binary array indicating selected centers)
    """
    S = instance.S
    P = instance.P
    B = instance.B
    D = instance.D
    
    # Preprocess distance bounds
    dist_min = calculate_dist_min(D)
    p_max = calculate_p_max(D, P)
    
    # Master problem
    master = Model(Gurobi.Optimizer)
    set_optimizer_attribute(master, "TimeLimit", time_limit)
    #set_optimizer_attribute(master, "MIPGap", 1e-4)
    set_optimizer_attribute(master, "LogFile", log_file)
    #set_optimizer_attribute(master, "OutputFlag", 0)
    
    @variable(master, y[1:S], Bin)
    @variable(master, z[j in 1:B], lower_bound=dist_min[j], upper_bound=p_max[j])
    
    # Objective
    @objective(master, Min, sum(z))
    
    # Select exactly P centers
    @constraint(master, sum(y) == P)
    
    # Type constraints
    k = instance.K
    @constraint(master, low_k[K in 1:k],
        instance.Lk[K] <= sum(y[i] for i in instance.Sk[K]))
    @constraint(master, upp_k[K in 1:k],
        sum(y[i] for i in instance.Sk[K]) <= instance.Uk[K])
    
    # Set up for lazy constraints
    set_optimizer_attribute(master, "LazyConstraints", 1)
    
    # Benders callback
    function benders_callback(cb_data, cb_where::Cint)
        if cb_where == GRB_CB_MIPSOL
            Gurobi.load_callback_variable_primal(cb_data, cb_where)
            ȳ = callback_value.(cb_data, master[:y])
            
            sub = Model(Gurobi.Optimizer)
            set_optimizer_attribute(sub, "OutputFlag", 0)
            set_optimizer_attribute(sub, "InfUnbdInfo", 0)
            
            # Add bounds to dual variables
            @variable(sub, u >= 0)
            @variable(sub, v[1:S] >= 0)
            
            # Initialize constraints
            dual_constraints = Dict{Int, ConstraintRef}()
            for i in 1:S
                dual_constraints[i] = @constraint(sub, u - v[i] <= 0)
            end
            
            cuts = []
            for j in 1:B
                # Update RHS of constraints
                for i in 1:S
                    set_normalized_rhs(dual_constraints[i], D[i, j])
                end
                
                # Set objective - High pareto density cut
                @objective(sub, Max, (1-1/S)*(u - sum(ȳ[i] * v[i] for i in 1:S)))
                
                optimize!(sub)
                
                if termination_status(sub) == MOI.OPTIMAL
                    ū = value(u)
                    v̄ = value.(v)
                    cut_expr = @build_constraint(z[j] >= ū - sum(y[i] * v̄[i] for i in 1:S))
                    push!(cuts, cut_expr)
                end
            end
            
            for cut in cuts
                MOI.submit(master, MOI.LazyConstraint(cb_data), cut)
            end
        end
    end
    
    MOI.set(master, Gurobi.CallbackFunction(), benders_callback)
    
    # Solve
    optimize!(master)
    
    status = termination_status(master)
    if status in [MOI.OPTIMAL, MOI.TIME_LIMIT]
        Y_solution = round.(Int, value.(master[:y]))
        return Y_solution
    else
        return nothing
    end
end

# Wrapper function to choose between strong and Benders formulation
function solve_multi_type_p_median_for_heuristic(instance, time_limit, log_file; use_benders=false)
    """
    Main entry point for the multi-type p-median problem
    
    Parameters:
    - instance: Instance object with all problem data
    - time_limit: Time limit in seconds
    - log_file: Path to Gurobi log file
    - use_benders: If true, use Benders decomposition; otherwise use strong formulation
    
    Returns:
    - Y: Binary vector of selected centers, or nothing if failed
    """
    if use_benders
        return solve_multi_type_p_median_benders(instance, time_limit, log_file)
    else
        return solve_multi_type_p_median_strong(instance, time_limit, log_file)
    end
end


function pdisp_heuristic(instance)
    Y_integer, time = pdisp_2(instance)
    Y_bool = zeros(Int, instance.S)
    Y_bool[Y_integer] .= 1
    return Y_bool
end

function solve_phase2_model(instance::Instance, y_fixed)
    B = instance.B
    S = instance.S
    D = instance.D
    V = instance.V
    μ = instance.μ
    R = instance.R
    β = instance.β
    m = 3 # activities
    ε = 1e-2
    
    model = Model()
    
    @variable(model, 0 <= x[1:S, 1:B] <= 1)
    
    @objective(model, Min, sum(D .* x))
    
    @constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)
    @constraint(model, risk[i in 1:S], sum(x[i, j] * R[j] for j in 1:B) <= β[i] * y_fixed[i])
    
    @constraint(
        model,
        tol_lb[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) >= (1-ε) * y_fixed[i] * μ[M][i]
    )
    @constraint(
        model,
        tol_ub[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) <= (1+ε) * y_fixed[i] * μ[M][i]
    )
    return model
end

function build_model_x_relaxed(instance::Instance)
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

    m = 3 # activities
    k = 4 # of branches

    model = Model()

    @variable(model, x[1:S, 1:B], lower_bound = 0, upper_bound = 1)
    @variable(model, y[1:S], Bin)

    @objective(model, Min, sum(D .* x))

    @constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)
    @constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y[i])
    @constraint(model, cardinality, sum(y) == P)
    @constraint(model, risk[i in 1:S], sum(x[i, j] * R[j] for j in 1:B) <= β[i])

    @constraint(
        model,
        tol_l[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) >= y[i] * ceil(Int, μ[M][i] * (1 - T[M]))
    )
    @constraint(
        model,
        tol_u[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) <= y[i] * floor(Int, μ[M][i] * (1 + T[M]))
    )

    @constraint(
        model,
        low_k[K in 1:k],
        Lk[K] <= sum(y[i] for i in Sk[K]),
    )

    @constraint(
        model,
        upp_k[K in 1:k],
        sum(y[i] for i in Sk[K]) <= Uk[K],
    )

    return model
end

function analyze_splits(x_continuous, S, B)
    split_count = 0
    split_branches = Int[]
    tolerance = 1e-6
    
    for j in 1:B
        non_zero_assignments = 0
        
        for i in 1:S
            if x_continuous[i, j] > tolerance
                non_zero_assignments += 1
            end
        end
        
        if non_zero_assignments > 1
            split_count += 1
            push!(split_branches, j)
        end
    end
    
    split_proportion = split_count / B
    
    return split_count, split_proportion, split_branches
end

function split_resolution_heuristic(x_continuous, S, B)
    x_binary = zeros(Int, S, B)
    
    for j in 1:B
        max_val = 0.0
        max_idx = 0
        
        for i in 1:S
            if x_continuous[i, j] > max_val
                max_val = x_continuous[i, j]
                max_idx = i
            end
        end
        
        x_binary[max_idx, j] = 1
    end
    
    return x_binary
end

function multi_pdp_min_instance(instance)
    n = instance.S
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    p = instance.P
    
    metric = Distances.Euclidean()
    d = trunc.(Int, Distances.pairwise(metric, instance.S_coords, dims=1))
    
    Dmax = maximum(d)
    k = length(Sk)
    
    model = Model(Gurobi.Optimizer)
    
    @variable(model, y[1:n], Bin)
    @variable(model, u >= 0)
    @variable(model, w[1:k], Int)
    
    @objective(model, Max, u)
    
    for i in 1:n
        for j in i+1:n
            @constraint(model, u <= d[i,j] + Dmax * (2 - (y[i] + y[j])))
        end
    end
    
    @constraint(model, sum(y) == p)
    
    for k_idx in 1:k
        @constraint(model, w[k_idx] == sum(y[i] for i in Sk[k_idx]))
    end
    
    for k_idx in 1:k
        @constraint(model, w[k_idx] >= Lk[k_idx])
        @constraint(model, w[k_idx] <= Uk[k_idx])
    end

    set_time_limit_sec(model, 300)
    
    return model, d
end

# Initialize Y using relaxed model approach
function init_y_relaxed(instance, time_limit, log_file)
    start_time = time()
    
    model = build_model_x_relaxed(instance)
    set_optimizer(model, Gurobi.Optimizer)
    set_time_limit_sec(model, time_limit)
    set_optimizer_attribute(model, "LogFile", log_file)
    
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.TIME_LIMIT
        y_sol = value.(model[:y])
        Y = round.(Int, y_sol)
        elapsed_time = time() - start_time
        return Y, elapsed_time, true
    else
        return nothing, time() - start_time, false
    end
end

# Initialize Y using p-dispersion heuristic
function init_y_pdisp(instance, time_limit)
    start_time = time()
    Y = pdisp_heuristic(instance)
    elapsed_time = time() - start_time
    return Y, elapsed_time, true
end

# Initialize Y using multi-p-dispersion
function init_y_multi_pdp(instance, time_limit, log_file)
    start_time = time()
    
    model, d = multi_pdp_min_instance(instance)
    set_optimizer_attribute(model, "LogFile", log_file)
    set_time_limit_sec(model, time_limit)
    
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.TIME_LIMIT
        y_sol = value.(model[:y])
        Y = round.(Int, y_sol)
        elapsed_time = time() - start_time
        return Y, elapsed_time, true
    else
        return nothing, time() - start_time, false
    end
end

# Initialize Y using multi-type p-median
function init_y_multi_p_median(instance, time_limit, log_file, use_benders)
    start_time = time()
    
    # Choose between strong formulation and Benders decomposition
    if use_benders
        Y = solve_multi_type_p_median_benders(instance, time_limit, log_file)
    else
        Y = solve_multi_type_p_median_strong(instance, time_limit, log_file)
    end
    
    elapsed_time = time() - start_time
    
    if Y !== nothing
        return Y, elapsed_time, true
    else
        return nothing, elapsed_time, false
    end
end

# Count constraint violations
function count_violations(instance, x_binary, y_fixed)
    B = instance.B
    S = instance.S
    V = instance.V
    μ = instance.μ
    R = instance.R
    β = instance.β
    m = 3
    ε = instance.T[1]
    
    activity_violations = 0
    risk_violations = 0
    
    for i in 1:S
        if y_fixed[i] == 1
            # Check activity constraints
            for M in 1:m
                activity_level = sum(x_binary[i, j] * V[M][j] for j in 1:B)
                lower_bound = (1-ε) * μ[M][i]
                upper_bound = (1+ε) * μ[M][i]
                
                if activity_level < ceil(lower_bound) || activity_level > floor(upper_bound)
                    activity_violations += 1
                end
            end
            
            # Check risk constraint
            risk_level = sum(x_binary[i, j] * R[j] for j in 1:B)
            if risk_level > β[i]
                risk_violations += 1
            end
        end
    end
    
    return activity_violations, risk_violations
end

# Modified local search to track improvements with configurable moves
function local_search_with_tracking(solution, targets_lower, targets_upper, move_config)
    start_time = time()
    oldSol = solution
    ls_iterations = 0
    successful_simple = 0
    successful_interchange = 0
    successful_deactivate = 0
    initial_weight = solution.Weight
    
    any_improvement = true
    max_iterations = get(move_config, "max_iterations", 100)
    
    while any_improvement && ls_iterations < max_iterations
        ls_iterations += 1
        any_improvement = false
        prev_weight = oldSol.Weight
        
        # Simple move
        if get(move_config, "enable_simple_move", true)
            sol_moved = simple_bu_improve_optimized(oldSol, targets_lower, targets_upper, :bf)
            if sol_moved.Weight < prev_weight
                any_improvement = true
                successful_simple += 1
                prev_weight = sol_moved.Weight
                oldSol = sol_moved
            end
        end
        
        # Interchange move
        if get(move_config, "enable_interchange_move", true)
            sol_interchange = interchange_bu_improve_optimized(oldSol, targets_lower, targets_upper, :bf)
            if sol_interchange.Weight < prev_weight
                any_improvement = true
                successful_interchange += 1
                prev_weight = sol_interchange.Weight
                oldSol = sol_interchange
            end
        end
        
        # Deactivate move
        if get(move_config, "enable_deactivate_move", true)
            sol_deactivate = deactivate_center_improve(oldSol, targets_lower, targets_upper)
            if sol_deactivate.Weight < prev_weight
                any_improvement = true
                successful_deactivate += 1
                prev_weight = sol_deactivate.Weight
                oldSol = sol_deactivate
            end
        end
    end
    
    ls_time = time() - start_time
    total_improvement = initial_weight - oldSol.Weight
    
    return oldSol, ls_time, ls_iterations, total_improvement, successful_simple, successful_interchange, successful_deactivate
end

function repair_activity_and_risk_constraints(instance::Instance, x_binary, y_fixed)
    B = instance.B
    S = instance.S
    D = instance.D
    V = instance.V
    μ = instance.μ
    R = instance.R
    β = instance.β
    m = 3 # activities
    ε = instance.T[1]
    
    # Make a copy of the solution to modify
    x_repaired = copy(x_binary)
    
    # Calculate current activity levels for each service and activity type
    activity_levels = zeros(S, m)
    for i in 1:S, M in 1:m
        activity_levels[i, M] = sum(x_repaired[i, j] * V[M][j] for j in 1:B)
    end
    
    # Calculate current risk levels for each service
    risk_levels = zeros(S)
    for i in 1:S
        risk_levels[i] = sum(x_repaired[i, j] * R[j] for j in 1:B)
    end
    
    # Check which services violate their activity or risk constraints
    activity_violations = Dict{Tuple{Int, Int}, Float64}() # (service, activity) => violation amount
    risk_violations = Dict{Int, Float64}() # service => violation amount
    
    for i in 1:S
        if y_fixed[i] == 1 # Only check active services
            # Check activity constraints
            for M in 1:m
                lower_bound = (1-ε) * μ[M][i]
                upper_bound = (1+ε) * μ[M][i]
                
                if activity_levels[i, M] < ceil(lower_bound)
                    activity_violations[(i, M)] = ceil(lower_bound) - activity_levels[i, M]
                elseif activity_levels[i, M] > floor(upper_bound)
                    activity_violations[(i, M)] = activity_levels[i, M] - floor(upper_bound)
                end
            end
            
            # Check risk constraint
            if risk_levels[i] > β[i]
                risk_violations[i] = risk_levels[i] - β[i]
            end
        end
    end
    
    if isempty(activity_violations) && isempty(risk_violations)
        println("No constraints violations to repair!")
        return x_repaired
    end
    
    println("Found $(length(activity_violations)) activity violations and $(length(risk_violations)) risk violations")
    
    # Maximum number of repair iterations to avoid infinite loops
    max_repair_iterations = 500
    iteration = 0
    
    # Continue until all violations are fixed or max iterations reached
    while ((!isempty(activity_violations) || !isempty(risk_violations)) && iteration < max_repair_iterations)
        iteration += 1
        
        # Determine which type of violation to fix first
        # Prioritize risk violations if they exist
        if !isempty(risk_violations)
            # Sort risk violations by magnitude (largest first)
            sorted_risk_violations = sort(collect(risk_violations), by=x->x[2], rev=true)
            service_i, violation_amount = first(sorted_risk_violations)
            
            println("Iteration $iteration: Fixing risk violation for center $service_i (excess: $violation_amount)")
            
            # Find branches currently assigned to this service, sorted by risk (highest first)
            current_branches = [(j, R[j]) for j in 1:B if x_repaired[service_i, j] == 1]
            sort!(current_branches, by=x->x[2], rev=true)
            
            # Try to reassign high-risk branches to other services
            reassignment_success = false
            
            for (branch_j, branch_risk) in current_branches
                # Find potential services to receive this branch
                receiving_services = [i for i in 1:S if i != service_i && y_fixed[i] == 1]
                
                valid_candidates = []
                for receiving_service in receiving_services
                    # Skip if this would violate risk constraint for receiving service
                    if risk_levels[receiving_service] + branch_risk > β[receiving_service]
                        continue
                    end
                    
                    # Check if reassignment would create/worsen activity violations
                    valid_assignment = true
                    
                    for M in 1:m
                        # Check if the receiving service would violate upper activity bound
                        new_receiver_level = activity_levels[receiving_service, M] + V[M][branch_j]
                        if new_receiver_level > floor((1+ε) * μ[M][receiving_service])
                            valid_assignment = false
                            break
                        end
                        
                        # Check if removing from current service would violate lower activity bound
                        new_current_level = activity_levels[service_i, M] - V[M][branch_j]
                        if new_current_level < ceil((1-ε) * μ[M][service_i])
                            valid_assignment = false
                            break
                        end
                    end
                    
                    if !valid_assignment
                        continue
                    end
                    
                    # Calculate cost impact
                    cost_change = D[receiving_service, branch_j] - D[service_i, branch_j]
                    
                    push!(valid_candidates, (receiving_service, cost_change))
                end
                
                # Sort candidates by cost change (ascending)
                sort!(valid_candidates, by=x->x[2])
                
                if !isempty(valid_candidates)
                    receiving_service, _ = first(valid_candidates)  # Fixed typo here
                    
                    # Reassign the branch
                    x_repaired[service_i, branch_j] = 0
                    x_repaired[receiving_service, branch_j] = 1
                    
                    # Update activity levels
                    for M in 1:m
                        activity_levels[service_i, M] -= V[M][branch_j]
                        activity_levels[receiving_service, M] += V[M][branch_j]
                    end
                    
                    # Update risk levels
                    risk_levels[service_i] -= branch_risk
                    risk_levels[receiving_service] += branch_risk
                    
                    println("  Reassigned BU $branch_j (risk: $branch_risk) from center $service_i to center $receiving_service")
                    reassignment_success = true
                    break
                end
            end
            
            if !reassignment_success
                println("  Could not find a suitable reassignment to reduce risk for center $service_i")
                # Mark as unresolvable to avoid getting stuck
                delete!(risk_violations, service_i)
            end
            
        else
            # No risk violations, handle activity violations
            # Sort violations by magnitude (largest first)
            sorted_activity_violations = sort(collect(activity_violations), by=x->x[2], rev=true)
            (service_i, activity_M), violation_amount = first(sorted_activity_violations)
            
            # Determine if we need to increase or decrease activity
            needs_increase = activity_levels[service_i, activity_M] < ceil((1-ε) * μ[activity_M][service_i])
            
            if needs_increase
                # Need to increase activity - look for branches to add
                println("Iteration $iteration: Fixing low activity violation for center $service_i, activity $activity_M")
                
                # Find branches currently assigned to other services
                other_service_branches = [(j, i) for j in 1:B, i in 1:S if x_repaired[i, j] == 1 && i != service_i]
                
                # Calculate the benefit of reassigning each branch
                candidates = []
                for (branch_j, current_service) in other_service_branches
                    # Skip if this would violate risk constraint
                    if risk_levels[service_i] + R[branch_j] > β[service_i]
                        continue
                    end
                    
                    # Check if adding this branch would cause violations for any activity measure
                    valid_assignment = true
                    for check_M in 1:m
                        new_level = activity_levels[service_i, check_M] + V[check_M][branch_j]
                        if new_level > floor((1+ε) * μ[check_M][service_i])
                            valid_assignment = false
                            break
                        end
                    end
                    
                    if !valid_assignment
                        continue
                    end
                    
                    # Check if removing would cause a violation for the current service
                    for check_M in 1:m
                        new_level = activity_levels[current_service, check_M] - V[check_M][branch_j]
                        if new_level < ceil((1-ε) * μ[check_M][current_service])
                            valid_assignment = false
                            break
                        end
                    end
                    
                    if !valid_assignment
                        continue
                    end
                    
                    # Check if removing would cause a risk violation for the current service
                    # This shouldn't be possible since we're lowering risk, but included for completeness
                    
                    # Calculate activity gain for the target activity
                    activity_gain = V[activity_M][branch_j]
                    
                    # Calculate cost impact
                    cost_change = D[service_i, branch_j] - D[current_service, branch_j]
                    
                    push!(candidates, (branch_j, current_service, activity_gain, cost_change))
                end
                
                # Sort candidates by activity gain (descending) and then by cost change (ascending)
                sort!(candidates, by=x->(x[3], -x[4]), rev=true)
                
                # Try the best candidate
                if !isempty(candidates)
                    branch_j, current_service, gain, _ = first(candidates)
                    
                    # Reassign the branch
                    x_repaired[service_i, branch_j] = 1
                    x_repaired[current_service, branch_j] = 0
                    
                    # Update activity levels for ALL activity types
                    for M in 1:m
                        activity_levels[service_i, M] += V[M][branch_j]
                        activity_levels[current_service, M] -= V[M][branch_j]
                    end
                    
                    # Update risk levels
                    risk_levels[service_i] += R[branch_j]
                    risk_levels[current_service] -= R[branch_j]
                    
                    println("  Reassigned BU $branch_j from center $current_service to center $service_i to increase activity $activity_M by $gain")
                else
                    println("  Could not find a suitable BU to increase activity $activity_M for center $service_i")
                    # Mark this violation as unresolvable for now to avoid getting stuck
                    delete!(activity_violations, (service_i, activity_M))
                end
            else
                # Need to decrease activity - look for branches to remove
                println("Iteration $iteration: Fixing high activity violation for center $service_i, activity $activity_M")
                
                # Find branches currently assigned to this service
                current_branches = [(j, V[activity_M][j]) for j in 1:B if x_repaired[service_i, j] == 1]
                # Sort by activity contribution (highest first)
                sort!(current_branches, by=x->x[2], rev=true)
                
                # Try to reassign high-activity branches to other services
                reassignment_success = false
                
                for (branch_j, branch_activity) in current_branches
                    # Find potential services to receive this branch
                    receiving_services = [i for i in 1:S if i != service_i && y_fixed[i] == 1]
                    
                    valid_candidates = []
                    for receiving_service in receiving_services
                        # Skip if this would violate risk constraint for receiving service
                        if risk_levels[receiving_service] + R[branch_j] > β[receiving_service]
                            continue
                        end
                        
                        # Check if reassignment would create/worsen activity violations
                        valid_assignment = true
                        
                        for M in 1:m
                            # Check if the receiving service would violate upper activity bound
                            new_receiver_level = activity_levels[receiving_service, M] + V[M][branch_j]
                            if new_receiver_level > floor((1+ε) * μ[M][receiving_service])
                                valid_assignment = false
                                break
                            end
                        end
                        
                        if !valid_assignment
                            continue
                        end
                        
                        # Calculate cost impact
                        cost_change = D[receiving_service, branch_j] - D[service_i, branch_j]
                        
                        push!(valid_candidates, (receiving_service, cost_change))
                    end
                    
                    # Sort candidates by cost change (ascending)
                    sort!(valid_candidates, by=x->x[2])
                    
                    if !isempty(valid_candidates)
                        receiving_service, _ = first(valid_candidates)  # Fixed typo here
                        
                        # Reassign the branch
                        x_repaired[service_i, branch_j] = 0
                        x_repaired[receiving_service, branch_j] = 1
                        
                        # Update activity levels
                        for M in 1:m
                            activity_levels[service_i, M] -= V[M][branch_j]
                            activity_levels[receiving_service, M] += V[M][branch_j]
                        end
                        
                        # Update risk levels
                        risk_levels[service_i] -= R[branch_j]
                        risk_levels[receiving_service] += R[branch_j]
                        
                        println("  Reassigned BU $branch_j from center $service_i to center $receiving_service to decrease activity $activity_M by $branch_activity")
                        reassignment_success = true
                        break
                    end
                end
                
                if !reassignment_success
                    println("  Could not find a suitable reassignment to decrease activity $activity_M for center $service_i")
                    # Mark as unresolvable to avoid getting stuck
                    delete!(activity_violations, (service_i, activity_M))
                end
            end
        end
        
        # Recalculate violations after the reassignment
        activity_violations = Dict{Tuple{Int, Int}, Float64}()
        risk_violations = Dict{Int, Float64}()
        
        for i in 1:S
            if y_fixed[i] == 1
                # Check activity constraints
                for M in 1:m
                    lower_bound = (1-ε) * μ[M][i]
                    upper_bound = (1+ε) * μ[M][i]
                    
                    if activity_levels[i, M] < ceil(lower_bound)
                        activity_violations[(i, M)] = ceil(lower_bound) - activity_levels[i, M]
                    elseif activity_levels[i, M] > floor(upper_bound)
                        activity_violations[(i, M)] = activity_levels[i, M] - floor(upper_bound)
                    end
                end
                
                # Check risk constraint
                if risk_levels[i] > β[i]
                    risk_violations[i] = risk_levels[i] - β[i]
                end
            end
        end
    end
    
    # Final report on violations
    if isempty(activity_violations) && isempty(risk_violations)
        println("Successfully repaired all constraints violations in $iteration iterations!")
    else
        println("Could not repair all violations after $iteration iterations:")
        println("  Remaining activity violations: $(length(activity_violations))")
        println("  Remaining risk violations: $(length(risk_violations))")
    end
    
    return x_repaired
end