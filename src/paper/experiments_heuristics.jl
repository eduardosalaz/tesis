using Gurobi, JuMP, JLD2
using Types
using MathOptInterface
const MOI = MathOptInterface
using Dates, TOML, CSV, JSON, DataFrames, Statistics
using LinearAlgebra, Distances, Plots
using Printf
include("generator.jl")
include("solver.jl") 
include("../heuristics/constructive.jl")
include("../math_prog/first_model.jl")
include("../heuristics/ls.jl")

# Include the model building and helper functions from the heuristic code
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
    set_optimizer_attribute(master, "MIPGap", 1e-4)
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

# Data structures for tracking results
struct ExperimentResult
    instance_id::String
    method::String
    move_config::String  # New field to track which moves were enabled
    B::Int
    S::Int
    P::Int
    # Timing
    init_time::Float64
    phase2_time::Float64
    repair_time::Float64
    ls_time::Float64
    total_time::Float64
    # Objective values
    phase2_obj::Float64
    split_resolved_obj::Float64
    post_repair_obj::Float64
    post_ls_obj::Float64
    # Split information
    num_split_bus::Int
    prop_split_bus::Float64
    # Violations
    activity_violations_initial::Int
    risk_violations_initial::Int
    repair_iterations::Int
    remaining_violations::Int
    # Local search
    ls_iterations::Int
    total_ls_improvement::Float64
    successful_moves_simple::Int
    successful_moves_interchange::Int
    successful_moves_deactivate::Int
    # Warmstart (optional)
    warmstart_used::Bool
    warmstart_time::Float64
    warmstart_obj::Float64
    # Status
    status::String
    error_message::String
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

# Modified repair function to track iterations
function repair_with_tracking(instance, x_binary, y_fixed)
    start_time = time()
    x_repaired = copy(x_binary)
    max_iterations = 200
    iteration = 0
    
    # Initial violation count
    initial_activity, initial_risk = count_violations(instance, x_binary, y_fixed)
    
    # Run repair (simplified version - use your full repair function)
    x_repaired = repair_activity_and_risk_constraints(instance, x_binary, y_fixed)
    
    # Count remaining violations
    final_activity, final_risk = count_violations(instance, x_repaired, y_fixed)
    
    repair_time = time() - start_time
    
    return x_repaired, repair_time, iteration, initial_activity + initial_risk, final_activity + final_risk
end

# Structure to hold construction phase results
struct ConstructionResult
    solution::Solution
    stats::NamedTuple
    targets_lower::Any
    targets_upper::Any
end

# Run only construction phases (no local search) - WITH ROBUST ERROR HANDLING
function run_construction_phases(instance, instance_file, method, config, output_dirs, experiment_log)
    instance_name = splitext(basename(instance_file))[1]
    B, S, P = instance.B, instance.S, instance.P
    size_str = "$(B)_$(S)_$(P)"
    
    println("    Running construction phases...")
    
    total_start_time = time()
    
    try
        # Phase 1: Initialize Y
        y_log_file = joinpath(output_dirs["logs"], "y_init_$(instance_name)_$(method).log")
        
        if method == "relaxed"
            Y, init_time, success = init_y_relaxed(instance, config["time_limits"]["exact_methods"], y_log_file)
        elseif method == "pdisp"
            Y, init_time, success = init_y_pdisp(instance, config["time_limits"]["exact_methods"])
        elseif method == "multi_pdp"
            Y, init_time, success = init_y_multi_pdp(instance, config["time_limits"]["exact_methods"], y_log_file)
        elseif method == "multi_p_median"
            Y, init_time, success = init_y_multi_p_median(instance, config["time_limits"]["exact_methods"], y_log_file, false)
        elseif method == "multi_p_median_benders"
            Y, init_time, success = init_y_multi_p_median(instance, config["time_limits"]["exact_methods"], y_log_file, true)
        else
            error("Unknown method: $method")
        end
        
        if !success || Y === nothing
            println("      Y initialization failed")
            return nothing
        end
        
        println("      Y initialized in $(round(init_time, digits=2))s")
        
        # Phase 2: Solve transportation problem
        phase2_start = time()
        phase2_log_file = joinpath(output_dirs["logs"], "phase2_$(instance_name)_$(method).log")
        
        model_transport = solve_phase2_model(instance, Y)
        set_optimizer(model_transport, Gurobi.Optimizer)
        set_optimizer_attribute(model_transport, "LogFile", phase2_log_file)
        set_time_limit_sec(model_transport, config["time_limits"]["phase2"])
        
        optimize!(model_transport)
        
        if termination_status(model_transport) != MOI.OPTIMAL && termination_status(model_transport) != MOI.TIME_LIMIT
            println("      Phase 2 optimization failed")
            return nothing
        end
        
        x_continuous = value.(model_transport[:x])
        phase2_obj = objective_value(model_transport)
        phase2_time = time() - phase2_start
        
        println("      Phase 2 solved in $(round(phase2_time, digits=2))s, obj: $(round(phase2_obj, digits=2))")
        
        # Analyze splits
        num_splits, prop_splits, _ = analyze_splits(x_continuous, instance.S, instance.B)
        println("      Splits: $num_splits BUs ($(round(prop_splits * 100, digits=1))%)")
        
        # Split resolution
        x_binary = split_resolution_heuristic(x_continuous, instance.S, instance.B)
        split_resolved_obj = dot(x_binary, instance.D)
        
        # Count initial violations
        initial_activity_viol, initial_risk_viol = count_violations(instance, x_binary, Y)
        
        # Repair constraints - WRAPPED IN TRY-CATCH
        repair_start = time()
        x_repaired = x_binary  # Default to unrepaired if repair fails
        repair_time = 0.0
        try
            x_repaired = repair_activity_and_risk_constraints(instance, x_binary, Y)
            repair_time = time() - repair_start
        catch e
            println("      WARNING: Repair failed with error: $(sprint(showerror, e))")
            println("      Continuing with unrepaired solution...")
            repair_time = time() - repair_start
        end
        
        # Count remaining violations
        final_activity_viol, final_risk_viol = count_violations(instance, x_repaired, Y)
        repair_iterations = 0  # Could track this in repair function if needed
        
        post_repair_obj = dot(x_repaired, instance.D)
        
        println("      Repair completed in $(round(repair_time, digits=2))s, obj: $(round(post_repair_obj, digits=2))")
        
        # Create solution object
        solution = Solution(instance, x_repaired, Y, post_repair_obj, 1)
        
        # Calculate targets for local search
        targets_lower, targets_upper = calculate_targets_optimized(instance)
        
        # Create stats tuple
        stats = (
            init_time = init_time,
            phase2_time = phase2_time,
            repair_time = repair_time,
            phase2_obj = phase2_obj,
            split_resolved_obj = split_resolved_obj,
            post_repair_obj = post_repair_obj,
            num_split_bus = num_splits,
            prop_split_bus = prop_splits,
            activity_violations_initial = initial_activity_viol,
            risk_violations_initial = initial_risk_viol,
            repair_iterations = repair_iterations,
            remaining_violations = final_activity_viol + final_risk_viol
        )
        
        # Log construction results
        log_message = """
        Construction completed for $instance_name with $method:
          Y init time: $(round(init_time, digits=2))s
          Phase 2 obj: $phase2_obj (time: $(round(phase2_time, digits=2))s)
          Split BUs: $(round(prop_splits * 100, digits=2))%
          Split resolved obj: $split_resolved_obj  
          Post repair obj: $post_repair_obj (time: $(round(repair_time, digits=2))s)
          Initial violations: $(initial_activity_viol + initial_risk_viol)
          Remaining violations: $(final_activity_viol + final_risk_viol)
        """
        open(experiment_log, "a") do f
            write(f, log_message * "\n")
        end
        
        return ConstructionResult(solution, stats, targets_lower, targets_upper)
        
    catch e
        error_msg = sprint(showerror, e)
        println("      ERROR in construction: $error_msg")
        
        open(experiment_log, "a") do f
            write(f, "ERROR in construction for $instance_name with $method: $error_msg\n")
        end
        
        return nothing
    end
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

# Define move configurations for experiments
function get_move_configurations()
    # Only testing 4 specific combinations as requested
    return Dict(
        "all_moves" => Dict(
            "enable_simple_move" => true,
            "enable_interchange_move" => true,
            "enable_deactivate_move" => true,
            "max_iterations" => 100
        ),
        "simple_interchange" => Dict(
            "enable_simple_move" => true,
            "enable_interchange_move" => true,
            "enable_deactivate_move" => false,
            "max_iterations" => 100
        ),
        "simple_deactivate" => Dict(
            "enable_simple_move" => true,
            "enable_interchange_move" => false,
            "enable_deactivate_move" => true,
            "max_iterations" => 100
        ),
        "interchange_deactivate" => Dict(
            "enable_simple_move" => false,
            "enable_interchange_move" => true,
            "enable_deactivate_move" => true,
            "max_iterations" => 100
        )
    )
end

# Main experiment runner with efficient construction-LS separation - FULLY ROBUST
function run_heuristic_experiments(config_file="config.toml")
    # Load configuration
    config = TOML.parsefile(config_file)
    
    # Create experiment directory
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    experiment_name = get(config["general"], "experiment_name", "heuristic_comparison")
    experiment_id = "exp_$(timestamp)_$(experiment_name)"
    experiment_dir = joinpath(config["general"]["base_path"], "experiments", experiment_id)
    
    # Create directory structure
    output_dirs = Dict(
        "experiment" => experiment_dir,
        "results" => joinpath(experiment_dir, "results"),
        "instances" => joinpath(experiment_dir, "instances"),
        "solutions" => joinpath(experiment_dir, "solutions"),
        "pre_ls_solutions" => joinpath(experiment_dir, "pre_ls_solutions"),  # New directory for pre-LS solutions
        "plots" => joinpath(experiment_dir, "plots"),
        "logs" => joinpath(experiment_dir, "logs")
    )
    
    for dir in values(output_dirs)
        mkpath(dir)
    end
    
    # Copy config
    cp(config_file, joinpath(experiment_dir, "config.toml"))
    
    # Create experiment log
    experiment_log = joinpath(experiment_dir, "experiment_log.txt")
    open(experiment_log, "w") do f
        write(f, """
        Heuristic Comparison Experiment with Move Configurations
        ID: $experiment_id
        Started: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        Configuration: $config_file
        
        """)
    end
    
    # Get list of instance files
    input_dir = config["general"]["input_directory"]
    instance_files = String[]
    
    for (root, dirs, files) in walkdir(input_dir)
        for file in files
            if endswith(file, ".jld2")
                push!(instance_files, joinpath(root, file))
            end
        end
    end
    
    if isempty(instance_files)
        error("No instance files found in $input_dir")
    end
    
    sort!(instance_files)
    println("Found $(length(instance_files)) instance files")
    
    # Determine which methods to test
    methods_to_test = String[]
    if get(config["methods"], "test_relaxed", true)
        push!(methods_to_test, "relaxed")
    end
    if get(config["methods"], "test_pdisp", true)
        push!(methods_to_test, "pdisp")
    end
    if get(config["methods"], "test_multi_pdp", true)
        push!(methods_to_test, "multi_pdp")
    end
    if get(config["methods"], "test_multi_p_median", false)
        push!(methods_to_test, "multi_p_median")
    end
    if get(config["methods"], "test_multi_p_median_benders", false)
        push!(methods_to_test, "multi_p_median_benders")
    end
    
    println("Testing methods: ", join(methods_to_test, ", "))
    
    # Get move configurations
    move_configs = get_move_configurations()
    
    # Determine which move configurations to test
    test_move_ablation = get(config["local_search"], "test_move_ablation", true)
    
    if test_move_ablation
        configs_to_test = collect(keys(move_configs))
        println("Testing move configurations: ", join(configs_to_test, ", "))
    else
        configs_to_test = ["all_moves"]
        println("Testing only with all moves enabled")
    end
    
    # Results storage
    all_results = ExperimentResult[]
    
    # PHASE 1: Build and save all pre-LS solutions
    println("\n" * "="^70)
    println("PHASE 1: CONSTRUCTIVE HEURISTICS")
    println("="^70)
    
    pre_ls_solutions = Dict{Tuple{String, String}, Any}()  # (instance_name, method) => solution data
    
    for (idx, instance_file) in enumerate(instance_files)
        instance_name = splitext(basename(instance_file))[1]
        println("\n=== Building pre-LS solutions for instance $idx/$(length(instance_files)): $instance_name ===")
        
        # Load instance - WRAPPED IN TRY-CATCH
        instance = nothing
        try
            instance = read_instance(instance_file)
        catch e
            println("ERROR: Failed to load instance $instance_name: $(sprint(showerror, e))")
            continue  # Skip to next instance
        end
        
        # Copy instance file
        size_str = "$(instance.B)_$(instance.S)_$(instance.P)"
        inst_dir = joinpath(output_dirs["instances"], size_str)
        mkpath(inst_dir)
        cp(instance_file, joinpath(inst_dir, basename(instance_file)))
        
        # Build pre-LS solution for each method
        for method in methods_to_test
            println("\n  Method: $method")
            
            # Run construction phases (once per method-instance pair) - ROBUST TO FAILURES
            construction_result = nothing
            try
                construction_result = run_construction_phases(
                    instance, instance_file, method, config, output_dirs, experiment_log
                )
            catch e
                println("    ERROR in construction: $(sprint(showerror, e))")
                println("    Skipping this method for instance $instance_name")
                continue  # Skip to next method
            end
            
            if construction_result !== nothing
                # Save the pre-LS solution
                pre_ls_sol_dir = joinpath(output_dirs["pre_ls_solutions"], size_str, method)
                mkpath(pre_ls_sol_dir)
                pre_ls_sol_file = joinpath(pre_ls_sol_dir, "pre_ls_$(instance_name)_$(method).jld2")
                
                # Save all necessary data for local search
                try
                    JLD2.save(pre_ls_sol_file,
                        "solution", construction_result.solution,
                        "construction_stats", construction_result.stats,
                        "targets_lower", construction_result.targets_lower,
                        "targets_upper", construction_result.targets_upper
                    )
                    
                    pre_ls_solutions[(instance_name, method)] = construction_result
                    
                    println("    Pre-LS objective: $(construction_result.stats.post_repair_obj)")
                    println("    Saved to: $(basename(pre_ls_sol_file))")
                catch e
                    println("    ERROR saving pre-LS solution: $(sprint(showerror, e))")
                end
            else
                println("    FAILED - Construction phases failed")
            end
        end
    end
    
    # PHASE 2: Apply local search configurations to saved solutions
    println("\n" * "="^70)
    println("PHASE 2: LOCAL SEARCH ABLATION")
    println("="^70)
    
    for (idx, instance_file) in enumerate(instance_files)
        instance_name = splitext(basename(instance_file))[1]
        
        # Load instance - WRAPPED IN TRY-CATCH
        instance = nothing
        try
            instance = read_instance(instance_file)
        catch e
            println("ERROR: Failed to load instance $instance_name: $(sprint(showerror, e))")
            continue  # Skip to next instance
        end
        
        size_str = "$(instance.B)_$(instance.S)_$(instance.P)"
        
        println("\n=== Applying LS to instance $idx/$(length(instance_files)): $instance_name ===")
        
        for method in methods_to_test
            # Check if we have a pre-LS solution for this instance-method pair
            if haskey(pre_ls_solutions, (instance_name, method))
                construction_data = pre_ls_solutions[(instance_name, method)]
                
                println("\n  Method: $method (pre-LS obj: $(construction_data.stats.post_repair_obj))")
                
                # Apply each move configuration to the SAME pre-LS solution
                for config_name in configs_to_test
                    println("    Testing moves: $config_name")
                    
                    try
                        # Create a copy of the solution for this LS configuration
                        sol_copy = deepcopy(construction_data.solution)
                        
                        # Run local search with specific move configuration
                        sol_after_ls, ls_time, ls_iters, total_improvement, succ_simple, succ_interchange, succ_deactivate = 
                            local_search_with_tracking(
                                sol_copy, 
                                construction_data.targets_lower, 
                                construction_data.targets_upper, 
                                move_configs[config_name]
                            )
                        
                        post_ls_obj = sol_after_ls.Weight
                        
                        # Save the final solution
                        sol_dir = joinpath(output_dirs["solutions"], size_str, method, config_name)
                        mkpath(sol_dir)
                        sol_file = joinpath(sol_dir, "sol_$(instance_name)_$(method)_$(config_name).jld2")
                        write_solution(sol_after_ls, sol_file)
                        
                        # Create result combining construction stats with LS results
                        result = ExperimentResult(
                            instance_name, method, config_name, instance.B, instance.S, instance.P,
                            construction_data.stats.init_time,
                            construction_data.stats.phase2_time,
                            construction_data.stats.repair_time,
                            ls_time,
                            construction_data.stats.init_time + construction_data.stats.phase2_time + 
                                construction_data.stats.repair_time + ls_time,
                            construction_data.stats.phase2_obj,
                            construction_data.stats.split_resolved_obj,
                            construction_data.stats.post_repair_obj,
                            post_ls_obj,
                            construction_data.stats.num_split_bus,
                            construction_data.stats.prop_split_bus,
                            construction_data.stats.activity_violations_initial,
                            construction_data.stats.risk_violations_initial,
                            construction_data.stats.repair_iterations,
                            construction_data.stats.remaining_violations,
                            ls_iters, total_improvement, succ_simple, succ_interchange, succ_deactivate,
                            false, 0.0, post_ls_obj,
                            "success", ""
                        )
                        
                        push!(all_results, result)
                        
                        println("      Post-LS obj: $post_ls_obj (improvement: $total_improvement)")
                        
                    catch e
                        println("      ERROR in local search: $(sprint(showerror, e))")
                        
                        # Create failed result for this configuration
                        result = ExperimentResult(
                            instance_name, method, config_name, instance.B, instance.S, instance.P,
                            construction_data.stats.init_time,
                            construction_data.stats.phase2_time,
                            construction_data.stats.repair_time,
                            0.0, 0.0,
                            construction_data.stats.phase2_obj,
                            construction_data.stats.split_resolved_obj,
                            construction_data.stats.post_repair_obj,
                            construction_data.stats.post_repair_obj,  # No LS improvement
                            construction_data.stats.num_split_bus,
                            construction_data.stats.prop_split_bus,
                            construction_data.stats.activity_violations_initial,
                            construction_data.stats.risk_violations_initial,
                            construction_data.stats.repair_iterations,
                            construction_data.stats.remaining_violations,
                            0, 0.0, 0, 0, 0,
                            false, 0.0, 0.0,
                            "ls_failed", sprint(showerror, e)
                        )
                        push!(all_results, result)
                    end
                end
            else
                # Construction failed for this instance-method pair
                # Create failed results for all move configurations
                for config_name in configs_to_test
                    result = ExperimentResult(
                        instance_name, method, config_name, instance.B, instance.S, instance.P,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0,
                        0, 0, 0, 0, 0, 0.0, 0, 0, 0,
                        false, 0.0, 0.0,
                        "failed_construction", "Construction phases failed"
                    )
                    push!(all_results, result)
                end
            end
        end
        
        # Save intermediate results after each instance
        try
            save_results(all_results, output_dirs["results"])
        catch e
            println("WARNING: Failed to save intermediate results: $(sprint(showerror, e))")
        end
    end
    
    # Generate final summary
    try
        generate_summary_with_moves(all_results, output_dirs["results"], experiment_log)
    catch e
        println("WARNING: Failed to generate summary: $(sprint(showerror, e))")
    end
    
    # Log completion
    completion_msg = "\n=== Experiment Complete at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) ==="
    println(completion_msg)
    open(experiment_log, "a") do f
        write(f, completion_msg * "\n")
    end
    
    println("\nResults saved to: $experiment_dir")
end

# Save results to CSV and JSON - Rest of the functions remain the same...
function save_results(results::Vector{ExperimentResult}, results_dir::String)
    # Convert to DataFrame for CSV
    df = DataFrame(
        instance_id = [r.instance_id for r in results],
        method = [r.method for r in results],
        move_config = [r.move_config for r in results],
        B = [r.B for r in results],
        S = [r.S for r in results],
        P = [r.P for r in results],
        init_time = [r.init_time for r in results],
        phase2_time = [r.phase2_time for r in results],
        repair_time = [r.repair_time for r in results],
        ls_time = [r.ls_time for r in results],
        total_time = [r.total_time for r in results],
        phase2_obj = [r.phase2_obj for r in results],
        split_resolved_obj = [r.split_resolved_obj for r in results],
        post_repair_obj = [r.post_repair_obj for r in results],
        post_ls_obj = [r.post_ls_obj for r in results],
        num_split_bus = [r.num_split_bus for r in results],
        prop_split_bus = [r.prop_split_bus for r in results],
        activity_violations_initial = [r.activity_violations_initial for r in results],
        risk_violations_initial = [r.risk_violations_initial for r in results],
        repair_iterations = [r.repair_iterations for r in results],
        remaining_violations = [r.remaining_violations for r in results],
        ls_iterations = [r.ls_iterations for r in results],
        total_ls_improvement = [r.total_ls_improvement for r in results],
        successful_moves_simple = [r.successful_moves_simple for r in results],
        successful_moves_interchange = [r.successful_moves_interchange for r in results],
        successful_moves_deactivate = [r.successful_moves_deactivate for r in results],
        warmstart_used = [r.warmstart_used for r in results],
        warmstart_time = [r.warmstart_time for r in results],
        warmstart_obj = [r.warmstart_obj for r in results],
        status = [r.status for r in results],
        error_message = [r.error_message for r in results]
    )
    
    CSV.write(joinpath(results_dir, "detailed_results.csv"), df)
    
    # Save as JSON
    json_results = [Dict(
        "instance_id" => r.instance_id,
        "method" => r.method,
        "move_config" => r.move_config,
        "parameters" => Dict("B" => r.B, "S" => r.S, "P" => r.P),
        "timing" => Dict(
            "init" => r.init_time,
            "phase2" => r.phase2_time,
            "repair" => r.repair_time,
            "local_search" => r.ls_time,
            "total" => r.total_time
        ),
        "objectives" => Dict(
            "phase2" => r.phase2_obj,
            "split_resolved" => r.split_resolved_obj,
            "post_repair" => r.post_repair_obj,
            "post_ls" => r.post_ls_obj
        ),
        "local_search" => Dict(
            "iterations" => r.ls_iterations,
            "total_improvement" => r.total_ls_improvement,
            "successful_moves" => Dict(
                "simple" => r.successful_moves_simple,
                "interchange" => r.successful_moves_interchange,
                "deactivate" => r.successful_moves_deactivate
            )
        ),
        "status" => r.status
    ) for r in results]
    
    open(joinpath(results_dir, "detailed_results.json"), "w") do f
        JSON.print(f, json_results, 2)
    end
end

# Generate summary with move analysis
function generate_summary_with_moves(results::Vector{ExperimentResult}, results_dir::String, log_file::String)
    # Group by method and move configuration
    methods = unique([r.method for r in results])
    move_configs = unique([r.move_config for r in results])
    
    summary_df = DataFrame()
    
    for method in methods
        for move_config in move_configs
            filtered_results = filter(r -> r.method == method && r.move_config == move_config && r.status == "success", results)
            
            if !isempty(filtered_results)
                summary_row = DataFrame(
                    method = method,
                    move_config = move_config,
                    num_instances = length(filtered_results),
                    avg_total_time = mean([r.total_time for r in filtered_results]),
                    avg_post_repair_obj = mean([r.post_repair_obj for r in filtered_results]),
                    avg_post_ls_obj = mean([r.post_ls_obj for r in filtered_results]),
                    avg_ls_improvement = mean([r.total_ls_improvement for r in filtered_results]),
                    avg_ls_iterations = mean([r.ls_iterations for r in filtered_results]),
                    avg_simple_moves = mean([r.successful_moves_simple for r in filtered_results]),
                    avg_interchange_moves = mean([r.successful_moves_interchange for r in filtered_results]),
                    avg_deactivate_moves = mean([r.successful_moves_deactivate for r in filtered_results])
                )
                summary_df = vcat(summary_df, summary_row)
            end
        end
    end
    
    CSV.write(joinpath(results_dir, "summary_by_method_and_moves.csv"), summary_df)
    
    # Create move contribution analysis
    if "all_moves" in move_configs
        move_contribution_df = DataFrame()
        
        for method in methods
            baseline_results = filter(r -> r.method == method && r.move_config == "all_moves" && r.status == "success", results)
            
            if !isempty(baseline_results)
                baseline_obj = mean([r.post_ls_obj for r in baseline_results])
                baseline_time = mean([r.ls_time for r in baseline_results])
                
                for move_config in move_configs
                    if move_config != "all_moves"
                        config_results = filter(r -> r.method == method && r.move_config == move_config && r.status == "success", results)
                        
                        if !isempty(config_results)
                            config_obj = mean([r.post_ls_obj for r in config_results])
                            config_time = mean([r.ls_time for r in config_results])
                            
                            contribution_row = DataFrame(
                                method = method,
                                move_config = move_config,
                                obj_degradation = config_obj - baseline_obj,
                                obj_degradation_pct = (config_obj - baseline_obj) / baseline_obj * 100,
                                time_saved = baseline_time - config_time,
                                time_saved_pct = (baseline_time - config_time) / baseline_time * 100
                            )
                            move_contribution_df = vcat(move_contribution_df, contribution_row)
                        end
                    end
                end
            end
        end
        
        if !isempty(move_contribution_df)
            CSV.write(joinpath(results_dir, "move_contribution_analysis.csv"), move_contribution_df)
        end
    end
    
    # Write summary to log
    summary_text = "\n\n=== SUMMARY STATISTICS ===\n"
    for method in methods
        summary_text *= "\n$method:\n"
        for move_config in move_configs
            method_config_results = filter(r -> r.method == method && r.move_config == move_config && r.status == "success", results)
            if !isempty(method_config_results)
                avg_obj = mean([r.post_ls_obj for r in method_config_results])
                avg_improvement = mean([r.total_ls_improvement for r in method_config_results])
                avg_time = mean([r.ls_time for r in method_config_results])
                summary_text *= "  $move_config: obj=$(round(avg_obj, digits=2)), improvement=$(round(avg_improvement, digits=2)), time=$(round(avg_time, digits=2))s\n"
            end
        end
    end
    
    println(summary_text)
    open(log_file, "a") do f
        write(f, summary_text)
    end
end

# Main function
function main()
    if length(ARGS) == 0
        println("Using default config file: config.toml")
        config_file = "config.toml"
    else
        config_file = ARGS[1]
    end
    
    run_heuristic_experiments(config_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end