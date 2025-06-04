using Gurobi, JuMP, JLD2
using Types
using MathOptInterface
const MOI = MathOptInterface
using Dates, TOML, CSV, JSON, DataFrames, Statistics
using LinearAlgebra, Distances, Plots
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

# The full repair_activity_and_risk_constraints function would go here
# For brevity, I'm including a simplified version
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
    max_repair_iterations = 200
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
                    receiving_service, _ = first(valid_candidates)
                    
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
                        receiving_service, _ = first(valid_candidates)
                        
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
    
    # Verify all constraints are satisfied
    # 1. Final check of activity constraints
    activity_satisfied = true
    for i in 1:S
        if y_fixed[i] == 1
            for M in 1:m
                lower_bound = (1-ε) * μ[M][i]
                upper_bound = (1+ε) * μ[M][i]
                
                if activity_levels[i, M] < ceil(lower_bound) || activity_levels[i, M] > floor(upper_bound)
                    activity_satisfied = false
                    println("  Activity constraint violated: Service $i, Activity $M, Level: $(activity_levels[i, M]), Bounds: [$(ceil(lower_bound)), $(floor(upper_bound))]")
                end
            end
        end
    end
    
    # 2. Final check of risk constraints
    risk_satisfied = true
    for i in 1:S
        if y_fixed[i] == 1
            if risk_levels[i] > β[i]
                risk_satisfied = false
                println("  Risk constraint violated: Service $i, Risk: $(risk_levels[i]), Threshold: $(β[i])")
            end
        end
    end
    
    if activity_satisfied && risk_satisfied
        println("All constraints are satisfied in the final solution!")
    else
        println("The final solution still has constraint violations.")
    end
    
    return x_repaired
end

# Data structures for tracking results
struct ExperimentResult
    instance_id::String
    method::String
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

# Modified local search to track improvements
function local_search_with_tracking(solution, targets_lower, targets_upper, config)
    start_time = time()
    oldSol = solution
    ls_iterations = 0
    successful_simple = 0
    successful_interchange = 0
    successful_deactivate = 0
    initial_weight = solution.Weight
    
    any_improvement = true
    max_iterations = get(config["local_search"], "max_iterations", 100)
    
    while any_improvement && ls_iterations < max_iterations
        ls_iterations += 1
        any_improvement = false
        prev_weight = oldSol.Weight
        
        # Simple move
        if get(config["local_search"], "enable_simple_move", true)
            sol_moved = simple_bu_improve_optimized(oldSol, targets_lower, targets_upper, :bf)
            if sol_moved.Weight < prev_weight
                any_improvement = true
                successful_simple += 1
                prev_weight = sol_moved.Weight
                oldSol = sol_moved
            end
        end
        
        # Interchange move
        if get(config["local_search"], "enable_interchange_move", true)
            sol_interchange = interchange_bu_improve_optimized(oldSol, targets_lower, targets_upper, :bf)
            if sol_interchange.Weight < prev_weight
                any_improvement = true
                successful_interchange += 1
                prev_weight = sol_interchange.Weight
                oldSol = sol_interchange
            end
        end
        
        # Deactivate move
        if get(config["local_search"], "enable_deactivate_move", true)
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

# Process single instance with one method
function process_instance_method(instance, instance_file, method, config, output_dirs, experiment_log)
    instance_name = splitext(basename(instance_file))[1]
    B, S, P = instance.B, instance.S, instance.P
    size_str = "$(B)_$(S)_$(P)"
    
    println("\n--- Processing $instance_name with method: $method ---")
    
    # Initialize result
    result = ExperimentResult(
        instance_name, method, B, S, P,
        0.0, 0.0, 0.0, 0.0, 0.0,  # timing
        0.0, 0.0, 0.0, 0.0,       # objectives
        0, 0.0,                    # splits
        0, 0, 0, 0,               # violations
        0, 0.0, 0, 0, 0,          # local search
        false, 0.0, 0.0,          # warmstart
        "failed", ""              # status
    )
    
    total_start_time = time()
    
    try
        # Phase 1: Initialize Y
        y_log_file = joinpath(output_dirs["logs"], "y_init_$(instance_name)_$(method).log")
        
        if method == "relaxed"
            Y, init_time, success = init_y_relaxed(instance, config["time_limits"]["y_init_relaxed"], y_log_file)
        elseif method == "pdisp"
            Y, init_time, success = init_y_pdisp(instance, config["time_limits"]["y_init_pdisp"])
        elseif method == "multi_pdp"
            Y, init_time, success = init_y_multi_pdp(instance, config["time_limits"]["y_init_multi_pdp"], y_log_file)
        else
            error("Unknown method: $method")
        end
        
        if !success || Y === nothing
            result = ExperimentResult(result.instance_id, result.method, result.B, result.S, result.P,
                init_time, result.phase2_time, result.repair_time, result.ls_time, init_time,
                result.phase2_obj, result.split_resolved_obj, result.post_repair_obj, result.post_ls_obj,
                result.num_split_bus, result.prop_split_bus,
                result.activity_violations_initial, result.risk_violations_initial,
                result.repair_iterations, result.remaining_violations,
                result.ls_iterations, result.total_ls_improvement,
                result.successful_moves_simple, result.successful_moves_interchange, result.successful_moves_deactivate,
                result.warmstart_used, result.warmstart_time, result.warmstart_obj,
                "failed_y_init", "Y initialization failed")
            return result
        end
        
        # Phase 2: Solve transportation problem
        phase2_start = time()
        phase2_log_file = joinpath(output_dirs["logs"], "phase2_$(instance_name)_$(method).log")
        
        model_transport = solve_phase2_model(instance, Y)
        set_optimizer(model_transport, Gurobi.Optimizer)
        set_optimizer_attribute(model_transport, "LogFile", phase2_log_file)
        set_time_limit_sec(model_transport, config["time_limits"]["phase2"])
        
        optimize!(model_transport)
        
        if termination_status(model_transport) != MOI.OPTIMAL && termination_status(model_transport) != MOI.TIME_LIMIT
            result = ExperimentResult(result.instance_id, result.method, result.B, result.S, result.P,
                init_time, time() - phase2_start, result.repair_time, result.ls_time, time() - total_start_time,
                result.phase2_obj, result.split_resolved_obj, result.post_repair_obj, result.post_ls_obj,
                result.num_split_bus, result.prop_split_bus,
                result.activity_violations_initial, result.risk_violations_initial,
                result.repair_iterations, result.remaining_violations,
                result.ls_iterations, result.total_ls_improvement,
                result.successful_moves_simple, result.successful_moves_interchange, result.successful_moves_deactivate,
                result.warmstart_used, result.warmstart_time, result.warmstart_obj,
                "failed_phase2", "Phase 2 optimization failed")
            return result
        end
        
        x_continuous = value.(model_transport[:x])
        phase2_obj = objective_value(model_transport)
        phase2_time = time() - phase2_start
        
        # Analyze splits
        num_splits, prop_splits, _ = analyze_splits(x_continuous, instance.S, instance.B)
        
        # Split resolution
        x_binary = split_resolution_heuristic(x_continuous, instance.S, instance.B)
        split_resolved_obj = dot(x_binary, instance.D)
        
        # Count initial violations
        initial_activity_viol, initial_risk_viol = count_violations(instance, x_binary, Y)
        
        # Repair constraints
        x_repaired, repair_time, repair_iters, initial_viols, remaining_viols = repair_with_tracking(instance, x_binary, Y)
        post_repair_obj = dot(x_repaired, instance.D)
        
        # Create solution for local search
        sol_before_ls = Solution(instance, x_repaired, Y, post_repair_obj, 1)
        
        # Calculate targets
        targets_lower, targets_upper = calculate_targets_optimized(instance)
        
        # Local search
        sol_after_ls, ls_time, ls_iters, total_improvement, succ_simple, succ_interchange, succ_deactivate = 
            local_search_with_tracking(sol_before_ls, targets_lower, targets_upper, config)
        
        post_ls_obj = sol_after_ls.Weight
        
        # Save solution
        sol_dir = joinpath(output_dirs["solutions"], size_str, method)
        if !isdir(sol_dir)
            mkpath(sol_dir)
        end
        sol_file = joinpath(sol_dir, "sol_$(instance_name)_$(method).jld2")
        write_solution(sol_after_ls, sol_file)
        
        # Plot solution
        plot_dir = joinpath(output_dirs["plots"], size_str, method)
        if !isdir(plot_dir)
            mkpath(plot_dir)
        end
        plot_file = joinpath(plot_dir, "plot_$(instance_name)_$(method).png")
        plot_solution(sol_after_ls, plot_file)
        
        # Optional warmstart
        warmstart_used = false
        warmstart_time = 0.0
        warmstart_obj = post_ls_obj
        
        if get(config["warmstart"], "enable", false)
            warmstart_start = time()
            warmstart_log_file = joinpath(output_dirs["logs"], "warmstart_$(instance_name)_$(method).log")
            
            model_warmstart = build_model(instance)
            set_start_value.(model_warmstart[:x], sol_after_ls.X)
            set_start_value.(model_warmstart[:y], sol_after_ls.Y)
            set_optimizer(model_warmstart, Gurobi.Optimizer)
            set_optimizer_attribute(model_warmstart, "LogFile", warmstart_log_file)
            set_optimizer_attribute(model_warmstart, "Method", config["warmstart"]["method"])
            set_time_limit_sec(model_warmstart, config["time_limits"]["warmstart"])
            
            optimize!(model_warmstart)
            
            if termination_status(model_warmstart) == MOI.OPTIMAL || termination_status(model_warmstart) == MOI.TIME_LIMIT
                warmstart_obj = objective_value(model_warmstart)
                warmstart_used = true
            end
            
            warmstart_time = time() - warmstart_start
        end
        
        total_time = time() - total_start_time
        
        # Create result
        result = ExperimentResult(
            instance_name, method, B, S, P,
            init_time, phase2_time, repair_time, ls_time, total_time,
            phase2_obj, split_resolved_obj, post_repair_obj, post_ls_obj,
            num_splits, prop_splits,
            initial_activity_viol, initial_risk_viol,
            repair_iters, remaining_viols,
            ls_iters, total_improvement, succ_simple, succ_interchange, succ_deactivate,
            warmstart_used, warmstart_time, warmstart_obj,
            "success", ""
        )
        
        # Log result
        log_message = """
        Completed $instance_name with $method:
          Phase 2 obj: $phase2_obj
          Split BUs: $(round(prop_splits * 100, digits=2))%
          Split resolved obj: $split_resolved_obj  
          Post repair obj: $post_repair_obj
          Post LS obj: $post_ls_obj
          Total time: $(round(total_time, digits=2))s
        """
        println(log_message)
        open(experiment_log, "a") do f
            write(f, log_message * "\n")
        end
        
    catch e
        error_msg = sprint(showerror, e)
        println("ERROR processing $instance_name with $method: $error_msg")
        
        result = ExperimentResult(
            instance_name, method, B, S, P,
            result.init_time, result.phase2_time, result.repair_time, result.ls_time, time() - total_start_time,
            result.phase2_obj, result.split_resolved_obj, result.post_repair_obj, result.post_ls_obj,
            result.num_split_bus, result.prop_split_bus,
            result.activity_violations_initial, result.risk_violations_initial,
            result.repair_iterations, result.remaining_violations,
            result.ls_iterations, result.total_ls_improvement,
            result.successful_moves_simple, result.successful_moves_interchange, result.successful_moves_deactivate,
            result.warmstart_used, result.warmstart_time, result.warmstart_obj,
            "error", error_msg
        )
    end
    
    return result
end

# Main experiment runner
function run_heuristic_experiments(config_file="config_heuristic.toml")
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
        Heuristic Comparison Experiment
        ID: $experiment_id
        Started: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        Configuration: $config_file
        
        """)
    end
    
    # Get list of instance files recursively from all subdirectories
    input_dir = config["general"]["input_directory"]
    instance_files = String[]
    
    # Recursively find all .jld2 files
    for (root, dirs, files) in walkdir(input_dir)
        for file in files
            if endswith(file, ".jld2")
                push!(instance_files, joinpath(root, file))
            end
        end
    end
    
    if isempty(instance_files)
        error("No instance files found in $input_dir or its subdirectories")
    end
    
    # Sort files for consistent processing order
    sort!(instance_files)
    
    println("Found $(length(instance_files)) instance files across all subdirectories")
    
    # Group by size_str for reporting
    instance_groups = Dict{String, Vector{String}}()
    for file in instance_files
        # Try to extract B_S_P from path or filename
        dir_name = basename(dirname(file))
        if occursin(r"\d+_\d+_\d+", dir_name)
            size_str = dir_name
        else
            # Try to extract from filename
            match_result = match(r"(\d+)_(\d+)_(\d+)", basename(file))
            if match_result !== nothing
                size_str = match_result.match
            else
                size_str = "unknown"
            end
        end
        
        if !haskey(instance_groups, size_str)
            instance_groups[size_str] = String[]
        end
        push!(instance_groups[size_str], file)
    end
    
    println("\nInstance distribution:")
    for (size_str, files) in sort(collect(instance_groups))
        println("  $size_str: $(length(files)) instances")
    end
    
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
    
    println("Testing methods: ", join(methods_to_test, ", "))
    
    # Results storage
    all_results = ExperimentResult[]
    
    # Process each instance
    for (idx, instance_file) in enumerate(instance_files)
        instance_name = basename(instance_file)
        println("\n=== Processing instance $idx/$(length(instance_files)): $instance_name ===")
        
        # Load instance
        instance = read_instance(instance_file)
        
        # Maintain directory structure when copying instance file
        size_str = "$(instance.B)_$(instance.S)_$(instance.P)"
        inst_dir = joinpath(output_dirs["instances"], size_str)
        mkpath(inst_dir)
        cp(instance_file, joinpath(inst_dir, instance_name))
        
        # Test each method
        for method in methods_to_test
            result = process_instance_method(instance, instance_file, method, config, output_dirs, experiment_log)
            push!(all_results, result)
            
            # Save intermediate results after each method to avoid data loss
            save_results(all_results, output_dirs["results"])
        end
    end
    
    # Generate final summary
    generate_summary(all_results, output_dirs["results"], experiment_log)
    
    # Log completion
    completion_msg = "\n=== Experiment Complete at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) ==="
    println(completion_msg)
    open(experiment_log, "a") do f
        write(f, completion_msg * "\n")
    end
    
    println("\nResults saved to: $experiment_dir")
end

# Save results to CSV and JSON
function save_results(results::Vector{ExperimentResult}, results_dir::String)
    # Convert to DataFrame for CSV
    df = DataFrame(
        instance_id = [r.instance_id for r in results],
        method = [r.method for r in results],
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
    
    # Save as JSON for detailed analysis
    json_results = [Dict(
        "instance_id" => r.instance_id,
        "method" => r.method,
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
            "post_ls" => r.post_ls_obj,
            "warmstart" => r.warmstart_obj
        ),
        "splits" => Dict(
            "count" => r.num_split_bus,
            "proportion" => r.prop_split_bus
        ),
        "violations" => Dict(
            "initial_activity" => r.activity_violations_initial,
            "initial_risk" => r.risk_violations_initial,
            "repair_iterations" => r.repair_iterations,
            "remaining" => r.remaining_violations
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
        "status" => r.status,
        "error" => r.error_message
    ) for r in results]
    
    open(joinpath(results_dir, "detailed_results.json"), "w") do f
        JSON.print(f, json_results, 2)
    end
end

# Generate summary statistics
function generate_summary(results::Vector{ExperimentResult}, results_dir::String, log_file::String)
    # Group by method
    methods = unique([r.method for r in results])
    
    summary_df = DataFrame()
    
    for method in methods
        method_results = filter(r -> r.method == method && r.status == "success", results)
        
        if !isempty(method_results)
            summary_row = DataFrame(
                method = method,
                num_instances = length(method_results),
                avg_total_time = mean([r.total_time for r in method_results]),
                avg_post_repair_obj = mean([r.post_repair_obj for r in method_results]),
                avg_post_ls_obj = mean([r.post_ls_obj for r in method_results]),
                avg_ls_improvement = mean([r.total_ls_improvement for r in method_results]),
                avg_ls_iterations = mean([r.ls_iterations for r in method_results]),
                avg_prop_split_bus = mean([r.prop_split_bus for r in method_results]),
                success_rate = length(method_results) / length(filter(r -> r.method == method, results))
            )
            summary_df = vcat(summary_df, summary_row)
        end
    end
    
    CSV.write(joinpath(results_dir, "summary_results.csv"), summary_df)
    
    # Also write summary to log
    summary_text = "\n\n=== SUMMARY STATISTICS ===\n"
    for row in eachrow(summary_df)
        summary_text *= "\nMethod: $(row.method)\n"
        summary_text *= "  Success rate: $(round(row.success_rate * 100, digits=1))%\n"
        summary_text *= "  Avg time: $(round(row.avg_total_time, digits=2))s\n"
        summary_text *= "  Avg final obj: $(round(row.avg_post_ls_obj, digits=2))\n"
        summary_text *= "  Avg LS improvement: $(round(row.avg_ls_improvement, digits=2))\n"
    end
    
    println(summary_text)
    open(log_file, "a") do f
        write(f, summary_text)
    end
end

# Main function
function main()
    if length(ARGS) == 0
        println("Usage: julia run_heuristic_experiments.jl config_heuristic.toml")
        return
    end
    
    config_file = ARGS[1]
    run_heuristic_experiments(config_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end