using Gurobi, JuMP, JLD2
using Types
using MathOptInterface
const MOI = MathOptInterface
using DelimitedFiles
using Plots
using CPLEX
include("../heuristics/constructive.jl")
include("first_model.jl")
include("../heuristics/ls.jl")

function pdisp_heuristic(instance)
    Y_integer, time = pdisp_2(instance)
    Y_bool = zeros(Int, instance.S)
    Y_bool[Y_integer] .= 1
    return Y_bool
end

# Phase 2: Fix Y values and solve transportation problem
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
    #@constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y_fixed[i])
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
    k = 4 #  of branches

    model = Model() # THIS IS WHERE THE FUN BEGINS

    #@variable(model, x[1:S, 1:B], Bin)
    @variable(model, x[1:S, 1:B], lower_bound = 0, upper_bound = 1)
    # num suc and num bu, Xᵢⱼ

    @variable(model, y[1:S], Bin)
    # Yᵢ

    @objective(model, Min, sum(D .* x))
    # Xᵢⱼ * Dᵢⱼ

    @constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)

    # ∑ᵢ∈S Xᵢⱼ = 1, ∀ j ∈ B

    @constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y[i])

    # Xᵢⱼ ≤ Yᵢ , ∀ i ∈ S, j ∈ B

    @constraint(model, cardinality, sum(y) == P)

    # ∑ i ∈ S Yᵢ = p

    @constraint(model, risk[i in 1:S], sum(x[i, j] * R[j] for j in 1:B) <= β[i])

    # ∑ j ∈ B Xᵢⱼ Rⱼ ≤ βᵢ, ∀ i ∈ S

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

    # Yᵢμₘⁱ(1-tᵐ) ≤ ∑i∈S Xᵢⱼ vⱼᵐ ≤ Yᵢμₘʲ(1+tᵐ) ∀ j ∈ B, m = 1 … 3

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

    # lₖ ≤ ∑i ∈ Sₖ Yᵢ ≤ uₖ, k = 1 … 4
    return model
end

function analyze_splits(x_sol::Matrix{Float64}, epsilon::Float64=1e-6)
    splits = 0
    S, B = size(x_sol)
    split_info = Dict()
    
    for j in 1:B
        assignments = [(i, x_sol[i,j]) for i in 1:S if abs(x_sol[i,j]) > epsilon]
        if length(assignments) > 1
            splits += 1
            println("\nBasic unit $j is split among $(length(assignments)) centers:")
            for (i, val) in sort(assignments, by=x->x[2], rev=true)
                println("  center $i: $(round(val, digits=4))")
            end
            split_info[j] = assignments
        end
    end
    
    println("\nTotal number of split BUs: $splits")
    return splits, split_info
end

function split_resolution_heuristic(x_continuous, S, B)
    # Initialize binary solution matrix
    x_binary = zeros(Int, S, B)
    
    # For each branch j, find the service i with the largest x[i,j] value
    for j in 1:B
        # Find index of maximum value in column j
        max_val = 0.0
        max_idx = 0
        
        for i in 1:S
            if x_continuous[i, j] > max_val
                max_val = x_continuous[i, j]
                max_idx = i
            end
        end
        
        # Assign branch j to service with max value
        x_binary[max_idx, j] = 1
    end
    
    return x_binary
end



# Optional: Function to check if the binary solution is feasible
function analyze_splits(x_continuous, S, B)
    # Count how many branches have split assignments
    split_count = 0
    
    # For detailed analysis, track which branches are split
    split_branches = Int[]
    
    # Numerical tolerance for considering a value non-zero
    # (helps avoid counting very small values due to numerical precision)
    tolerance = 1e-6
    
    for j in 1:B
        # Count how many services have non-zero assignment for this branch
        non_zero_assignments = 0
        
        for i in 1:S
            if x_continuous[i, j] > tolerance
                non_zero_assignments += 1
            end
        end
        
        # If more than one service is assigned to this branch, it's a split
        if non_zero_assignments > 1
            split_count += 1
            push!(split_branches, j)
        end
    end
    
    # Calculate proportion of splits
    split_proportion = split_count / B
    
    println("Number of split BUs: ", split_count, " out of ", B, " total BUs")
    println("Proportion of split BUs: ", round(split_proportion * 100, digits=2), "%")
    println("Split BUs: ", split_branches)
    
    # Optionally, analyze the severity of splits
    if split_count > 0
        #println("\nAnalyzing split severity:")
        for j in split_branches
            # Get all non-zero assignments for this branch
            assignments = [(i, x_continuous[i, j]) for i in 1:S if x_continuous[i, j] > tolerance]
            # Sort by assignment value (descending)
            sort!(assignments, by = x -> x[2], rev = true)
            
            #println("  BU $j assignments:")
            #for (i, val) in assignments
            #    println("    Center $i: $(round(val, digits=4))")
            #end
        end
    end
    
    return split_count, split_proportion, split_branches
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

function multi_pdp_min_instance(instance)
    # Extract data from instance
    n = instance.S
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    p = instance.P
    
    # Calculate distance matrix
    metric = Distances.Euclidean()
    d = trunc.(Int, Distances.pairwise(metric, instance.S_coords, dims=1))
    
    # Maximum distance
    Dmax = maximum(d)
    
    # Number of types
    k = length(Sk)
    
    # Create model
    model = Model(Gurobi.Optimizer)
    
    # Variables
    @variable(model, y[1:n], Bin)  # 1 if point i is selected
    @variable(model, u >= 0)       # Minimum distance between any two selected points
    @variable(model, w[1:k], Int)  # Number of points of type k selected
    
    # Objective: Maximize the minimum distance
    @objective(model, Max, u)
    
    # Constraint: Linearization of minimum distance
    for i in 1:n
        for j in i+1:n
            @constraint(model, u <= d[i,j] + Dmax * (2 - (y[i] + y[j])))
        end
    end

    #@constraint(model, u >= minimum(d[i,j]) for i in 1:n for j in 1:n)
    
    # Constraint: Select exactly p points
    @constraint(model, sum(y) == p)
    
    # Constraint: Relationship between y and w variables
    for k_idx in 1:k
        @constraint(model, w[k_idx] == sum(y[i] for i in Sk[k_idx]))
    end
    
    # Constraints: Type balance
    for k_idx in 1:k
        @constraint(model, w[k_idx] >= Lk[k_idx])
        @constraint(model, w[k_idx] <= Uk[k_idx])
    end

    set_time_limit_sec(model, 300)
    
    return model, d
end

function main()
    instance = read_instance(ARGS[1])
    # phase 1
    #model, d = multi_pdp_min_instance(instance)
    #optimize!(model)
    model = build_model_x_relaxed(instance)
    set_optimizer(model, Gurobi.Optimizer)
    set_time_limit_sec(model, 180)
    optimize!(model)
    y_sol = value.(model[:y])
    #y_sol = pdisp_heuristic(instance)
    Y2 = round.(Int, y_sol)
    model_transport = solve_phase2_model(instance, Y2)
    set_optimizer(model_transport, Gurobi.Optimizer)
    optimize!(model_transport)
    x_continuous = value.(model_transport[:x])
    analyze_splits(x_continuous, instance.S, instance.B)
    x_binary = split_resolution_heuristic(x_continuous, instance.S, instance.B)
    weight = dot(x_binary, instance.D)
    sol_binary = Solution(instance, x_binary, Y2, weight, 1)
    println(isFactible(sol_binary, true))
    targets_lower, targets_upper = calculate_targets(instance)
    targets_lower_op, targets_upper_op = calculate_targets_optimized(instance)
    x_repaired = repair_activity_and_risk_constraints(instance, x_binary, Y2)
    weight = dot(x_repaired, instance.D)
    oldSol = Solution(instance, x_repaired, Y2, weight, 1)
    loop = 0
    any_improvement = true

    println(oldSol.Weight)
    start_time = time()
    while any_improvement
        loop += 1
        prev_weight = oldSol.Weight
        any_improvement = false  # Reset the flag at the start
        # First improvement function
        println("simple optimized bf")
        sol_moved_bu = simple_bu_improve_optimized(oldSol, targets_lower_op, targets_upper_op, :bf)
        new_weight_moved = sol_moved_bu.Weight
        println(new_weight_moved)

        if new_weight_moved < prev_weight
            any_improvement = true
            prev_weight = new_weight_moved
            oldSol = sol_moved_bu
        end

        # Second improvement function
        println("interchange optimized bf")
        sol_interchanged_bu = interchange_bu_improve_optimized(oldSol, targets_lower_op, targets_upper_op, :bf)
        new_weight_moved = sol_interchanged_bu.Weight
        println(new_weight_moved)

        if new_weight_moved < prev_weight
            any_improvement = true
            prev_weight = new_weight_moved
            oldSol = sol_interchanged_bu
        end
        println("---------------------")

        # Third improvement function
        println("deactivate ")
        sol_deactivated_center = deactivate_center_improve(oldSol, targets_lower, targets_upper)
        new_weight_moved = sol_deactivated_center.Weight
        println(new_weight_moved)

        if new_weight_moved < prev_weight
            any_improvement = true
            prev_weight = new_weight_moved
            oldSol = sol_deactivated_center
        end
        println("---------------------")

        @info any_improvement
    end
    end_time = time()
    println("LS time: ", end_time - start_time)
    println(oldSol.Weight)
    plot_solution(oldSol, "solution_laih_3000_600_40_1.png")
    #analyze_splits(x_solution)
    model_warmstart = build_model(instance)
    x = model_warmstart[:x]
    y = model_warmstart[:y]
    set_start_value.(x, oldSol.X)
    set_start_value.(y, oldSol.Y)
    set_optimizer(model_warmstart, Gurobi.Optimizer)
    #set_optimizer_attribute(model, "Method", 1)
    set_time_limit_sec(model_warmstart, 1800)
    optimize!(model_warmstart)
end
main()