using Gurobi, JuMP, JLD2
using Types
using MathOptInterface
const MOI = MathOptInterface
using DelimitedFiles
using Plots
using CPLEX
include("../heuristics/constructive.jl")
include("first_model.jl")

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
        println("\nAnalyzing split severity:")
        for j in split_branches
            # Get all non-zero assignments for this branch
            assignments = [(i, x_continuous[i, j]) for i in 1:S if x_continuous[i, j] > tolerance]
            # Sort by assignment value (descending)
            sort!(assignments, by = x -> x[2], rev = true)
            
            println("  BU $j assignments:")
            for (i, val) in assignments
                println("    Center $i: $(round(val, digits=4))")
            end
        end
    end
    
    return split_count, split_proportion, split_branches
end

function repair_activity_constraints(instance::Instance, x_binary, y_fixed)
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
    
    # Check which services violate their activity constraints
    violations = Dict{Tuple{Int, Int}, Float64}() # (service, activity) => violation amount
    
    for i in 1:S
        if y_fixed[i] == 1 # Only check active services
            for M in 1:m
                lower_bound = (1-ε) * μ[M][i]
                upper_bound = (1+ε) * μ[M][i]
                
                if activity_levels[i, M] < ceil(lower_bound)
                    violations[(i, M)] = ceil(lower_bound) - activity_levels[i, M]
                elseif activity_levels[i, M] > floor(upper_bound)
                    violations[(i, M)] = activity_levels[i, M] - floor(upper_bound)
                end
            end
        end
    end
    
    if isempty(violations)
        println("No activity violations to repair!")
        return x_repaired
    end
    
    println("Found $(length(violations)) activity violations")
    
    # Sort violations by magnitude (largest first)
    sorted_violations = sort(collect(violations), by=x->x[2], rev=true)
    
    # Maximum number of repair iterations to avoid infinite loops
    max_repair_iterations = 100
    iteration = 0
    
    while !isempty(violations) && iteration < max_repair_iterations
        iteration += 1
        
        # Take the largest violation
        (service_i, activity_M), violation_amount = first(sorted_violations)
        
        # Determine if we need to increase or decrease activity
        needs_increase = activity_levels[service_i, activity_M] < ceil((1-ε) * μ[activity_M][service_i])
        
        if needs_increase
            # Need to increase activity - look for branches to add
            # Find branches currently assigned to other services
            other_service_branches = [(j, i) for j in 1:B, i in 1:S if x_repaired[i, j] == 1 && i != service_i]
            
            # Calculate the benefit of reassigning each branch
            candidates = []
            for (branch_j, current_service) in other_service_branches
                # Skip if this would violate risk constraint
                if sum(x_repaired[service_i, j] * R[j] for j in 1:B) + R[branch_j] > β[service_i]
                    continue
                end
                
                # FIX: Check if adding this branch would cause violations for any activity measure
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
                
                # Calculate the activity gain for the deficient service
                activity_gain = V[activity_M][branch_j]
                
                # Calculate the activity loss for the current service
                activity_loss = V[activity_M][branch_j]
                
                # Check if this reassignment would cause a violation for the current service
                current_level = activity_levels[current_service, activity_M]
                new_level = current_level - activity_loss
                
                # Skip if it would cause a new violation
                if new_level < ceil((1-ε) * μ[activity_M][current_service])
                    continue
                end
                
                # Calculate net benefit (weighted by how much each service is out of bounds)
                current_violation = get(violations, (current_service, activity_M), 0.0)
                benefit = activity_gain - (current_violation > 0 ? activity_loss : 0)
                
                # Calculate cost impact
                cost_change = D[service_i, branch_j] - D[current_service, branch_j]
                
                push!(candidates, (branch_j, current_service, benefit, cost_change))
            end
            
            # Sort candidates by benefit (descending) and then by cost change (ascending)
            sort!(candidates, by=x->(x[3], -x[4]), rev=true)
            
            # Try the best candidate
            if !isempty(candidates)
                branch_j, current_service, _, _ = first(candidates)
                
                # Reassign the branch
                x_repaired[service_i, branch_j] = 1
                x_repaired[current_service, branch_j] = 0
                
                # Update activity levels for ALL activity types
                for M in 1:m
                    activity_levels[service_i, M] += V[M][branch_j]
                    activity_levels[current_service, M] -= V[M][branch_j]
                end
                
                println("Iteration $iteration: Reassigned BU $branch_j from center $current_service to center $service_i to increase activity $activity_M")
            else
                println("Iteration $iteration: Could not find a suitable BU to increase activity $activity_M for center $service_i")
                # Mark this violation as unresolvable for now to avoid getting stuck
                delete!(violations, (service_i, activity_M))
                sorted_violations = sort(collect(violations), by=x->x[2], rev=true)
                continue
            end
        else
            # Need to decrease activity - look for branches to remove
            # Find branches currently assigned to this service
            current_branches = [j for j in 1:B if x_repaired[service_i, j] == 1]
            
            # Calculate the benefit of reassigning each branch
            candidates = []
            for branch_j in current_branches
                # Calculate the activity reduction
                activity_reduction = V[activity_M][branch_j]
                
                # Find potential services to receive this branch
                receiving_services = [i for i in 1:S if i != service_i && y_fixed[i] == 1]
                
                for receiving_service in receiving_services
                    # Skip if this would violate risk constraint
                    if sum(x_repaired[receiving_service, j] * R[j] for j in 1:B) + R[branch_j] > β[receiving_service]
                        continue
                    end
                    
                    # FIX: Check if adding this branch would cause violations for any activity measure
                    valid_assignment = true
                    for check_M in 1:m
                        new_level = activity_levels[receiving_service, check_M] + V[check_M][branch_j]
                        if new_level > floor((1+ε) * μ[check_M][receiving_service])
                            valid_assignment = false
                            break
                        end
                    end
                    
                    if !valid_assignment
                        continue
                    end
                    
                    # Check if this reassignment would cause a violation for the receiving service
                    receiving_level = activity_levels[receiving_service, activity_M]
                    new_level = receiving_level + activity_reduction
                    
                    # Calculate net benefit
                    receiving_violation = get(violations, (receiving_service, activity_M), 0.0)
                    benefit = activity_reduction - (receiving_violation > 0 ? activity_reduction : 0)
                    
                    # Calculate cost impact
                    cost_change = D[receiving_service, branch_j] - D[service_i, branch_j]
                    
                    push!(candidates, (branch_j, receiving_service, benefit, cost_change))
                end
            end
            
            # Sort candidates by benefit (descending) and then by cost change (ascending)
            sort!(candidates, by=x->(x[3], -x[4]), rev=true)
            
            # Try the best candidate
            if !isempty(candidates)
                branch_j, receiving_service, _, _ = first(candidates)
                
                # Reassign the branch
                x_repaired[service_i, branch_j] = 0
                x_repaired[receiving_service, branch_j] = 1
                
                # Update activity levels for ALL activity types
                for M in 1:m
                    activity_levels[service_i, M] -= V[M][branch_j]
                    activity_levels[receiving_service, M] += V[M][branch_j]
                end
                
                println("Iteration $iteration: Reassigned BU $branch_j from center $service_i to center $receiving_service to decrease activity $activity_M")
            else
                println("Iteration $iteration: Could not find a suitable center to receive BU from center $service_i to decrease activity $activity_M")
                # Mark this violation as unresolvable for now to avoid getting stuck
                delete!(violations, (service_i, activity_M))
                sorted_violations = sort(collect(violations), by=x->x[2], rev=true)
                continue
            end
        end
        
        # Recalculate violations after the reassignment
        violations = Dict{Tuple{Int, Int}, Float64}()
        for i in 1:S
            if y_fixed[i] == 1
                for M in 1:m
                    lower_bound = (1-ε) * μ[M][i]
                    upper_bound = (1+ε) * μ[M][i]
                    
                    if activity_levels[i, M] < ceil(lower_bound)
                        violations[(i, M)] = ceil(lower_bound) - activity_levels[i, M]
                    elseif activity_levels[i, M] > floor(upper_bound)
                        violations[(i, M)] = activity_levels[i, M] - floor(upper_bound)
                    end
                end
            end
        end
        
        sorted_violations = sort(collect(violations), by=x->x[2], rev=true)
    end
    
    # Final check to see if we resolved all violations
    if isempty(violations)
        println("Successfully repaired all activity violations in $iteration iterations!")
    else
        println("Could not repair all violations after $iteration iterations. Remaining violations: $(length(violations))")
    end
    
    # Verify all constraints are satisfied
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

    set_time_limit_sec(model, 600)
    
    return model, d
end

function main()
    instance = read_instance(ARGS[1])
    #Y, time = pdisp_2(instance)
    #println(Y)
    #Y_bool = zeros(Int, instance.S)
    #for idx in Y
    #    Y_bool[idx] = 1
    #end
    #plot_ys2(instance, Y_bool, "pdisp_plotted_all.png")
    #println(Y_bool)
    #model_original = build_model_x_relaxed(instance)
    #set_optimizer(model_original, Gurobi.Optimizer)
    #set_time_limit_sec(model_original, 600)
    #write_to_file(model_original, "original_model_grb.lp")
    #optimize!(model_original)
    #x_solution = value.(model_original[:x])
    #y_sol = value.(model_original[:y])
    #X = round.(Int, x_solution)
    model, d = multi_pdp_min_instance(instance)
    optimize!(model)
    y_sol = value.(model[:y])
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
    x_repaired = repair_activity_constraints(instance, x_binary, Y2)
    #analyze_splits(x_solution)
    model_warmstart = build_model(instance)
    x = model_warmstart[:x]
    y = model_warmstart[:y]
    set_start_value.(x, x_repaired)
    set_start_value.(y, Y2)
    set_optimizer(model_warmstart, Gurobi.Optimizer)
    set_time_limit_sec(model_warmstart, 1800)
    optimize!(model_warmstart)
end
main()