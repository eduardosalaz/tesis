# Multi-Type P-Median Problem
# Select P centers to minimize total distance from branches to nearest centers
# While respecting type constraints

using JuMP
using Gurobi
using LinearAlgebra
using Printf

# Helper functions for Benders preprocessing
function calculate_dist_min(D::Matrix)
    return [minimum(filter(x -> x != typemax(eltype(D)), D[:,j])) for j in axes(D,2)]
end

function calculate_p_max(D::Matrix, p_value::Int)
    n = size(D, 2)
    if p_value > n
        return fill(nothing, n)
    end
    
    return [partialsort(filter(x -> x != typemax(eltype(D)), D[:,j]), p_value, rev=true) for j in axes(D,2)]
end

"""
    solve_multi_type_p_median(instance; kwargs...)

Solve multi-type p-median problem to select P centers that minimize
total distance from all branches to their nearest selected center.

This is the RIGHT model for facility location serving customers!

Arguments:
- instance: Your standard instance with:
  - S: number of candidate centers
  - P: number of centers to select
  - B: number of branches to serve
  - S_coords: Sx2 matrix of center coordinates
  - B_coords: Bx2 matrix of branch coordinates (if available)
  - D: SxB distance matrix (center i to branch j)
  - Sk: vector of vectors, where Sk[k] contains indices of type k centers
  - Lk: minimum number of type k centers
  - Uk: maximum number of type k centers

Returns:
- Y: binary vector of selected centers
- info: solution information
"""
function solve_multi_type_p_median(instance;
                                  time_limit=600,
                                  gap_tolerance=0.0001,
                                  verbose=true,
                                  formulation=:strong)  # :strong or :benders
    
    S = instance.S  # number of centers
    P = instance.P  # number to select
    B = instance.B  # number of branches
    n_types = length(instance.Sk)
    
    if verbose
        println("\n" * "="^70)
        println("MULTI-TYPE P-MEDIAN PROBLEM")
        println("="^70)
        println("Candidate centers: $S")
        println("Centers to select: $P")
        println("Branches to serve: $B")
        println("Types: $n_types")
        println("Formulation: $formulation")
    end
    
    # Get distance matrix
    if hasfield(typeof(instance), :D) && size(instance.D) == (S, B)
        D = instance.D
        if verbose
            println("Using provided distance matrix D")
        end
    else
        # Compute from coordinates
        if verbose
            println("Computing distance matrix from coordinates...")
        end
        D = compute_distance_matrix_centers_branches(instance)
    end
    
    if formulation == :strong
        return solve_strong_formulation(instance, D, S, P, B, n_types, time_limit, gap_tolerance, verbose)
    elseif formulation == :benders
        return solve_benders_formulation(instance, D, S, P, B, n_types, time_limit, gap_tolerance, verbose)
    else
        error("Unknown formulation: $formulation. Use :strong or :benders")
    end
end

function solve_strong_formulation(instance, D, S, P, B, n_types, time_limit, gap_tolerance, verbose)
    # Create model
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", time_limit)
    set_optimizer_attribute(model, "MIPGap", gap_tolerance)
    if !verbose
        set_silent(model)
    end
    
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
    @constraint(
        model,
        low_k[K in 1:k],
        instance.Lk[K] <= sum(y[i] for i in instance.Sk[K]),
    )

    @constraint(
        model,
        upp_k[K in 1:k],
        sum(y[i] for i in instance.Sk[K]) <= instance.Uk[K],
    )
    
    # Solve
    start_time = time()
    optimize!(model)
    solve_time = time() - start_time
    
    return extract_solution(model, y, S, P, B, n_types, instance, solve_time, verbose, D, :strong)
end

function solve_benders_formulation(instance, D, S, P, B, n_types, time_limit, gap_tolerance, verbose)
    # Preprocess distance bounds
    dist_min = calculate_dist_min(D)
    p_max = calculate_p_max(D, P)
    
    # Master problem
    master = Model(Gurobi.Optimizer)
    set_optimizer_attribute(master, "TimeLimit", time_limit)
    set_optimizer_attribute(master, "MIPGap", gap_tolerance)
    if !verbose
        set_silent(master)
    end
    
    @variable(master, y[1:S], Bin)
    @variable(master, z[j in 1:B], lower_bound=dist_min[j], upper_bound=p_max[j])
    
    # Objective
    @objective(master, Min, sum(z))
    
    # Select exactly P centers
    @constraint(master, sum(y) == P)
    
    # Type constraints
    k = instance.K
    @constraint(
        master,
        low_k[K in 1:k],
        instance.Lk[K] <= sum(y[i] for i in instance.Sk[K]),
    )

    @constraint(
        master,
        upp_k[K in 1:k],
        sum(y[i] for i in instance.Sk[K]) <= instance.Uk[K],
    )
    
    # Set up for lazy constraints
    set_optimizer_attribute(master, "LazyConstraints", 1)
    #set_optimizer_attribute(master, "BranchDir", 1)
    #set_optimizer_attribute(master, "MIPGapAbs", 0.9999)
    set_optimizer_attribute(master, "StartNodeLimit", 2000000)
    
    # Benders callback
    function benders_callback(cb_data, cb_where::Cint)
        if cb_where == GRB_CB_MIPSOL
            Gurobi.load_callback_variable_primal(cb_data, cb_where)
            ȳ = callback_value.(cb_data, master[:y])
            
            sub = Model(Gurobi.Optimizer)
            set_optimizer_attribute(sub, "OutputFlag", 0)
            set_optimizer_attribute(sub, "InfUnbdInfo", 0)
            
            # Add bounds to dual variables
            @variable(sub, u >= 0)  # Add lower bound to u
            @variable(sub, v[1:S] >= 0)  # Add lower bound to v
            
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
                else
                    if verbose
                        println("Warning: Subproblem for j=$j is not optimal. Status: $(termination_status(sub))")
                    end
                end
            end
            
            for cut in cuts
                MOI.submit(master, MOI.LazyConstraint(cb_data), cut)
            end
        end
    end
    
    MOI.set(master, Gurobi.CallbackFunction(), benders_callback)
    
    # Solve
    start_time = time()
    optimize!(master)
    solve_time = time() - start_time
    
    return extract_solution(master, y, S, P, B, n_types, instance, solve_time, verbose, D, :benders)
end

function extract_solution(model, y, S, P, B, n_types, instance, solve_time, verbose, D, formulation)
    # Extract solution
    status = termination_status(model)
    
    if status in [MOI.OPTIMAL, MOI.TIME_LIMIT]
        Y_solution = round.(Int, value.(y))
        objective_val = objective_value(model)
        gap = relative_gap(model)
        
        if verbose
            println("\n" * "-"^70)
            println("SOLUTION FOUND")
            println("-"^70)
            println("Status: $status")
            println("Objective (total distance): $(round(objective_val, digits=2))")
            println("Average distance per branch: $(round(objective_val/B, digits=2))")
            println("Gap: $(round(gap*100, digits=2))%")
            println("Time: $(round(solve_time, digits=2))s")
            
            # Type distribution
            println("\nSelected centers by type:")
            for k in 1:n_types
                count = sum(Y_solution[i] for i in instance.Sk[k])
                percent = count / P * 100
                @printf("  Type %d: %2d centers (%5.1f%%) [bounds: %d-%d]\n",
                        k, count, percent, instance.Lk[k], instance.Uk[k])
            end
            
            # Coverage analysis for strong formulation
            if formulation == :strong && has_values(model)
                x_val = value.(model[:x])
                
                # Average branches per center
                branches_per_center = [sum(x_val[i,:]) for i in 1:S if Y_solution[i] > 0]
                println("\nWorkload distribution:")
                println("  Branches per center: min=$(minimum(branches_per_center)), " *
                       "avg=$(round(mean(branches_per_center), digits=1)), " *
                       "max=$(maximum(branches_per_center))")
                
                # Distance distribution
                distances = [D[i,j] for i in 1:S, j in 1:B if x_val[i,j] > 0.5]
                println("\nDistance statistics:")
                println("  Min: $(round(minimum(distances), digits=2))")
                println("  Avg: $(round(mean(distances), digits=2))")
                println("  Max: $(round(maximum(distances), digits=2))")
                println("  90th percentile: $(round(quantile(distances, 0.9), digits=2))")
            end
        end
        
        info = Dict(
            "objective" => objective_val,
            "average_distance" => objective_val / B,
            "gap" => gap,
            "time" => solve_time,
            "status" => status,
            "formulation" => formulation
        )
        
        return Y_solution, info
        
    else
        error("Failed to solve: $status")
    end
end

# Helper: Compute distances from centers to branches
function compute_distance_matrix_centers_branches(instance)
    S = instance.S
    B = instance.B
    
    # Check what coordinates we have
    if hasfield(typeof(instance), :B_coords)
        # We have both center and branch coordinates
        D = zeros(S, B)
        for i in 1:S
            for j in 1:B
                D[i,j] = norm(instance.S_coords[i,:] - instance.B_coords[j,:])
            end
        end
        return D
    else
        # Generate branch locations or use heuristic
        error("Need either D matrix or B_coords in instance")
    end
end

# Comparison with p-median (no types) and p-dispersion
function compare_objectives(instance; time_limit=120)
    println("="^80)
    println("COMPARING DIFFERENT OBJECTIVES")
    println("="^80)
    
    # 1. Multi-type P-median (serves branches well)
    println("\n1. MULTI-TYPE P-MEDIAN (minimize total branch distance)")
    Y_median, info_median = solve_multi_type_p_median(instance, 
                                                     time_limit=time_limit,
                                                     verbose=false)
    println("   Total distance: $(round(info_median["objective"], digits=2))")
    println("   Avg per branch: $(round(info_median["average_distance"], digits=2))")
    
    # 2. P-sum-dispersion (spreads centers)
    println("\n2. P-SUM-DISPERSION (maximize center separation)")
    include("p_sum_dispersion_plugin.jl")
    Y_dispersion, info_disp = solve_p_sum_dispersion(instance,
                                                    time_limit=time_limit,
                                                    verbose=false)
    
    # Evaluate p-sum solution on p-median objective
    total_dist_disp = evaluate_median_objective(Y_dispersion, instance)
    println("   Total distance: $(round(total_dist_disp, digits=2))")
    println("   Avg per branch: $(round(total_dist_disp/instance.B, digits=2))")
    
    # Compare
    println("\n3. COMPARISON")
    improvement = (total_dist_disp - info_median["objective"]) / total_dist_disp * 100
    println("   P-median is $(round(improvement, digits=1))% better for serving branches")
    
    return Y_median, Y_dispersion
end

# Evaluate median objective for any solution
function evaluate_median_objective(Y, instance)
    total_distance = 0.0
    
    for j in 1:instance.B
        # Find nearest selected center
        min_dist = Inf
        for i in 1:instance.S
            if Y[i] > 0 && instance.D[i,j] < min_dist
                min_dist = instance.D[i,j]
            end
        end
        total_distance += min_dist
    end
    
    return total_distance
end

# Heuristic for large instances
function p_median_greedy_heuristic(instance; verbose=false)
    S = instance.S
    P = instance.P
    B = instance.B
    D = instance.D
    
    # Start with empty
    selected = Int[]
    remaining = collect(1:S)
    
    # Greedy: repeatedly add center that reduces total distance most
    for iter in 1:P
        best_reduction = Inf
        best_center = -1
        
        # Try each remaining center
        for i in remaining
            # Calculate total distance if we add center i
            total_dist = 0.0
            for j in 1:B
                # Distance to nearest among selected ∪ {i}
                min_dist = D[i,j]
                for k in selected
                    min_dist = min(min_dist, D[k,j])
                end
                total_dist += min_dist
            end
            
            if total_dist < best_reduction
                best_reduction = total_dist
                best_center = i
            end
        end
        
        # Add best center
        push!(selected, best_center)
        filter!(x -> x != best_center, remaining)
        
        if verbose && iter % 5 == 0
            println("Added center $iter: total distance = $(round(best_reduction, digits=2))")
        end
    end
    
    # Convert to binary
    Y = zeros(Int, S)
    Y[selected] .= 1
    
    return Y
end

# Test function
function test_p_median()
    # Create test instance with branches
    Random.seed!(123)
    
    S = 100  # centers
    P = 10   # select
    B = 200  # branches
    
    # Center locations (in 4 regions for types)
    S_coords = zeros(S, 2)
    for i in 1:25
        S_coords[i, :] = randn(2) * 5 .+ [-20, -20]
    end
    for i in 26:50
        S_coords[i, :] = randn(2) * 5 .+ [20, -20]
    end
    for i in 51:75
        S_coords[i, :] = randn(2) * 5 .+ [-20, 20]
    end
    for i in 76:100
        S_coords[i, :] = randn(2) * 5 .+ [20, 20]
    end
    
    # Branch locations (more concentrated in populated areas)
    B_coords = vcat(
        randn(80, 2) * 8 .+ [0, 0],      # 40% central
        randn(40, 2) * 5 .+ [-20, -20],  # 20% region 1
        randn(40, 2) * 5 .+ [20, -20],   # 20% region 2
        randn(20, 2) * 5 .+ [-20, 20],   # 10% region 3
        randn(20, 2) * 5 .+ [20, 20]     # 10% region 4
    )
    
    # Compute distance matrix
    D = zeros(S, B)
    for i in 1:S
        for j in 1:B
            D[i,j] = norm(S_coords[i,:] - B_coords[j,:])
        end
    end
    
    instance = (
        S = S,
        P = P,
        B = B,
        S_coords = S_coords,
        B_coords = B_coords,
        D = D,
        Sk = [1:25, 26:50, 51:75, 76:100],
        Lk = [1, 1, 1, 1],  # at least 1 of each
        Uk = [4, 4, 4, 4]   # at most 4 of each
    )
    
    # Solve
    Y, info = solve_multi_type_p_median(instance, time_limit=60)
    
    return Y, info, instance
end