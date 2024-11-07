using JuMP, Gurobi
using LinearAlgebra, Random
using DelimitedFiles


# Helper functions for solution validation and gap calculation
function calculate_gap(obj_lb::Float64, actual_wd::Float64)
    if abs(obj_lb) < 1e-10
        return Inf
    end
    return (actual_wd - obj_lb) / abs(obj_lb)
end

function get_closest_sites(distances::Matrix{Float64}, h_values::Vector{Int})
    n = size(distances, 1)
    H = Vector{Vector{Int}}(undef, n)
    for i in 1:n
        # Get indices sorted by distance
        sorted_indices = sortperm(distances[i, :])
        # Take h_i closest, including self
        H[i] = sorted_indices[1:h_values[i]]
    end
    return H
end



"""
    solve_pmedian_beamr(distances, p, demands; initial_h, delta, max_gap)

Solve the p-median problem using the BEAMR (Both Exact and Approximate Model Representation) approach.

# Arguments
- `distances::Matrix{Float64}`: Matrix of distances between nodes
- `p::Int`: Number of facilities to locate
- `demands::Vector{Float64}`: Vector of demand weights for each node
- `initial_h::Int=5`: Initial value for h_i parameters
- `delta::Int=5`: Increment size for h_i adjustments
- `max_gap::Float64=0.0`: Maximum acceptable optimality gap (0.0 for exact solution)

# Returns
Dictionary containing:
- "objective": Final objective value
- "facilities": Selected facility locations
- "h_values": Final h_i values
- "is_exact": Boolean indicating if solution is exact
- "gap": Optimality gap (0.0 for exact solutions)
"""
function solve_beamr_model(distances::Matrix{Float64}, p::Int,
    h_values::Vector{Int}, demands::Vector{Float64};
    time_limit::Float64=Inf,
    gap_tolerance::Float64=1e-3)
    n = size(distances, 1)

    # Get h_i+1 closest facility distances
    # Ensure h_values don't exceed n-1
    h_values_capped = min.(h_values, n - 1)

    # Get h_i+1 closest facility distances
    h_plus_one_dists = zeros(n)
    for i in 1:n
        sorted_dists = sort(distances[i, :])
        h_plus_one_dists[i] = sorted_dists[h_values_capped[i]+1]
    end

    # Get H_i sets
    H = get_closest_sites(distances, h_values)

    # Create model
    model = Model(Gurobi.Optimizer)

    # Set time limit and gap tolerance
    set_time_limit_sec(model, time_limit)
    set_optimizer_attribute(model, "MIPGap", gap_tolerance)
    set_silent(model)

    # Variables and constraints remain the same...
    @variable(model, x[i=1:n, j=1:n], Bin)
    @variable(model, f[1:n], Bin)

    for i in 1:n
        for j in 1:n
            if !(j in H[i])
                fix(x[i, j], 0)
            end
        end
    end

    @objective(model, Min,
        sum(demands[i] * distances[i, j] * x[i, j] for i in 1:n for j in H[i]) +
        sum(demands[i] * h_plus_one_dists[i] * f[i] for i in 1:n))

    for i in 1:n
        @constraint(model, sum(x[i, j] for j in H[i]) + f[i] == 1)
    end

    @constraint(model, sum(x[j, j] for j in 1:n) == p)

    for i in 1:n
        for j in H[i]
            if i != j
                @constraint(model, x[i, j] <= x[j, j])
            end
        end
    end

    optimize!(model)

    # Get solution status
    status = termination_status(model)
    println(status)
    println(JuMP.has_values(model))
    final = JuMP.has_values(model) && (status == OPTIMAL || status == TIME_LIMIT || status == SOLUTION_LIMIT)


    # Return solution information
    return Dict(
        "model" => model,
        "status" => status,
        "has_solution" => final,
        "objective" => final ? objective_value(model) : Inf,
        "facilities" => final ? [j for j in 1:n if value(x[j, j]) > 0.5] : Int[],
        "f_values" => final ? [value(f[i]) > 0.5 for i in 1:n] : Bool[],
        "mip_gap" => final ? relative_gap(model) : Inf,
        "solve_time" => solve_time(model)
    )
end

"""
Process B: Use Teitz-Bart heuristic to find initial h_i values
"""
function process_b_optimized(distances::Matrix{Float64}, p::Int,
    demands::Vector{Float64}, initial_h::Int, delta::Int)
    n = size(distances, 1)
    h_values = fill(initial_h, n)

    # Pre-allocate arrays
    new_facilities = Vector{Int}(undef, p)
    assignments = Vector{Int}(undef, n)
    beyond_h = BitVector(undef, n)
    needs_increase = BitVector(undef, n)

    while true
        println("\nProcess B iteration with h_i values range: ", extrema(h_values))
        # Add check to prevent h_values from exceeding n-1
        if any(h_values .>= n - 1)
            println("Warning: h_values reaching maximum allowed value (n-1)")
            h_values = min.(h_values, n - 1)
            break
        end

        H = get_closest_sites_optimized(distances, h_values)
        facilities = randperm(n)[1:p]
        improved = true

        while improved
            improved = false

            fill!(assignments, 0)
            fill!(beyond_h, false)

            @inbounds for i in 1:n
                min_dist = Inf
                best_j = 0
                for j in facilities
                    if distances[i, j] < min_dist && H[i, j]
                        min_dist = distances[i, j]
                        best_j = j
                    end
                end
                if best_j == 0
                    beyond_h[i] = true
                else
                    assignments[i] = best_j
                end
            end

            current_beyond_count = sum(beyond_h)

            @inbounds for j in setdiff(1:n, facilities)
                for k in 1:p
                    copyto!(new_facilities, facilities)
                    new_facilities[k] = j

                    new_beyond = 0
                    for i in 1:n
                        has_close = false
                        for f in new_facilities
                            if H[i, f]
                                has_close = true
                                break
                            end
                        end
                        if !has_close
                            new_beyond += 1
                        end
                    end

                    if new_beyond < current_beyond_count
                        copyto!(facilities, new_facilities)
                        improved = true
                    end
                end
            end
        end

        fill!(needs_increase, false)
        @inbounds for i in 1:n
            has_close = false
            for j in facilities
                if H[i, j]
                    has_close = true
                    break
                end
            end
            needs_increase[i] = !has_close
        end

        if !any(needs_increase)
            println("Process B converged!")
            break
        end

        @inbounds for i in 1:n
            if needs_increase[i]
                h_values[i] = min(h_values[i] + delta, n - 1)
            end
        end
    end

    # Final increase before Process A as documented in paper
    h_values .= min.(h_values .+ delta, n - 1)

    return h_values
end

function process_a(distances::Matrix{Float64}, p::Int, h_values::Vector{Int},
    demands::Vector{Float64}, delta;
    time_limit_per_model::Float64=180.0,
    total_time_limit::Float64=3600.0,
    max_gap::Float64=1e-3,
    gap_tolerance::Float64=1e-3)

    n = size(distances, 1)
    iteration = 1
    start_time = time()
    best_solution = nothing
    best_gap = Inf

    while true
        if time() - start_time >= total_time_limit
            println("Reached total time limit!")
            break
        end
        println("\nProcess A iteration $iteration")
        println("Current h_i range: ", extrema(h_values))

        remaining_time = min(time_limit_per_model,
            total_time_limit - (time() - start_time))

        result = solve_beamr_model(distances, p, h_values, demands,
            time_limit=remaining_time,
            gap_tolerance=gap_tolerance)

        if !result["has_solution"]
            println("No feasible solution found!")
            break
        end

        # Calculate actual objective value
        actual_wd = 0.0
        for i in 1:n
            min_dist = minimum(distances[i, j] for j in result["facilities"])
            actual_wd += demands[i] * min_dist
        end

        # Calculate gap using robust function
        gap = calculate_gap(result["objective"], actual_wd)



        # Track best solution
        if gap < best_gap
            best_gap = gap
            best_solution = Dict(
                "objective" => actual_wd,
                "facilities" => result["facilities"],
                "h_values" => copy(h_values),
                "is_exact" => !any(result["f_values"]),
                "gap" => gap,
                "iterations" => iteration
            )
        end

        println("Current solution:")
        println("- Lower bound from BEAMR: ", round(result["objective"], digits=4))
        println("- Actual weighted distance: ", round(actual_wd, digits=4))
        println("- Gap: ", round(gap * 100, digits=4), "%")
        println("- Nodes with f_i > 0: ", sum(result["f_values"]))

        if !any(result["f_values"])
            println("Found exact solution - all f_i = 0!")
            return best_solution
        end

        if gap <= max_gap
            println("Reached acceptable gap!")
            return best_solution
        end

        num_increased = 0
        for i in 1:n
            if result["f_values"][i]
                h_values[i] = min(h_values[i] + delta, n - 1)
                num_increased += 1
            end
        end
        println("Increased h_i for $num_increased nodes")

        iteration += 1
    end

    return best_solution !== nothing ? best_solution : Dict(
        "status" => "no_solution_found",
        "total_time" => time() - start_time
    )
end


"""
Get the h_i closest sites to each demand node i
Returns a vector of vectors where H[i] contains indices of h_i closest sites to node i
"""
function get_closest_sites_optimized(distances::Matrix{Float64}, h_values::Vector{Int64})
    n = size(distances, 1)
    H = BitMatrix(zeros(Bool, n, n))

    @inbounds for i in 1:n
        # Get indices sorted by distance
        sorted_indices = sortperm(distances[i, :])
        # Take h_i closest sites
        for k in 1:h_values[i]
            H[i, sorted_indices[k]] = true
        end
    end
    return H
end

"""
Main BEAMR solver implementing both Process A and B
"""
function solve_pmedian_beamr(distances::Matrix{Float64}, p::Int,
    demands::Vector{Float64}=ones(size(distances, 1));
    initial_h::Int=5, delta::Int=5, max_gap::Float64=0.001)
    n = size(distances, 1)
    gap_tolerance = 1e-3
    println("\nStarting BEAMR algorithm for $n nodes, p=$p")

    # Process B: Get initial h_i values
    h_values = process_b_optimized(distances, p, demands, initial_h, delta)

    # Process A: Iteratively refine solution
    # Process A with time limits
    solution = process_a(distances, p, h_values, demands, delta,
        time_limit_per_model=180.0,
        total_time_limit=3600.0,
        max_gap=max_gap,
        gap_tolerance=gap_tolerance)
end
"""
Equations from paper Section 2:
Min Z = ∑(i=1 to n)∑(j=1 to n) ai*dij*xij      [eq. 1]
s.t.
∑(j=1 to n) xij = 1                  for each i [eq. 2]
∑(j=1 to n) xjj = p                            [eq. 3]
xij ≤ xjj        for each i,j where i≠j        [eq. 4]
xij = 0,1        for each i,j                  [eq. 5]
"""
function solve_pmedian_classic(distances::Matrix{Float64}, p::Int,
    demands::Vector{Float64}=ones(size(distances, 1)))
    n = size(distances, 1)

    # Create model
    model = Model(Gurobi.Optimizer)
    set_time_limit_sec(model, 3600.0) # 30 minutos
    set_silent(model)

    # Variables
    @variable(model, x[1:n, 1:n], Bin)

    # Objective (equation 1)
    @objective(model, Min,
        sum(demands[i] * distances[i, j] * x[i, j] for i in 1:n, j in 1:n))

    # Assignment constraints (equation 2)
    for i in 1:n
        @constraint(model, sum(x[i, j] for j in 1:n) == 1)
    end

    # Number of facilities constraint (equation 3)
    @constraint(model, sum(x[j, j] for j in 1:n) == p)

    # Facility-assignment constraints (equation 4)
    for i in 1:n
        for j in 1:n
            if i != j
                @constraint(model, x[i, j] <= x[j, j])
            end
        end
    end

    optimize!(model)

    return Dict(
        "objective" => objective_value(model),
        "facilities" => [j for j in 1:n if value(x[j, j]) > 0.5],
        "model_stats" => Dict(
            "num_variables" => num_variables(model),
            #"num_constraints" => num_constraints(model),
            "solve_time" => solve_time(model)
        )
    )
end

function read_pmedian_file(filename::String)
    # Read all lines from the file
    lines = readlines(filename)

    # Parse first line for problem parameters
    n, m, p = parse.(Int, split(lines[1]))

    # Initialize cost matrix with Inf
    c = fill(Inf, n, n)

    # Set diagonal elements to 0
    for i in 1:n
        c[i, i] = 0
    end

    # Process edge data
    for line in lines[2:end]
        # Parse each edge line
        i, j, cost = parse.(Int, split(line))
        # Set costs both ways (symmetric)
        c[i, j] = cost
        c[j, i] = cost
    end

    # Apply Floyd's algorithm
    floyd_warshall!(c)

    return n, p, c
end

function floyd_warshall!(dist::Matrix{Float64})
    n = size(dist, 1)

    for k in 1:n
        for i in 1:n
            for j in 1:n
                if dist[i, k] != Inf && dist[k, j] != Inf
                    new_dist = dist[i, k] + dist[k, j]
                    if new_dist < dist[i, j]
                        dist[i, j] = new_dist
                    end
                end
            end
        end
    end
end

# Example usage
function process_pmedian_example(filename::String)
    # Read and process the file
    n, p, cost_matrix = read_pmedian_file(filename)

    println("Problem parameters:")
    println("Number of vertices: $n")
    println("Number of facilities to locate: $p")
    println("\nFinal distance matrix (first 5x5 corner shown):")

    # Display a small portion of the result
    display(cost_matrix[1:min(5, n), 1:min(5, n)])

    return n, p, cost_matrix
end

"""
Compare classic p-median formulation with BEAMR
"""
function compare_pmedian_formulations(distances::Matrix{Float64}, p::Int,
    demands::Vector{Float64}=ones(size(distances, 1));
    initial_h::Int=5, delta::Int=5)
    n = size(distances, 1)
    println("\nComparing p-median formulations for $n nodes, p=$p")

    # Solve classic formulation
    println("\nSolving classic ReVelle-Swain formulation...")
    t1 = time()
    classic_sol = solve_pmedian_classic(distances, p, demands)
    classic_time = time() - t1

    # Solve BEAMR
    println("\nSolving BEAMR formulation...")
    t2 = time()
    beamr_sol = solve_pmedian_beamr(distances, p, demands,
        initial_h=initial_h, delta=delta, max_gap=0.001)
    beamr_time = time() - t2

    # Calculate model size reduction
    classic_vars = classic_sol["model_stats"]["num_variables"]
    #classic_cons = classic_sol["model_stats"]["num_constraints"]

    # Print comparison
    println("\nResults Comparison:")
    println("==================")
    println("\nObjective Values:")
    println("Classic: ", round(classic_sol["objective"], digits=2))
    println("BEAMR:  ", round(beamr_sol["objective"], digits=2))

    println("\nModel Sizes:")
    println("Classic: $classic_vars variables, constraints")
    println("BEAMR:  Varies by iteration (see Process A output above)")

    println("\nSolution Times:")
    println("Classic: ", round(classic_time, digits=2), " seconds")
    println("BEAMR:  ", round(beamr_time, digits=2), " seconds")

    println("\nFacility Locations:")
    println("Classic: ", sort(classic_sol["facilities"]))
    println("BEAMR:  ", sort(beamr_sol["facilities"]))

    # Compare solutions
    classic_set = Set(classic_sol["facilities"])
    beamr_set = Set(beamr_sol["facilities"])

    if classic_set == beamr_set
        println("\nSolutions are identical!")
    else
        println("\nSolutions differ! Number of different locations: ",
            length(symdiff(classic_set, beamr_set)))
    end

    return Dict(
        "classic" => classic_sol,
        "beamr" => beamr_sol,
        "times" => Dict("classic" => classic_time, "beamr" => beamr_time)
    )
end

# Example usage
function test_beamr()
    # Read pmed file
    n, p, distances = read_pmedian_file(ARGS[1])
    println("Matrix dimensions: ", size(distances))  # Will print (100, 100)
    println("First 5x5 submatrix:")
    println(distances[1:5, 1:5])
    #println(distances)
    #writedlm("distanciaaaas.txt", [distances])
    demands = ones(n)

    # Solve with BEAMR
    results = compare_pmedian_formulations(distances, p, demands)
    return results

    return solution
end
test_beamr()