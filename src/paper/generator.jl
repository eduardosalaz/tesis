using Distances, JLD2, Random, TOML
using DelimitedFiles  # For writedlm
using Types

function read_config(config_file="config.toml")
    isfile(config_file) || error("Configuration file not found: $config_file")
    return TOML.parsefile(config_file)
end

function generate_coords(B, S)
    BU_coords = rand(5:10000, (B, 2))
    S_coords = rand(5:10000, (S, 2))
    return BU_coords, S_coords
end

function generate_dist(BU_coords, S_coords, B, S)
    metrica = Euclidean()
    mat = zeros(S, B)
    for i in 1:S
        for j in 1:B
            distancia = metrica(BU_coords[j, :], S_coords[i, :])
            mat[i, j] = distancia
        end
    end
    @debug "Wrote distance matrix"
    return round.(Int, mat)
end

# Generate base data that is NOT dependent on P
function generate_base_data(B, S, config)
    # Read model parameters from config
    K = config["model"]["K"]
    M = config["model"]["M"]
    
    # Generate coordinates and distance matrix
    BU_coords, S_coords = generate_coords(B, S)
    dist_mat = generate_dist(BU_coords, S_coords, B, S)
    
    # Generate activity values V (independent of P)
    V = generate_activity_values(B, M, config)
    
    # Generate risk values R (independent of P)
    R = generate_risk_values(B, config)
    
    # Generate facility types
    percentages = get(config["instance_sets"]["small"], "percentages", [0.4, 0.3, 0.2, 0.1])
    Sk = generate_facility_types(S, K, percentages)
    
    return BU_coords, S_coords, dist_mat, V, R, Sk, percentages
end

# Generate activity values (independent of P)
function generate_activity_values(B, M, config)
    # Get ranges from config or use defaults
    ranges = [
        get(config["activity_ranges"], "act1_range", [1, 10]),
        get(config["activity_ranges"], "act2_range", [1000, 10000]),
        get(config["activity_ranges"], "act3_range", [1000, 5000])
    ]
    
    # Initialize V with M arrays, each of length B
    V = [zeros(Int64, B) for _ in 1:M]
    
    # Fill each array in V with random numbers according to their ranges
    for m in 1:M
        if m <= length(ranges)
            min_val, max_val = ranges[m]
            for b in 1:B
                V[m][b] = rand(min_val:max_val)
            end
        end
    end
    
    return V
end

# Generate risk values (independent of P)
function generate_risk_values(B, config)
    min_risk, max_risk = get(config["risk_ranges"], "risk_range", [30, 60])
    return rand(min_risk:max_risk, B)
end

# Generate facility types (independent of P)
function generate_facility_types(S, K, percentages)
    # Calculate how many facilities should be of each type
    counts_per_type = round.(Int, S .* percentages)
    
    # Adjust for rounding errors to ensure sum equals S
    while sum(counts_per_type) != S
        if sum(counts_per_type) < S
            idx = argmax(percentages .- counts_per_type/S)
            counts_per_type[idx] += 1
        else
            idx = argmin(percentages .- counts_per_type/S)
            counts_per_type[idx] -= 1
        end
    end
    
    # Create the assignment vector
    Ks = Int64[]
    for k in 1:K
        append!(Ks, fill(k, counts_per_type[k]))
    end
    shuffle!(Ks)  # Randomize the order
    
    # Create Sk (facilities of each type)
    Sk = [findall(x -> x == k, Ks) for k in 1:K]
    
    return Sk
end

# Generate P-dependent parameters
function generate_p_dependent_params(B, S, P, V, R, Sk, percentages, config)
    # Calculate Lk and Uk (depend on P)
    Lk, Uk = calculate_bounds(P, percentages, S, Sk)
    
    # Calculate μ values (depend on P)
    μ = calculate_mu_values(S, P, V, config["model"]["M"])
    
    # Calculate β values (depend on P)
    β = calculate_beta_values(S, P, R, B, config)
    
    # Get tolerance values from config
    T = config["model"]["T"]
    
    return Lk, Uk, μ, β, T
end

# Calculate Lk and Uk bounds (dependent on P)
function calculate_bounds(P, percentages, S, Sk)
    # Calculate bounds based on P
    # For each type k, we want between (percentage-0.05)*p and (percentage+0.05)*p centers
    Lk = round.(Int, (percentages .- 0.05) .* P)
    Uk = round.(Int, (percentages .+ 0.05) .* P)
    
    # Ensure Lk isn't negative and Uk doesn't exceed available facilities
    counts_per_type = [length(s) for s in Sk]
    Lk = max.(0, Lk)
    Uk = min.(Uk, counts_per_type)
    
    return Lk, Uk
end

# Calculate μ values (dependent on P)
function calculate_mu_values(S, P, V, M)
    # Initialize μ with M arrays, each of length S
    μ = [zeros(Int64, S) for _ in 1:M]
    
    # Calculate μ values based on sum of corresponding V array and P
    for m in 1:M
        sum_vals = sum(V[m])
        base_value = round(Int, sum_vals / P)
        for s in 1:S
            μ[m][s] = base_value
        end
    end
    
    return μ
end

# Calculate β values (dependent on P)
function calculate_beta_values(S, P, R, B, config)
    min_risk, max_risk = get(config["risk_ranges"], "risk_range", [30, 60])
    τ = config["model"]["tau"]
    β = zeros(Int64, S)
    
    for s in 1:S
        base_value = round(Int, (((min_risk + max_risk) / 2 * B) / P) * (1 + τ))
        β[s] = base_value
    end
    
    return β
end

# Create complete instance with base data and P-dependent parameters
function create_instance(B, S, P, config, base_data=nothing)
    # If base data is not provided, generate it
    if base_data === nothing
        BU_coords, S_coords, dist_mat, V, R, Sk, percentages = generate_base_data(B, S, config)
    else
        BU_coords, S_coords, dist_mat, V, R, Sk, percentages = base_data
    end
    
    # Generate P-dependent parameters
    Lk, Uk, μ, β, T = generate_p_dependent_params(B, S, P, V, R, Sk, percentages, config)
    
    # Create instance
    K = config["model"]["K"]
    M = config["model"]["M"]
    
    instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat, Sk, Lk, Uk, V, μ, T, R, β)
    
    return instance
end

# Main function to process instances based on configuration
function generate_instances(config_file="config.toml")
    config = read_config(config_file)
    base_path = get(config["general"], "base_path", ".")
    prefix = get(config["general"], "prefix", "005_newk")
    
    for (set_name, set_config) in config["instance_sets"]
        B = set_config["B"]
        S = set_config["S"]
        
        # Generate base data (independent of P)
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        
        # Create an instance for each P value
        for P in set_config["P_values"]
            size_str = "$(B)_$(S)_$(P)"
            println("Generating instance with B=$B, S=$S, P=$P")
            
            # Create directory structure
            inst_dir_path = joinpath(base_path, "instances", size_str)
            if !isdir(inst_dir_path)
                mkpath(inst_dir_path)
            end
            
            # Generate instance parameters specific to this P value
            parameters = generate_params(B, S, P)
            
            # Create instance
            K = 4  # Fixed values as in original code
            M = 3
            instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat, parameters...)
            
            # Save instance
            file_name = "inst_1_$(prefix)_$(size_str).jld2"
            full_path = joinpath(inst_dir_path, file_name)
            write_instance(instance, full_path)
            
            # Plot instance
            plot_name = replace(file_name, ".jld2" => ".png")
            full_plot = joinpath(inst_dir_path, plot_name)
            plot_instance(instance, full_plot)
            
            println("Saved instance to $full_path")
        end
    end
    println("Instance generation complete.")
end

# Function that maintains compatibility with the original interface
function generate_solve(size_str, number_str; config_file="config.toml")
    config = read_config(config_file)
    K = 4
    M = 3
    
    # Parse size string and number
    B, S, P = parse.(Int, split(size_str, "_"))
    number = parse(Int, number_str)
    println("Processing: B=$B, S=$S, P=$P, instances=$number")
    
    # Get configuration values
    base_path = get(config["general"], "base_path", ".")
    prefix = get(config["general"], "prefix", "005_newk")
    time_limit = get(config["general"], "time_limit", 1800.0)
    method = get(config["general"], "method", 1)
    
    # Create directory structure
    inst_dir_path = joinpath(base_path, "insts_new", size_str)
    failed_dir_path = joinpath(base_path, "insts_new", size_str, "failed_instances")
    sol_dir_path = joinpath(base_path, "out_new", "solutions", size_str)
    plot_out_dir_path = joinpath(base_path, "out_new", "plots", size_str)
    plot_sol_dir_path = joinpath(plot_out_dir_path, "solutions")
    plot_inst_dir_path = joinpath(plot_out_dir_path, "instances")
    logs_dir_path = joinpath(base_path, "logs_new")
    
    # Create directories if they don't exist
    for dir in [inst_dir_path, failed_dir_path, sol_dir_path, 
                plot_out_dir_path, plot_sol_dir_path, plot_inst_dir_path, logs_dir_path]
        if !isdir(dir)
            mkpath(dir)
        end
    end
    
    # Generate instances
    for i in 1:number
        # Generate base data (independent of P)
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        
        # Generate parameters (dependent on P)
        parameters = generate_params(B, S, P)
        
        # Create instance
        instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat, parameters...)
        println("Built instance $i")
        
        # Build and solve model
        model = build_model(instance)
        println("Built model $i")
        
        # Configure log files
        log_file = joinpath(logs_dir_path, "log_$(size_str)_$(prefix)_$(time_limit)s_$i.txt")
        results_file = joinpath(logs_dir_path, "results_$(size_str)_$(prefix)_$(time_limit)s_$i.txt")
        
        # Solve model
        X, Y, obj_val, solve_time = optimize_model(model, i, log_file, results_file, time_limit, method)
        
        # File paths
        file_inst_path = "inst_$(i)_$(prefix)_$(size_str).jld2"
        file_sol_path = "sol_$(i)_$(prefix)_$(size_str).jld2"
        
        if obj_val == 0
            @error "Instance $i not solved"
            # Save failed instance
            full_failed_inst_path = joinpath(failed_dir_path, file_inst_path)
            write_instance(instance, full_failed_inst_path)
        else
            # Save solution
            solution = Solution(instance, X, Y, obj_val, solve_time)
            full_sol_path = joinpath(sol_dir_path, file_sol_path)
            full_inst_path = joinpath(inst_dir_path, file_inst_path)
            
            # Plot paths
            plot_sol_path = joinpath(plot_sol_dir_path, "sol_$(i)_$(prefix)_$(size_str).png")
            plot_inst_path = joinpath(plot_inst_dir_path, "inst_$(i)_$(prefix)_$(size_str).png")
            
            # Save everything
            write_instance(instance, full_inst_path)
            write_solution(solution, full_sol_path)
            plot_instance(instance, plot_inst_path)
            plot_solution(solution, plot_sol_path)
        end
    end
    println("Processing complete for $size_str")
end

# Entry point that handles command line arguments
function main()
    if length(ARGS) == 0
        println("Usage:")
        println("  julia generator.jl config.toml                  # Generate all instances from config")
        println("  julia generator.jl size_str number [config.toml] # Generate and solve specific instances")
        println("Example: julia generator.jl 625_200_30 5 config.toml")
        return
    end
    
    if length(ARGS) >= 2
        # Traditional interface: size_str and number
        size_str = ARGS[1]
        number = ARGS[2]
        config_file = length(ARGS) >= 3 ? ARGS[3] : "config.toml"
        generate_solve(size_str, number, config_file=config_file)
    else
        # Config-based interface
        config_file = ARGS[1]
        generate_instances(config_file)
    end
end

# Uncomment to run when script is executed directly
# main()