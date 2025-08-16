# run_y_initialization_comparison.jl
# Experiment 1: Compare different Y initialization methods with full local search

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
include("heuristic_common_functions.jl")  # Shared functions between experiments

# Data structure for tracking results
struct InitializationResult
    instance_id::String
    instance_number::Int
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
    remaining_violations::Int
    # Local search
    ls_iterations::Int
    total_ls_improvement::Float64
    successful_moves_simple::Int
    successful_moves_interchange::Int
    successful_moves_deactivate::Int
    # Status
    status::String
    error_message::String
end

# Process single instance with one Y initialization method
function process_instance_initialization(instance, instance_file, instance_number, method, config, output_dirs, experiment_log)
    instance_name = splitext(basename(instance_file))[1]
    B, S, P = instance.B, instance.S, instance.P
    size_str = "$(B)_$(S)_$(P)"
    
    println("\n--- Processing instance #$instance_number: $instance_name with method: $method ---")
    
    # Initialize result
    result = InitializationResult(
        instance_name, instance_number, method, B, S, P,
        0.0, 0.0, 0.0, 0.0, 0.0,  # timing
        0.0, 0.0, 0.0, 0.0,       # objectives
        0, 0.0,                    # splits
        0, 0, 0,                   # violations
        0, 0.0, 0, 0, 0,          # local search
        "failed", ""              # status
    )
    
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
            return InitializationResult(
                instance_name, instance_number, method, B, S, P,
                init_time, 0.0, 0.0, 0.0, init_time,
                0.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 0, 0, 0,
                "failed_y_init", "Y initialization failed"
            )
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
            return InitializationResult(
                instance_name, instance_number, method, B, S, P,
                init_time, time() - phase2_start, 0.0, 0.0, time() - total_start_time,
                0.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 0, 0, 0,
                "failed_phase2", "Phase 2 optimization failed"
            )
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
        repair_start = time()
        x_repaired = repair_activity_and_risk_constraints(instance, x_binary, Y)
        repair_time = time() - repair_start
        
        # Count remaining violations
        final_activity_viol, final_risk_viol = count_violations(instance, x_repaired, Y)
        
        post_repair_obj = dot(x_repaired, instance.D)
        
        # Create solution for local search
        sol_before_ls = Solution(instance, x_repaired, Y, post_repair_obj, 1)
        
        # Calculate targets
        targets_lower, targets_upper = calculate_targets_optimized(instance)
        
        # Local search with ALL moves enabled
        ls_config = Dict(
            "enable_simple_move" => true,
            "enable_interchange_move" => true,
            "enable_deactivate_move" => true,
            "max_iterations" => config["local_search"]["max_iterations"]
        )
        
        sol_after_ls, ls_time, ls_iters, total_improvement, succ_simple, succ_interchange, succ_deactivate = 
            local_search_with_tracking(sol_before_ls, targets_lower, targets_upper, ls_config)
        
        post_ls_obj = sol_after_ls.Weight
        
        # Save solution
        sol_dir = joinpath(output_dirs["solutions"], size_str, method)
        if !isdir(sol_dir)
            mkpath(sol_dir)
        end
        sol_file = joinpath(sol_dir, "sol_$(instance_name)_$(method).jld2")
        write_solution(sol_after_ls, sol_file)
        
        total_time = time() - total_start_time
        
        # Create result
        result = InitializationResult(
            instance_name, instance_number, method, B, S, P,
            init_time, phase2_time, repair_time, ls_time, total_time,
            phase2_obj, split_resolved_obj, post_repair_obj, post_ls_obj,
            num_splits, prop_splits,
            initial_activity_viol, initial_risk_viol,
            final_activity_viol + final_risk_viol,
            ls_iters, total_improvement, succ_simple, succ_interchange, succ_deactivate,
            "success", ""
        )
        
        # Log result
        log_message = """
        Completed instance #$instance_number ($instance_name) with $method:
          Phase 2 obj: $phase2_obj
          Split BUs: $(round(prop_splits * 100, digits=2))%
          Post repair obj: $post_repair_obj
          Post LS obj: $post_ls_obj (improvement: $total_improvement)
          Total time: $(round(total_time, digits=2))s
        """
        println(log_message)
        open(experiment_log, "a") do f
            write(f, log_message * "\n")
        end
        
    catch e
        error_msg = sprint(showerror, e)
        println("ERROR processing instance #$instance_number ($instance_name) with $method: $error_msg")
        
        result = InitializationResult(
            instance_name, instance_number, method, B, S, P,
            0.0, 0.0, 0.0, 0.0, time() - total_start_time,
            0.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 0, 0, 0,
            "error", error_msg
        )
    end
    
    return result
end

# Main experiment runner
function run_y_initialization_experiment(config_file="config_y_init.toml")
    # Load configuration
    config = TOML.parsefile(config_file)
    
    # Create experiment directory
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    experiment_name = get(config["general"], "experiment_name", "y_initialization_comparison")
    experiment_id = "exp_$(timestamp)_$(experiment_name)"
    experiment_dir = joinpath(config["general"]["base_path"], "experiments", experiment_id)
    
    # Create directory structure
    output_dirs = Dict(
        "experiment" => experiment_dir,
        "results" => joinpath(experiment_dir, "results"),
        "instances" => joinpath(experiment_dir, "instances"),
        "solutions" => joinpath(experiment_dir, "solutions"),
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
        Y Initialization Comparison Experiment
        ID: $experiment_id
        Started: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        Configuration: $config_file
        
        """)
    end
    
    # Get list of ALL instance files
    input_dir = config["general"]["input_directory"]
    all_instance_files = String[]
    
    for (root, dirs, files) in walkdir(input_dir)
        for file in files
            if endswith(file, ".jld2")
                push!(all_instance_files, joinpath(root, file))
            end
        end
    end
    
    if isempty(all_instance_files)
        error("No instance files found in $input_dir")
    end
    
    sort!(all_instance_files)
    
    # Filter instances based on IDs if specified
    instance_ids = get(config["general"], "instance_ids", Int[])
    
    if !isempty(instance_ids)
        # Extract instance number from filename (e.g., inst_1_... -> 1)
        filtered_files = String[]
        for file in all_instance_files
            # Try to extract instance number from filename
            match_result = match(r"inst_(\d+)_", basename(file))
            if match_result !== nothing
                inst_num = parse(Int, match_result[1])
                if inst_num in instance_ids
                    push!(filtered_files, file)
                end
            end
        end
        instance_files = filtered_files
        println("Processing $(length(instance_files)) instances (filtered from $(length(all_instance_files)) total)")
        println("Instance IDs: ", instance_ids)
    else
        instance_files = all_instance_files
        println("Processing all $(length(instance_files)) instances")
    end
    
    if isempty(instance_files)
        error("No instances matched the specified IDs")
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
    if get(config["methods"], "test_multi_p_median", true)
        push!(methods_to_test, "multi_p_median")
    end
    if get(config["methods"], "test_multi_p_median_benders", true)
        push!(methods_to_test, "multi_p_median_benders")
    end
    
    println("Testing methods: ", join(methods_to_test, ", "))
    
    # Results storage
    all_results = InitializationResult[]
    
    # Process each instance
    for (idx, instance_file) in enumerate(instance_files)
        instance_name = basename(instance_file)
        
        # Extract instance number for tracking
        match_result = match(r"inst_(\d+)_", instance_name)
        inst_num = match_result !== nothing ? parse(Int, match_result[1]) : idx
        
        println("\n=== Processing instance $idx/$(length(instance_files)): #$inst_num - $instance_name ===")
        
        # Load instance
        instance = read_instance(instance_file)
        
        # Copy instance file
        size_str = "$(instance.B)_$(instance.S)_$(instance.P)"
        inst_dir = joinpath(output_dirs["instances"], size_str)
        mkpath(inst_dir)
        cp(instance_file, joinpath(inst_dir, instance_name))
        
        # Test each method
        for method in methods_to_test
            result = process_instance_initialization(instance, instance_file, inst_num, method, config, output_dirs, experiment_log)
            push!(all_results, result)
            
            # Save intermediate results after each method
            save_initialization_results(all_results, output_dirs["results"])
        end
    end
    
    # Generate final summary
    generate_initialization_summary(all_results, output_dirs["results"], experiment_log)
    
    # Log completion
    completion_msg = "\n=== Experiment Complete at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) ==="
    println(completion_msg)
    open(experiment_log, "a") do f
        write(f, completion_msg * "\n")
    end
    
    println("\nResults saved to: $experiment_dir")
end

# Save results to CSV and JSON
function save_initialization_results(results::Vector{InitializationResult}, results_dir::String)
    # Convert to DataFrame for CSV
    df = DataFrame(
        instance_id = [r.instance_id for r in results],
        instance_number = [r.instance_number for r in results],
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
        remaining_violations = [r.remaining_violations for r in results],
        ls_iterations = [r.ls_iterations for r in results],
        total_ls_improvement = [r.total_ls_improvement for r in results],
        successful_moves_simple = [r.successful_moves_simple for r in results],
        successful_moves_interchange = [r.successful_moves_interchange for r in results],
        successful_moves_deactivate = [r.successful_moves_deactivate for r in results],
        status = [r.status for r in results],
        error_message = [r.error_message for r in results]
    )
    
    CSV.write(joinpath(results_dir, "y_initialization_results.csv"), df)
    
    # Save as JSON for detailed analysis
    open(joinpath(results_dir, "y_initialization_results.json"), "w") do f
        JSON.print(f, results, 2)
    end
end

# Generate summary statistics
function generate_initialization_summary(results::Vector{InitializationResult}, results_dir::String, log_file::String)
    # Filter successful results
    successful_results = filter(r -> r.status == "success", results)
    
    if isempty(successful_results)
        println("Warning: No successful results to summarize")
        return
    end
    
    # Group by method
    methods = unique([r.method for r in successful_results])
    
    summary_df = DataFrame()
    
    for method in methods
        method_results = filter(r -> r.method == method, successful_results)
        
        if !isempty(method_results)
            summary_row = DataFrame(
                method = method,
                num_instances = length(method_results),
                avg_total_time = mean([r.total_time for r in method_results]),
                std_total_time = std([r.total_time for r in method_results]),
                avg_init_time = mean([r.init_time for r in method_results]),
                avg_phase2_time = mean([r.phase2_time for r in method_results]),
                avg_repair_time = mean([r.repair_time for r in method_results]),
                avg_ls_time = mean([r.ls_time for r in method_results]),
                avg_post_ls_obj = mean([r.post_ls_obj for r in method_results]),
                std_post_ls_obj = std([r.post_ls_obj for r in method_results]),
                avg_ls_improvement = mean([r.total_ls_improvement for r in method_results]),
                avg_ls_iterations = mean([r.ls_iterations for r in method_results]),
                avg_prop_split_bus = mean([r.prop_split_bus for r in method_results]),
                avg_violations_remaining = mean([r.remaining_violations for r in method_results]),
                success_rate = length(method_results) / length(filter(r -> r.method == method, results))
            )
            summary_df = vcat(summary_df, summary_row)
        end
    end
    
    # Sort by average objective (best first)
    sort!(summary_df, :avg_post_ls_obj)
    
    CSV.write(joinpath(results_dir, "y_initialization_summary.csv"), summary_df)
    
    # Identify best method
    if nrow(summary_df) > 0
        best_method = summary_df[1, :method]
        best_obj = summary_df[1, :avg_post_ls_obj]
        
        # Write recommendation for next experiment
        recommendation = Dict(
            "best_method" => best_method,
            "avg_objective" => best_obj,
            "recommendation" => "Use '$best_method' for local search ablation experiment"
        )
        
        open(joinpath(results_dir, "best_method.json"), "w") do f
            JSON.print(f, recommendation, 2)
        end
    end
    
    # Write summary to log
    summary_text = "\n\n=== SUMMARY STATISTICS ===\n"
    for row in eachrow(summary_df)
        summary_text *= "\nMethod: $(row.method)\n"
        summary_text *= "  Success rate: $(round(row.success_rate * 100, digits=1))%\n"
        summary_text *= "  Avg time: $(round(row.avg_total_time, digits=2))s (std: $(round(row.std_total_time, digits=2)))\n"
        summary_text *= "  Avg final obj: $(round(row.avg_post_ls_obj, digits=2)) (std: $(round(row.std_post_ls_obj, digits=2)))\n"
        summary_text *= "  Avg LS improvement: $(round(row.avg_ls_improvement, digits=2))\n"
        summary_text *= "  Avg split proportion: $(round(row.avg_prop_split_bus * 100, digits=2))%\n"
    end
    
    if nrow(summary_df) > 0
        summary_text *= "\n\n=== RECOMMENDATION ===\n"
        summary_text *= "Best method: $(summary_df[1, :method]) with avg objective $(round(summary_df[1, :avg_post_ls_obj], digits=2))\n"
    end
    
    println(summary_text)
    open(log_file, "a") do f
        write(f, summary_text)
    end
end

# Main function
function main()
    if length(ARGS) == 0
        println("Using default config file: config_y_init.toml")
        config_file = "config_y_init.toml"
    else
        config_file = ARGS[1]
    end
    
    run_y_initialization_experiment(config_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end