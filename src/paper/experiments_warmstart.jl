# run_warmstart_experiment.jl
# Experiment: Heuristic + Gurobi Warmstart
# Runs constructive heuristic (multi_p_median_benders + all LS moves),
# then uses solution as warmstart for Gurobi with adjusted time limit

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
include("heuristic_common_functions.jl")

# Data structure for tracking results
struct WarmstartResult
    instance_id::String
    instance_number::Int
    B::Int
    S::Int
    P::Int
    # Heuristic phase timing
    y_init_time::Float64
    phase2_time::Float64
    repair_time::Float64
    ls_time::Float64
    total_heuristic_time::Float64
    # Heuristic objectives
    phase2_obj::Float64
    split_resolved_obj::Float64
    post_repair_obj::Float64
    heuristic_obj::Float64  # Final heuristic objective (post-LS)
    # Heuristic details
    num_split_bus::Int
    prop_split_bus::Float64
    ls_iterations::Int
    ls_improvement::Float64
    # Gurobi phase
    gurobi_time_limit::Float64
    gurobi_actual_time::Float64
    gurobi_obj::Float64
    gurobi_bound::Float64
    gurobi_gap::Float64
    gurobi_status::String
    # Improvement analysis
    total_time::Float64
    warmstart_improvement::Float64  # heuristic_obj - gurobi_obj (positive = improvement)
    warmstart_improvement_pct::Float64
    # Status
    heuristic_status::String
    overall_status::String
    error_message::String
end

# Run complete heuristic (construction + LS)
function run_heuristic(instance, instance_name, config, logs_dir)
    """
    Runs the complete heuristic: Y initialization, phase 2, repair, and local search.
    Returns (solution, Y, stats) or (nothing, nothing, nothing) if failed.
    """
    
    println("  Running heuristic construction...")
    
    heuristic_start = time()
    
    try
        # Phase 1: Initialize Y using multi_p_median_benders
        y_init_start = time()
        y_log_file = joinpath(logs_dir, "y_init_$(instance_name).log")
        
        Y, y_init_time, success = init_y_multi_p_median(
            instance, 
            config["time_limits"]["exact_methods"], 
            y_log_file, 
            false  # use_benders = true
        )
        
        if !success || Y === nothing
            println("    Y initialization failed")
            return nothing, nothing, (
                y_init_time = y_init_time,
                phase2_time = 0.0,
                repair_time = 0.0,
                ls_time = 0.0,
                total_heuristic_time = time() - heuristic_start,
                phase2_obj = 0.0,
                split_resolved_obj = 0.0,
                post_repair_obj = 0.0,
                heuristic_obj = 0.0,
                num_split_bus = 0,
                prop_split_bus = 0.0,
                ls_iterations = 0,
                ls_improvement = 0.0,
                status = "failed_y_init"
            )
        end
        
        # Phase 2: Solve transportation problem
        phase2_start = time()
        phase2_log_file = joinpath(logs_dir, "phase2_$(instance_name).log")
        
        model_transport = solve_phase2_model(instance, Y)
        set_optimizer(model_transport, Gurobi.Optimizer)
        set_optimizer_attribute(model_transport, "LogFile", phase2_log_file)
        set_time_limit_sec(model_transport, config["time_limits"]["phase2"])
        #set_optimizer_attribute(model_transport, "OutputFlag", 0)  # Suppress output
        
        optimize!(model_transport)
        
        if termination_status(model_transport) != MOI.OPTIMAL && 
           termination_status(model_transport) != MOI.TIME_LIMIT
            println("    Phase 2 optimization failed")
            phase2_time = time() - phase2_start
            return nothing, nothing, (
                y_init_time = y_init_time,
                phase2_time = phase2_time,
                repair_time = 0.0,
                ls_time = 0.0,
                total_heuristic_time = time() - heuristic_start,
                phase2_obj = 0.0,
                split_resolved_obj = 0.0,
                post_repair_obj = 0.0,
                heuristic_obj = 0.0,
                num_split_bus = 0,
                prop_split_bus = 0.0,
                ls_iterations = 0,
                ls_improvement = 0.0,
                status = "failed_phase2"
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
        
        # Repair constraints
        repair_start = time()
        x_repaired = repair_activity_and_risk_constraints(instance, x_binary, Y)
        repair_time = time() - repair_start
        post_repair_obj = dot(x_repaired, instance.D)
        
        # Create solution for local search
        sol_before_ls = Solution(instance, x_repaired, Y, post_repair_obj, 1)
        
        # Calculate targets
        targets_lower, targets_upper = calculate_targets_optimized(instance)
        
        # Local search with ALL moves enabled
        ls_start = time()
        ls_config = Dict(
            "enable_simple_move" => config["local_search"]["enable_simple_move"],
            "enable_interchange_move" => config["local_search"]["enable_interchange_move"],
            "enable_deactivate_move" => config["local_search"]["enable_deactivate_move"],
            "max_iterations" => config["local_search"]["max_iterations"]
        )
        
        sol_after_ls, ls_time, ls_iters, total_improvement, _, _, _ = 
            local_search_with_tracking(sol_before_ls, targets_lower, targets_upper, ls_config)
        
        heuristic_obj = sol_after_ls.Weight
        total_heuristic_time = time() - heuristic_start
        
        stats = (
            y_init_time = y_init_time,
            phase2_time = phase2_time,
            repair_time = repair_time,
            ls_time = ls_time,
            total_heuristic_time = total_heuristic_time,
            phase2_obj = phase2_obj,
            split_resolved_obj = split_resolved_obj,
            post_repair_obj = post_repair_obj,
            heuristic_obj = heuristic_obj,
            num_split_bus = num_splits,
            prop_split_bus = prop_splits,
            ls_iterations = ls_iters,
            ls_improvement = total_improvement,
            status = "success"
        )
        
        println("    Heuristic complete: obj = $(round(heuristic_obj, digits=2)), time = $(round(total_heuristic_time, digits=2))s")
        
        return sol_after_ls, Y, stats
        
    catch e
        error_msg = sprint(showerror, e)
        println("    ERROR in heuristic: $error_msg")
        return nothing, nothing, (
            y_init_time = 0.0,
            phase2_time = 0.0,
            repair_time = 0.0,
            ls_time = 0.0,
            total_heuristic_time = time() - heuristic_start,
            phase2_obj = 0.0,
            split_resolved_obj = 0.0,
            post_repair_obj = 0.0,
            heuristic_obj = 0.0,
            num_split_bus = 0,
            prop_split_bus = 0.0,
            ls_iterations = 0,
            ls_improvement = 0.0,
            status = "error"
        )
    end
end

# Solve with Gurobi using warmstart
function solve_with_warmstart(instance, heuristic_solution, Y, time_limit, config, logs_dir, instance_name)
    """
    Builds Gurobi model, applies warmstart, and solves.
    Returns (obj, bound, gap, actual_time, status, X_solution, Y_solution)
    """
    
    println("  Solving with Gurobi warmstart (time limit: $(round(time_limit, digits=2))s)...")
    
    if time_limit <= 0
        println("    Warning: No time remaining for Gurobi!")
        return (heuristic_solution.Weight, heuristic_solution.Weight, 0.0, 0.0, "no_time", 
                heuristic_solution.X, Y)
    end
    
    try
        # Build the full model
        model = build_model(instance)        
        # Apply warmstart
        X_heuristic = heuristic_solution.X
        Y_heuristic = Y
        
        # Set variable start values
        x = model[:x]
        y = model[:y]
        set_start_value.(x, X_heuristic)
        set_start_value.(y, Y_heuristic)

        # Set optimizer
        set_optimizer(model, Gurobi.Optimizer)
        
        # Configure Gurobi
        gurobi_log_file = joinpath(logs_dir, "gurobi_$(instance_name).log")
        set_optimizer_attribute(model, "LogFile", gurobi_log_file)
        set_optimizer_attribute(model, "TimeLimit", time_limit)
        set_optimizer_attribute(model, "Method", config["gurobi"]["method"])
        #set_optimizer_attribute(model, "MIPFocus", config["gurobi"]["mip_focus"])
        
        println("    Warmstart applied, optimizing...")
        
        # Solve
        gurobi_start = time()
        optimize!(model)
        gurobi_time = time() - gurobi_start
        
        # Extract results
        status = string(termination_status(model))
        obj = objective_value(model)
        bound = objective_bound(model)
        gap = abs(obj - bound) / abs(obj) * 100.0
        
        # Extract solution
        X_solution = value.(x)
        Y_solution = value.(y)
        
        println("    Gurobi complete: obj = $(round(obj, digits=2)), bound = $(round(bound, digits=2)), gap = $(round(gap, digits=2))%, time = $(round(gurobi_time, digits=2))s")
        
        return obj, bound, gap, gurobi_time, status, X_solution, Y_solution
        
    catch e
        error_msg = sprint(showerror, e)
        println("    ERROR in Gurobi solve: $error_msg")
        return (heuristic_solution.Weight, heuristic_solution.Weight, 0.0, 0.0, "error",
                heuristic_solution.X, Y)
    end
end

# Process single instance
function process_instance_warmstart(instance, instance_file, instance_number, config, output_dirs, experiment_log)
    instance_name = splitext(basename(instance_file))[1]
    B, S, P = instance.B, instance.S, instance.P
    size_str = "$(B)_$(S)_$(P)"
    
    println("\n=== Processing instance #$instance_number: $instance_name ===")
    
    total_start_time = time()
    total_time_limit = config["time_limits"]["total_time_limit"]
    
    # Run heuristic
    heuristic_solution, Y, heuristic_stats = run_heuristic(
        instance, instance_name, config, output_dirs["logs"]
    )
    
    # Check if heuristic failed
    if heuristic_solution === nothing
        println("  Heuristic failed - skipping instance")
        
        result = WarmstartResult(
            instance_name, instance_number, B, S, P,
            heuristic_stats.y_init_time,
            heuristic_stats.phase2_time,
            heuristic_stats.repair_time,
            heuristic_stats.ls_time,
            heuristic_stats.total_heuristic_time,
            heuristic_stats.phase2_obj,
            heuristic_stats.split_resolved_obj,
            heuristic_stats.post_repair_obj,
            heuristic_stats.heuristic_obj,
            heuristic_stats.num_split_bus,
            heuristic_stats.prop_split_bus,
            heuristic_stats.ls_iterations,
            heuristic_stats.ls_improvement,
            0.0, 0.0, 0.0, 0.0, 0.0, "not_run",
            heuristic_stats.total_heuristic_time,
            0.0, 0.0,
            heuristic_stats.status,
            "heuristic_failed",
            ""
        )
        
        open(experiment_log, "a") do f
            write(f, "Instance #$instance_number ($instance_name): Heuristic failed ($(heuristic_stats.status))\n")
        end
        
        return result
    end
    
    # Calculate remaining time for Gurobi
    gurobi_time_limit = total_time_limit - heuristic_stats.total_heuristic_time
    
    # Solve with warmstart
    gurobi_obj, gurobi_bound, gurobi_gap, gurobi_time, gurobi_status, X_gurobi, Y_gurobi = solve_with_warmstart(
        instance, heuristic_solution, Y, gurobi_time_limit, config, 
        output_dirs["logs"], instance_name
    )
    
    total_time = time() - total_start_time
    
    # Calculate improvements (positive = improvement for minimization)
    warmstart_improvement = heuristic_stats.heuristic_obj - gurobi_obj
    warmstart_improvement_pct = (warmstart_improvement / heuristic_stats.heuristic_obj) * 100.0
    
    # Create result
    result = WarmstartResult(
        instance_name, instance_number, B, S, P,
        heuristic_stats.y_init_time,
        heuristic_stats.phase2_time,
        heuristic_stats.repair_time,
        heuristic_stats.ls_time,
        heuristic_stats.total_heuristic_time,
        heuristic_stats.phase2_obj,
        heuristic_stats.split_resolved_obj,
        heuristic_stats.post_repair_obj,
        heuristic_stats.heuristic_obj,
        heuristic_stats.num_split_bus,
        heuristic_stats.prop_split_bus,
        heuristic_stats.ls_iterations,
        heuristic_stats.ls_improvement,
        gurobi_time_limit,
        gurobi_time,
        gurobi_obj,
        gurobi_bound,
        gurobi_gap,
        gurobi_status,
        total_time,
        warmstart_improvement,
        warmstart_improvement_pct,
        "success",
        "success",
        ""
    )
    
    # Save heuristic solution
    sol_dir = joinpath(output_dirs["solutions"], size_str, "heuristic")
    mkpath(sol_dir)
    heuristic_sol_file = joinpath(sol_dir, "heuristic_$(instance_name).jld2")
    write_solution(heuristic_solution, heuristic_sol_file)
    
    # Save Gurobi solution
    gurobi_sol_dir = joinpath(output_dirs["solutions"], size_str, "gurobi")
    mkpath(gurobi_sol_dir)
    gurobi_sol_file = joinpath(gurobi_sol_dir, "gurobi_ws_$(instance_name).jld2")
    gurobi_solution = Solution(instance, X_gurobi, Y_gurobi, gurobi_obj, 1)
    write_solution(gurobi_solution, gurobi_sol_file)
    
    # Log result
    log_message = """
    Completed instance #$instance_number ($instance_name):
      Heuristic obj: $(round(heuristic_stats.heuristic_obj, digits=2)) in $(round(heuristic_stats.total_heuristic_time, digits=2))s
      Gurobi obj: $(round(gurobi_obj, digits=2)) (bound: $(round(gurobi_bound, digits=2)), gap: $(round(gurobi_gap, digits=2))%) in $(round(gurobi_time, digits=2))s (limit: $(round(gurobi_time_limit, digits=2))s)
      Improvement: $(round(warmstart_improvement, digits=2)) ($(round(warmstart_improvement_pct, digits=2))%)
      Total time: $(round(total_time, digits=2))s
    """
    println(log_message)
    open(experiment_log, "a") do f
        write(f, log_message * "\n")
    end
    
    return result
end

# Main experiment runner
function run_warmstart_experiment(config_file="config_warmstart.toml")
    # Load configuration
    config = TOML.parsefile(config_file)
    
    # Create experiment directory
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    experiment_name = get(config["general"], "experiment_name", "warmstart_experiment")
    y_method = config["method"]["y_initialization"]
    experiment_id = "exp_$(timestamp)_$(experiment_name)_$(y_method)"
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
        Warmstart Experiment (Heuristic + Gurobi)
        Y Initialization: $(y_method)
        ID: $experiment_id
        Started: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        Configuration: $config_file
        Total time limit: $(config["time_limits"]["total_time_limit"])s per instance
        
        """)
    end
    
    # Get instance files
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
    
    # Filter instances based on IDs
    instance_ids = get(config["general"], "instance_ids", Int[])
    
    if !isempty(instance_ids)
        filtered_files = String[]
        for file in all_instance_files
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
    
    # Results storage
    all_results = WarmstartResult[]
    
    # Process each instance
    for (idx, instance_file) in enumerate(instance_files)
        instance_name = basename(instance_file)
        
        # Extract instance number
        match_result = match(r"inst_(\d+)_", instance_name)
        inst_num = match_result !== nothing ? parse(Int, match_result[1]) : idx
        
        # Load instance
        instance = read_instance(instance_file)
        
        # Copy instance file
        size_str = "$(instance.B)_$(instance.S)_$(instance.P)"
        inst_dir = joinpath(output_dirs["instances"], size_str)
        mkpath(inst_dir)
        cp(instance_file, joinpath(inst_dir, basename(instance_file)))
        
        # Process instance
        result = process_instance_warmstart(instance, instance_file, inst_num, config, output_dirs, experiment_log)
        push!(all_results, result)
        
        # Save intermediate results
        save_warmstart_results(all_results, output_dirs["results"])
    end
    
    # Generate final summary
    generate_warmstart_summary(all_results, output_dirs["results"], experiment_log)
    
    # Log completion
    completion_msg = "\n=== Experiment Complete at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) ==="
    println(completion_msg)
    open(experiment_log, "a") do f
        write(f, completion_msg * "\n")
    end
    
    println("\nResults saved to: $experiment_dir")
end

# Save results to CSV and JSON
function save_warmstart_results(results::Vector{WarmstartResult}, results_dir::String)
    # Convert to DataFrame
    df = DataFrame(
        instance_id = [r.instance_id for r in results],
        instance_number = [r.instance_number for r in results],
        B = [r.B for r in results],
        S = [r.S for r in results],
        P = [r.P for r in results],
        y_init_time = [r.y_init_time for r in results],
        phase2_time = [r.phase2_time for r in results],
        repair_time = [r.repair_time for r in results],
        ls_time = [r.ls_time for r in results],
        total_heuristic_time = [r.total_heuristic_time for r in results],
        heuristic_obj = [r.heuristic_obj for r in results],
        ls_iterations = [r.ls_iterations for r in results],
        ls_improvement = [r.ls_improvement for r in results],
        gurobi_time_limit = [r.gurobi_time_limit for r in results],
        gurobi_actual_time = [r.gurobi_actual_time for r in results],
        gurobi_obj = [r.gurobi_obj for r in results],
        gurobi_bound = [r.gurobi_bound for r in results],
        gurobi_gap = [r.gurobi_gap for r in results],
        total_time = [r.total_time for r in results],
        warmstart_improvement = [r.warmstart_improvement for r in results],
        warmstart_improvement_pct = [r.warmstart_improvement_pct for r in results],
        heuristic_status = [r.heuristic_status for r in results],
        overall_status = [r.overall_status for r in results]
    )
    
    CSV.write(joinpath(results_dir, "warmstart_results.csv"), df)
    
    # Save as JSON
    open(joinpath(results_dir, "warmstart_results.json"), "w") do f
        JSON.print(f, results, 2)
    end
end

# Generate summary statistics
function generate_warmstart_summary(results::Vector{WarmstartResult}, results_dir::String, log_file::String)
    # Filter successful results
    successful_results = filter(r -> r.overall_status == "success", results)
    
    if isempty(successful_results)
        println("Warning: No successful results to summarize")
        return
    end
    
    # Overall statistics
    summary_df = DataFrame(
        metric = String[],
        value = Float64[]
    )
    
    push!(summary_df, ("num_instances", length(successful_results)))
    push!(summary_df, ("avg_heuristic_time", mean([r.total_heuristic_time for r in successful_results])))
    push!(summary_df, ("avg_heuristic_obj", mean([r.heuristic_obj for r in successful_results])))
    push!(summary_df, ("avg_gurobi_time", mean([r.gurobi_actual_time for r in successful_results])))
    push!(summary_df, ("avg_gurobi_obj", mean([r.gurobi_obj for r in successful_results])))
    push!(summary_df, ("avg_gurobi_bound", mean([r.gurobi_bound for r in successful_results])))
    push!(summary_df, ("avg_total_time", mean([r.total_time for r in successful_results])))
    push!(summary_df, ("avg_warmstart_improvement", mean([r.warmstart_improvement for r in successful_results])))
    push!(summary_df, ("avg_warmstart_improvement_pct", mean([r.warmstart_improvement_pct for r in successful_results])))
    push!(summary_df, ("avg_final_gap", mean([r.gurobi_gap for r in successful_results])))
    
    # Count improvements (positive improvement = Gurobi better than heuristic)
    num_improved = count(r -> r.warmstart_improvement > 0.01, successful_results)
    num_same = count(r -> abs(r.warmstart_improvement) <= 0.01, successful_results)
    num_worse = count(r -> r.warmstart_improvement < -0.01, successful_results)
    
    push!(summary_df, ("num_improved", num_improved))
    push!(summary_df, ("num_same", num_same))
    push!(summary_df, ("num_worse", num_worse))
    
    CSV.write(joinpath(results_dir, "warmstart_summary.csv"), summary_df)
    
    # Time breakdown
    time_breakdown = DataFrame(
        phase = ["Y Initialization", "Phase 2", "Repair", "Local Search", "Total Heuristic", "Gurobi", "Total"],
        avg_time = [
            mean([r.y_init_time for r in successful_results]),
            mean([r.phase2_time for r in successful_results]),
            mean([r.repair_time for r in successful_results]),
            mean([r.ls_time for r in successful_results]),
            mean([r.total_heuristic_time for r in successful_results]),
            mean([r.gurobi_actual_time for r in successful_results]),
            mean([r.total_time for r in successful_results])
        ]
    )
    
    CSV.write(joinpath(results_dir, "time_breakdown.csv"), time_breakdown)
    
    # Write summary to log
    summary_text = "\n\n=== WARMSTART EXPERIMENT SUMMARY ===\n"
    summary_text *= "Total instances processed: $(length(results))\n"
    summary_text *= "Successful: $(length(successful_results))\n"
    summary_text *= "Failed: $(length(results) - length(successful_results))\n\n"
    
    summary_text *= "=== AVERAGE RESULTS ===\n"
    summary_text *= "Heuristic objective: $(round(mean([r.heuristic_obj for r in successful_results]), digits=2))\n"
    summary_text *= "Heuristic time: $(round(mean([r.total_heuristic_time for r in successful_results]), digits=2))s\n"
    summary_text *= "Gurobi objective: $(round(mean([r.gurobi_obj for r in successful_results]), digits=2))\n"
    summary_text *= "Gurobi bound: $(round(mean([r.gurobi_bound for r in successful_results]), digits=2))\n"
    summary_text *= "Gurobi time: $(round(mean([r.gurobi_actual_time for r in successful_results]), digits=2))s\n"
    summary_text *= "Total time: $(round(mean([r.total_time for r in successful_results]), digits=2))s\n"
    summary_text *= "Warmstart improvement: $(round(mean([r.warmstart_improvement for r in successful_results]), digits=2))\n"
    summary_text *= "Warmstart improvement %: $(round(mean([r.warmstart_improvement_pct for r in successful_results]), digits=2))%\n"
    summary_text *= "Final gap: $(round(mean([r.gurobi_gap for r in successful_results]), digits=2))%\n\n"
    
    summary_text *= "=== IMPROVEMENT BREAKDOWN ===\n"
    summary_text *= "Improved by Gurobi: $num_improved ($(round(num_improved/length(successful_results)*100, digits=1))%)\n"
    summary_text *= "Same as heuristic: $num_same ($(round(num_same/length(successful_results)*100, digits=1))%)\n"
    summary_text *= "Worse than heuristic: $num_worse ($(round(num_worse/length(successful_results)*100, digits=1))%)\n\n"
    
    summary_text *= "=== TIME BREAKDOWN ===\n"
    for row in eachrow(time_breakdown)
        summary_text *= "$(row.phase): $(round(row.avg_time, digits=2))s\n"
    end
    
    println(summary_text)
    open(log_file, "a") do f
        write(f, summary_text)
    end
end

# Main function
function main()
    if length(ARGS) == 0
        println("Using default config file: config_warmstart.toml")
        config_file = "config_warmstart.toml"
    else
        config_file = ARGS[1]
    end
    
    run_warmstart_experiment(config_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end