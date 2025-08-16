# run_local_search_ablation.jl
# Experiment 2: Test local search move combinations with the best Y initialization method

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
struct AblationResult
    instance_id::String
    instance_number::Int
    y_method::String  # Which Y initialization was used
    move_config::String
    B::Int
    S::Int
    P::Int
    # Construction timing (same for all move configs)
    init_time::Float64
    phase2_time::Float64
    repair_time::Float64
    construction_time::Float64
    # Local search timing
    ls_time::Float64
    total_time::Float64
    # Objective values
    phase2_obj::Float64
    split_resolved_obj::Float64
    post_repair_obj::Float64
    post_ls_obj::Float64
    # Local search details
    ls_iterations::Int
    total_ls_improvement::Float64
    successful_moves_simple::Int
    successful_moves_interchange::Int
    successful_moves_deactivate::Int
    # Status
    status::String
    error_message::String
end

# Structure to hold pre-LS solution data
struct PreLSSolution
    instance_name::String
    instance_number::Int
    solution::Solution
    Y::Vector{Int}
    targets_lower::Any
    targets_upper::Any
    construction_stats::NamedTuple
end

# Define move configurations for ablation study
function get_ablation_move_configs()
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

# Build pre-LS solution for a single instance
function build_pre_ls_solution(instance, instance_file, instance_number, y_method, config, output_dirs, experiment_log)
    instance_name = splitext(basename(instance_file))[1]
    B, S, P = instance.B, instance.S, instance.P
    
    println("  Building pre-LS solution for instance #$instance_number")
    
    total_start_time = time()
    
    try
        # Phase 1: Initialize Y using the specified method
        y_log_file = joinpath(output_dirs["logs"], "y_init_$(instance_name).log")
        
        if y_method == "relaxed"
            Y, init_time, success = init_y_relaxed(instance, config["time_limits"]["exact_methods"], y_log_file)
        elseif y_method == "pdisp"
            Y, init_time, success = init_y_pdisp(instance, config["time_limits"]["exact_methods"])
        elseif y_method == "multi_pdp"
            Y, init_time, success = init_y_multi_pdp(instance, config["time_limits"]["exact_methods"], y_log_file)
        elseif y_method == "multi_p_median"
            Y, init_time, success = init_y_multi_p_median(instance, config["time_limits"]["exact_methods"], y_log_file, false)
        elseif y_method == "multi_p_median_benders"
            Y, init_time, success = init_y_multi_p_median(instance, config["time_limits"]["exact_methods"], y_log_file, true)
        else
            error("Unknown method: $y_method")
        end
        
        if !success || Y === nothing
            println("    Y initialization failed")
            return nothing
        end
        
        # Phase 2: Solve transportation problem
        phase2_start = time()
        phase2_log_file = joinpath(output_dirs["logs"], "phase2_$(instance_name).log")
        
        model_transport = solve_phase2_model(instance, Y)
        set_optimizer(model_transport, Gurobi.Optimizer)
        set_optimizer_attribute(model_transport, "LogFile", phase2_log_file)
        set_time_limit_sec(model_transport, config["time_limits"]["phase2"])
        
        optimize!(model_transport)
        
        if termination_status(model_transport) != MOI.OPTIMAL && termination_status(model_transport) != MOI.TIME_LIMIT
            println("    Phase 2 optimization failed")
            return nothing
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
        
        # Create solution object
        solution = Solution(instance, x_repaired, Y, post_repair_obj, 1)
        
        # Calculate targets for local search
        targets_lower, targets_upper = calculate_targets_optimized(instance)
        
        construction_time = time() - total_start_time
        
        # Create stats tuple
        stats = (
            init_time = init_time,
            phase2_time = phase2_time,
            repair_time = repair_time,
            construction_time = construction_time,
            phase2_obj = phase2_obj,
            split_resolved_obj = split_resolved_obj,
            post_repair_obj = post_repair_obj,
            num_split_bus = num_splits,
            prop_split_bus = prop_splits,
            activity_violations_initial = initial_activity_viol,
            risk_violations_initial = initial_risk_viol,
            remaining_violations = final_activity_viol + final_risk_viol
        )
        
        println("    Pre-LS solution built: obj = $(round(post_repair_obj, digits=2)), time = $(round(construction_time, digits=2))s")
        
        return PreLSSolution(instance_name, instance_number, solution, Y, targets_lower, targets_upper, stats)
        
    catch e
        error_msg = sprint(showerror, e)
        println("    ERROR in construction: $error_msg")
        return nothing
    end
end

# Main experiment runner
function run_local_search_ablation(config_file="config_ls_ablation.toml")
    # Load configuration
    config = TOML.parsefile(config_file)
    
    # Get the Y initialization method to use
    y_method = config["method"]["y_initialization"]
    println("Using Y initialization method: $y_method")
    
    # Create experiment directory
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    experiment_name = get(config["general"], "experiment_name", "local_search_ablation")
    experiment_id = "exp_$(timestamp)_$(experiment_name)_$(y_method)"
    experiment_dir = joinpath(config["general"]["base_path"], "experiments", experiment_id)
    
    # Create directory structure
    output_dirs = Dict(
        "experiment" => experiment_dir,
        "results" => joinpath(experiment_dir, "results"),
        "instances" => joinpath(experiment_dir, "instances"),
        "pre_ls_solutions" => joinpath(experiment_dir, "pre_ls_solutions"),
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
        Local Search Ablation Experiment
        Y Initialization Method: $y_method
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
    
    # Get move configurations
    move_configs = get_ablation_move_configs()
    configs_to_test = ["all_moves", "simple_interchange", "simple_deactivate", "interchange_deactivate"]
    
    println("Testing move configurations: ", join(configs_to_test, ", "))
    
    # Results storage
    all_results = AblationResult[]
    
    # PHASE 1: Build all pre-LS solutions
    println("\n" * "="^70)
    println("PHASE 1: BUILDING PRE-LS SOLUTIONS WITH $y_method")
    println("="^70)
    
    pre_ls_solutions = Dict{String, PreLSSolution}()  # instance_name => solution
    
    for (idx, instance_file) in enumerate(instance_files)
        instance_name = splitext(basename(instance_file))[1]
        
        # Extract instance number
        match_result = match(r"inst_(\d+)_", instance_name)
        inst_num = match_result !== nothing ? parse(Int, match_result[1]) : idx
        
        println("\n=== Building pre-LS for instance $idx/$(length(instance_files)): #$inst_num - $instance_name ===")
        
        # Load instance
        instance = read_instance(instance_file)
        
        # Copy instance file
        size_str = "$(instance.B)_$(instance.S)_$(instance.P)"
        inst_dir = joinpath(output_dirs["instances"], size_str)
        mkpath(inst_dir)
        cp(instance_file, joinpath(inst_dir, basename(instance_file)))
        
        # Build pre-LS solution
        pre_ls = build_pre_ls_solution(instance, instance_file, inst_num, y_method, config, output_dirs, experiment_log)
        
        if pre_ls !== nothing
            # Save pre-LS solution
            pre_ls_dir = joinpath(output_dirs["pre_ls_solutions"], size_str)
            mkpath(pre_ls_dir)
            pre_ls_file = joinpath(pre_ls_dir, "pre_ls_$(instance_name).jld2")
            
            JLD2.save(pre_ls_file,
                "solution", pre_ls.solution,
                "Y", pre_ls.Y,
                "targets_lower", pre_ls.targets_lower,
                "targets_upper", pre_ls.targets_upper,
                "stats", pre_ls.construction_stats
            )
            
            pre_ls_solutions[instance_name] = pre_ls
        else
            println("  Failed to build pre-LS solution")
        end
    end
    
    println("\nSuccessfully built $(length(pre_ls_solutions))/$(length(instance_files)) pre-LS solutions")
    
    # PHASE 2: Apply different move configurations
    println("\n" * "="^70)
    println("PHASE 2: TESTING LOCAL SEARCH MOVE CONFIGURATIONS")
    println("="^70)
    
    for (instance_name, pre_ls) in pre_ls_solutions
        println("\n=== Testing LS moves for instance #$(pre_ls.instance_number): $instance_name ===")
        println("  Pre-LS objective: $(round(pre_ls.construction_stats.post_repair_obj, digits=2))")
        
        # Find instance file to get dimensions
        instance_file = ""
        for file in instance_files
            if splitext(basename(file))[1] == instance_name
                instance_file = file
                break
            end
        end
        
        instance = read_instance(instance_file)
        size_str = "$(instance.B)_$(instance.S)_$(instance.P)"
        
        for config_name in configs_to_test
            println("  Testing: $config_name")
            
            # Create a copy of the solution
            sol_copy = deepcopy(pre_ls.solution)
            
            # Run local search
            ls_start = time()
            sol_after_ls, ls_time, ls_iters, total_improvement, succ_simple, succ_interchange, succ_deactivate = 
                local_search_with_tracking(sol_copy, pre_ls.targets_lower, pre_ls.targets_upper, move_configs[config_name])
            
            post_ls_obj = sol_after_ls.Weight
            
            # Save solution
            sol_dir = joinpath(output_dirs["solutions"], size_str, config_name)
            mkpath(sol_dir)
            sol_file = joinpath(sol_dir, "sol_$(instance_name)_$(config_name).jld2")
            write_solution(sol_after_ls, sol_file)
            
            # Create result
            result = AblationResult(
                instance_name,
                pre_ls.instance_number,
                y_method,
                config_name,
                instance.B, instance.S, instance.P,
                pre_ls.construction_stats.init_time,
                pre_ls.construction_stats.phase2_time,
                pre_ls.construction_stats.repair_time,
                pre_ls.construction_stats.construction_time,
                ls_time,
                pre_ls.construction_stats.construction_time + ls_time,
                pre_ls.construction_stats.phase2_obj,
                pre_ls.construction_stats.split_resolved_obj,
                pre_ls.construction_stats.post_repair_obj,
                post_ls_obj,
                ls_iters,
                total_improvement,
                succ_simple,
                succ_interchange,
                succ_deactivate,
                "success",
                ""
            )
            
            push!(all_results, result)
            
            println("    Post-LS obj: $(round(post_ls_obj, digits=2)) (improvement: $(round(total_improvement, digits=2)))")
            println("    Moves: S=$succ_simple, I=$succ_interchange, D=$succ_deactivate")
        end
        
        # Save intermediate results
        save_ablation_results(all_results, output_dirs["results"])
    end
    
    # Generate final summary
    generate_ablation_summary(all_results, output_dirs["results"], experiment_log)
    
    # Log completion
    completion_msg = "\n=== Experiment Complete at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) ==="
    println(completion_msg)
    open(experiment_log, "a") do f
        write(f, completion_msg * "\n")
    end
    
    println("\nResults saved to: $experiment_dir")
end

# Save results
function save_ablation_results(results::Vector{AblationResult}, results_dir::String)
    # Convert to DataFrame
    df = DataFrame(
        instance_id = [r.instance_id for r in results],
        instance_number = [r.instance_number for r in results],
        y_method = [r.y_method for r in results],
        move_config = [r.move_config for r in results],
        B = [r.B for r in results],
        S = [r.S for r in results],
        P = [r.P for r in results],
        construction_time = [r.construction_time for r in results],
        ls_time = [r.ls_time for r in results],
        total_time = [r.total_time for r in results],
        post_repair_obj = [r.post_repair_obj for r in results],
        post_ls_obj = [r.post_ls_obj for r in results],
        ls_iterations = [r.ls_iterations for r in results],
        total_ls_improvement = [r.total_ls_improvement for r in results],
        successful_moves_simple = [r.successful_moves_simple for r in results],
        successful_moves_interchange = [r.successful_moves_interchange for r in results],
        successful_moves_deactivate = [r.successful_moves_deactivate for r in results],
        status = [r.status for r in results]
    )
    
    CSV.write(joinpath(results_dir, "ls_ablation_results.csv"), df)
    
    # Save as JSON
    open(joinpath(results_dir, "ls_ablation_results.json"), "w") do f
        JSON.print(f, results, 2)
    end
end

# Generate summary
function generate_ablation_summary(results::Vector{AblationResult}, results_dir::String, log_file::String)
    successful_results = filter(r -> r.status == "success", results)
    
    if isempty(successful_results)
        println("Warning: No successful results to summarize")
        return
    end
    
    # Group by move configuration
    move_configs = unique([r.move_config for r in successful_results])
    
    summary_df = DataFrame()
    
    # Get baseline (all_moves) performance
    baseline_results = filter(r -> r.move_config == "all_moves", successful_results)
    baseline_obj = mean([r.post_ls_obj for r in baseline_results])
    baseline_time = mean([r.ls_time for r in baseline_results])
    
    for config in move_configs
        config_results = filter(r -> r.move_config == config, successful_results)
        
        if !isempty(config_results)
            avg_obj = mean([r.post_ls_obj for r in config_results])
            avg_time = mean([r.ls_time for r in config_results])
            
            summary_row = DataFrame(
                move_config = config,
                num_instances = length(config_results),
                avg_ls_time = avg_time,
                avg_post_ls_obj = avg_obj,
                avg_ls_improvement = mean([r.total_ls_improvement for r in config_results]),
                avg_ls_iterations = mean([r.ls_iterations for r in config_results]),
                avg_simple_moves = mean([r.successful_moves_simple for r in config_results]),
                avg_interchange_moves = mean([r.successful_moves_interchange for r in config_results]),
                avg_deactivate_moves = mean([r.successful_moves_deactivate for r in config_results]),
                obj_degradation = avg_obj - baseline_obj,
                obj_degradation_pct = (avg_obj - baseline_obj) / baseline_obj * 100,
                time_saved = baseline_time - avg_time,
                time_saved_pct = (baseline_time - avg_time) / baseline_time * 100
            )
            summary_df = vcat(summary_df, summary_row)
        end
    end
    
    CSV.write(joinpath(results_dir, "ls_ablation_summary.csv"), summary_df)
    
    # Move contribution analysis
    move_contribution = DataFrame(
        move_disabled = String[],
        avg_obj_increase = Float64[],
        avg_obj_increase_pct = Float64[],
        avg_time_saved = Float64[],
        avg_time_saved_pct = Float64[]
    )
    
    # Analyze impact of disabling each move
    if "simple_interchange" in move_configs  # Deactivate disabled
        config_results = filter(r -> r.move_config == "simple_interchange", successful_results)
        avg_obj = mean([r.post_ls_obj for r in config_results])
        avg_time = mean([r.ls_time for r in config_results])
        push!(move_contribution, ("deactivate", avg_obj - baseline_obj, 
              (avg_obj - baseline_obj) / baseline_obj * 100,
              baseline_time - avg_time, (baseline_time - avg_time) / baseline_time * 100))
    end
    
    if "simple_deactivate" in move_configs  # Interchange disabled
        config_results = filter(r -> r.move_config == "simple_deactivate", successful_results)
        avg_obj = mean([r.post_ls_obj for r in config_results])
        avg_time = mean([r.ls_time for r in config_results])
        push!(move_contribution, ("interchange", avg_obj - baseline_obj,
              (avg_obj - baseline_obj) / baseline_obj * 100,
              baseline_time - avg_time, (baseline_time - avg_time) / baseline_time * 100))
    end
    
    if "interchange_deactivate" in move_configs  # Simple disabled
        config_results = filter(r -> r.move_config == "interchange_deactivate", successful_results)
        avg_obj = mean([r.post_ls_obj for r in config_results])
        avg_time = mean([r.ls_time for r in config_results])
        push!(move_contribution, ("simple", avg_obj - baseline_obj,
              (avg_obj - baseline_obj) / baseline_obj * 100,
              baseline_time - avg_time, (baseline_time - avg_time) / baseline_time * 100))
    end
    
    sort!(move_contribution, :avg_obj_increase_pct, rev=true)  # Most important first
    CSV.write(joinpath(results_dir, "move_contribution.csv"), move_contribution)
    
    # Write summary to log
    summary_text = "\n\n=== LOCAL SEARCH ABLATION SUMMARY ===\n"
    summary_text *= "Y initialization method: $(results[1].y_method)\n\n"
    
    for row in eachrow(summary_df)
        summary_text *= "\nMove config: $(row.move_config)\n"
        summary_text *= "  Avg LS time: $(round(row.avg_ls_time, digits=2))s\n"
        summary_text *= "  Avg final obj: $(round(row.avg_post_ls_obj, digits=2))\n"
        summary_text *= "  Avg improvement: $(round(row.avg_ls_improvement, digits=2))\n"
        if row.move_config != "all_moves"
            summary_text *= "  Degradation vs baseline: $(round(row.obj_degradation_pct, digits=2))%\n"
            summary_text *= "  Time saved vs baseline: $(round(row.time_saved_pct, digits=1))%\n"
        end
    end
    
    summary_text *= "\n\n=== MOVE IMPORTANCE RANKING ===\n"
    for (i, row) in enumerate(eachrow(move_contribution))
        summary_text *= "$i. $(row.move_disabled): "
        summary_text *= "$(round(row.avg_obj_increase_pct, digits=2))% degradation when disabled\n"
    end
    
    println(summary_text)
    open(log_file, "a") do f
        write(f, summary_text)
    end
end

# Main function
function main()
    if length(ARGS) == 0
        println("Using default config file: config_ls_ablation.toml")
        config_file = "config_ls_ablation.toml"
    else
        config_file = ARGS[1]
    end
    
    run_local_search_ablation(config_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end