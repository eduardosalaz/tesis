include("generator.jl")
include("solver.jl")
using Dates

function run_experiment(config_file="config.toml")
    # Check if config file exists
    if !isfile(config_file)
        error("Configuration file not found: $config_file")
    end
    
    # Read configuration
    config = read_config(config_file)
    base_path = get(config["general"], "base_path", ".")
    prefix = get(config["general"], "prefix", "005_newk")
    
    # Create a timestamp for this experiment run
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    experiment_id = "exp_$(timestamp)_$(prefix)"
    
    # Create experiment directory
    experiment_dir = joinpath(base_path, "experiments", experiment_id)
    if !isdir(experiment_dir)
        mkpath(experiment_dir)
    end
    
    # Copy config file to experiment directory
    config_copy_path = joinpath(experiment_dir, "config.toml")
    cp(config_file, config_copy_path, force=true)
    
    println("=== Starting Experiment: $experiment_id ===")
    println("Configuration saved to: $config_copy_path")
    
    # Create experiment log file
    log_file = joinpath(experiment_dir, "experiment_log.txt")
    open(log_file, "w") do f
        write(f, "Experiment ID: $experiment_id\n")
        write(f, "Started at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n")
        write(f, "Configuration file: $config_file\n\n")
    end
    
    # Process each instance set
    for (set_name, set_config) in config["instance_sets"]
        B = set_config["B"]
        S = set_config["S"]
        P_values = set_config["P_values"]
        num_instances = get(set_config, "num_instances", 1)
        
        set_info = "\n=== Processing Instance Set: $set_name (B=$B, S=$S) ==="
        println(set_info)
        open(log_file, "a") do f
            write(f, set_info * "\n")
        end
        
        # Create base directories
        base_data_dir = joinpath(experiment_dir, "base_data", "$(B)_$(S)")
        if !isdir(base_data_dir)
            mkpath(base_data_dir)
        end
        
        # For each instance number, generate unique base data
        for i in 1:num_instances
            instance_info = "\n--- Generating base data for instance $i ---"
            println(instance_info)
            open(log_file, "a") do f
                write(f, instance_info * "\n")
            end
            
            # Generate unique base data for this instance
            base_data = generate_base_data(B, S, config)
            BU_coords, S_coords, dist_mat, V, R, Sk, percentages = base_data
            
            # Save base data for this instance
            base_data_file = joinpath(base_data_dir, "base_data_$(B)_$(S)_instance_$i.jld2")
            @save base_data_file BU_coords S_coords dist_mat V R Sk percentages
            
            # Process each P value for this instance
            for P in P_values
                size_str = "$(B)_$(S)_$(P)"
                p_info = "Processing instance $i with P=$P"
                println(p_info)
                open(log_file, "a") do f
                    write(f, p_info * "\n")
                end
                
                # Create directories for this P value
                inst_dir_path = joinpath(experiment_dir, "insts", size_str)
                failed_dir_path = joinpath(experiment_dir, "insts", size_str, "failed_instances")
                sol_dir_path = joinpath(experiment_dir, "solutions", size_str)
                plot_out_dir_path = joinpath(experiment_dir, "plots", size_str)
                plot_sol_dir_path = joinpath(plot_out_dir_path, "solutions")
                plot_inst_dir_path = joinpath(plot_out_dir_path, "instances")
                logs_dir_path = joinpath(experiment_dir, "logs")
                
                for dir in [inst_dir_path, failed_dir_path, sol_dir_path, 
                            plot_out_dir_path, plot_sol_dir_path, plot_inst_dir_path, logs_dir_path]
                    if !isdir(dir)
                        mkpath(dir)
                    end
                end
                
                # Create an instance with this base data and P-specific parameters
                instance = create_instance(B, S, P, config, base_data)
                
                # Save instance
                file_inst_path = "inst_$(i)_$(prefix)_$(size_str).jld2"
                full_inst_path = joinpath(inst_dir_path, file_inst_path)
                write_instance(instance, full_inst_path)
                
                # Plot instance
                plot_inst_path = joinpath(plot_inst_dir_path, "inst_$(i)_$(prefix)_$(size_str).png")
                plot_instance(instance, plot_inst_path)
                
                solve_info = "Building model for instance $i with P=$P..."
                println(solve_info)
                open(log_file, "a") do f
                    write(f, solve_info * "\n")
                end
                
                model = build_model(instance)
                
                # Configure solver
                time_limit = get(config["general"], "time_limit", 1800.0)
                method = get(config["general"], "method", 1)
                solver_log_file = joinpath(logs_dir_path, "log_$(size_str)_$(prefix)_$(time_limit)s_$i.txt")
                results_file = joinpath(logs_dir_path, "results_$(size_str)_$(prefix)_$(time_limit)s_$i.txt")
                
                # Solve model
                println("Solving model for instance $i with P=$P...")
                X, Y, obj_val, solve_time = optimize_model(model, i, solver_log_file, results_file, time_limit, method)
                
                if obj_val == 0
                    error_info = "Instance $i with P=$P not solved"
                    @error error_info
                    open(log_file, "a") do f
                        write(f, "ERROR: $error_info\n")
                    end
                    
                    # Save failed instance
                    full_failed_inst_path = joinpath(failed_dir_path, file_inst_path)
                    write_instance(instance, full_failed_inst_path)
                else
                    # Save solution
                    file_sol_path = "sol_$(i)_$(prefix)_$(size_str).jld2"
                    full_sol_path = joinpath(sol_dir_path, file_sol_path)
                    solution = Solution(instance, X, Y, obj_val, solve_time)
                    write_solution(solution, full_sol_path)
                    
                    # Plot solution
                    plot_sol_path = joinpath(plot_sol_dir_path, "sol_$(i)_$(prefix)_$(size_str).png")
                    plot_solution(solution, plot_sol_path)
                    
                    result_info = "Solved instance $i with P=$P, objective value $obj_val in $solve_time seconds"
                    println(result_info)
                    open(log_file, "a") do f
                        write(f, result_info * "\n")
                    end
                end
            end
        end
    end
    
    # Record experiment completion
    completion_info = "\n=== Experiment Complete at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) ==="
    println(completion_info)
    open(log_file, "a") do f
        write(f, completion_info * "\n")
    end
    
    println("Results stored in: $experiment_dir")
    println("Configuration saved to: $config_copy_path")
end

# Combined main function that supports both interfaces
function main()
    if length(ARGS) == 0
        println("Usage:")
        println("  julia run_experiment.jl config.toml              # Run full experiment from config")
        println("  julia run_experiment.jl size_str number config   # Legacy interface")
        return
    end
    
    if length(ARGS) >= 2 && contains(ARGS[1], "_")
        # Traditional interface (assuming size_str contains underscores)
        size_str = ARGS[1]
        number = ARGS[2]
        config_file = length(ARGS) >= 3 ? ARGS[3] : "config.toml"
        
        println("Using legacy interface with size=$size_str, instances=$number")
        generate_solve(size_str, number, config_file=config_file)
    else
        # Config-based interface
        config_file = ARGS[1]
        run_experiment(config_file)
    end
end

# Execute main function only once
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end