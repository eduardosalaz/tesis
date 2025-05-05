using JuMP, Gurobi, JLD2, DelimitedFiles, TOML
using Types
include("../generator/generator.jl")

function read_config(config_file="config.toml")
    isfile(config_file) || error("Configuration file not found: $config_file")
    return TOML.parsefile(config_file)
end

function optimize_model(model::Model, number, log_file, results_file, time_limit=1800.0, method=1; verbose=true, solver=Gurobi)
    set_optimizer(model, solver.Optimizer)
    set_time_limit_sec(model, time_limit) # Default: 30 minutes
    
    if !verbose
        set_silent(model)
    end
    
    set_optimizer_attribute(model, "LogFile", log_file)
    set_optimizer_attribute(model, "Method", method)
    
    optimize!(model)
    
    tiempo = MOI.get(model, MOI.SolveTimeSec())
    time_int = trunc(Int, tiempo)
    
    if primal_status(model) != MOI.FEASIBLE_POINT
        @error "Feasible point not reached"
        obj_value = 0 
        X = [0 0; 0 0]
        Y = [0, 0]
    else
        obj_value = trunc(Int, objective_value(model))
        X = value.(model[:x])
        X = round.(Int, X)
        Y = value.(model[:y])
        Y = round.(Int, Y)
        obj_val = objective_value(model)
        bb_val = dual_objective_value(model)
        gap = relative_gap(model)
        
        writedlm(results_file, [obj_val, bb_val, gap, time_int])
        println("Solved in $time_int seconds")
    end
    
    if termination_status(model) != MOI.OPTIMAL
        @warn "Optimum not reached"
    end
    
    return X, Y, obj_value, time_int
end

function build_model(instance::Instance)
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
    
    model = Model() # THIS IS WHERE THE FUN BEGINS
    @variable(model, x[1:S, 1:B], Bin)
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

# Function to solve instances from configuration file
function solve_from_config(config_file="config.toml")
    config = read_config(config_file)
    base_path = get(config["general"], "base_path", ".")
    prefix = get(config["general"], "prefix", "005_newk")
    time_limit = get(config["general"], "time_limit", 1800.0)
    method = get(config["general"], "method", 1)
    
    for (set_name, set_config) in config["instance_sets"]
        B = set_config["B"]
        S = set_config["S"]
        
        for P in set_config["P_values"]
            size_str = "$(B)_$(S)_$(P)"
            println("Processing: B=$B, S=$S, P=$P")
            
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
            
            # Find existing instances to solve
            instance_files = filter(file -> startswith(file, "inst_") && endswith(file, ".jld2"), 
                                  readdir(inst_dir_path))
            
            if isempty(instance_files)
                println("No instances found in $inst_dir_path, generating one...")
                # Generate base data (independent of P)
                BU_coords, S_coords = generate_coords(B, S)
                dist_mat = generate_dist(BU_coords, S_coords, B, S)
                
                # Generate parameters (dependent on P)
                parameters = generate_params(B, S, P)
                
                # Create instance
                instance = Instance(B, S, 4, 3, P, BU_coords, S_coords, dist_mat, parameters...)
                
                # Save instance
                file_name = "inst_1_$(prefix)_$(size_str).jld2"
                full_path = joinpath(inst_dir_path, file_name)
                write_instance(instance, full_path)
                
                instance_files = [file_name]
            end
            
            # Process each instance
            for inst_file in instance_files
                i = parse(Int, split(inst_file, "_")[2])  # Extract instance number
                
                # Load instance
                full_inst_path = joinpath(inst_dir_path, inst_file)
                instance = read_instance(full_inst_path)
                
                # Check if solution already exists
                sol_file = replace(inst_file, "inst_" => "sol_")
                full_sol_path = joinpath(sol_dir_path, sol_file)
                
                if isfile(full_sol_path)
                    println("Solution for instance $i already exists, skipping...")
                    continue
                end
                
                println("Building model for instance $i...")
                model = build_model(instance)
                
                # Configure log files
                log_file = joinpath(logs_dir_path, "log_$(size_str)_$(prefix)_$(time_limit)s_$i.txt")
                results_file = joinpath(logs_dir_path, "results_$(size_str)_$(prefix)_$(time_limit)s_$i.txt")
                
                # Solve model
                X, Y, obj_val, solve_time = optimize_model(model, i, log_file, results_file, time_limit, method)
                
                if obj_val == 0
                    @error "Instance $i not solved"
                    # Save failed instance
                    failed_file = joinpath(failed_dir_path, inst_file)
                    write_instance(instance, failed_file)
                else
                    # Save solution
                    solution = Solution(instance, X, Y, obj_val, solve_time)
                    write_solution(solution, full_sol_path)
                    
                    # Plot solution
                    plot_sol_path = joinpath(plot_sol_dir_path, replace(sol_file, ".jld2" => ".png"))
                    plot_inst_path = joinpath(plot_inst_dir_path, replace(inst_file, ".jld2" => ".png"))
                    
                    plot_instance(instance, plot_inst_path)
                    plot_solution(solution, plot_sol_path)
                    
                    println("Solved instance $i with objective value $obj_val")
                end
            end
        end
    end
    println("All solving complete.")
end

# Command-line interface
function main()
    if length(ARGS) == 0
        println("Usage:")
        println("  julia solver.jl config.toml             # Solve all instances from config")
        println("  julia solver.jl size_str number [config.toml] # Solve specific instances")
        println("Example: julia solver.jl 625_200_30 5 config.toml")
        return
    end
    
    if length(ARGS) >= 2
        # Traditional interface
        size_str = ARGS[1]
        number = ARGS[2]
        config_file = length(ARGS) >= 3 ? ARGS[3] : "config.toml"
        generate_solve(size_str, number, config_file=config_file)  # Use the function from generator.jl
    else
        # Config-based interface
        config_file = ARGS[1]
        solve_from_config(config_file)
    end
end

# Uncomment to run when script is executed directly
# main()