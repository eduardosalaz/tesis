using Types
include("first_model.jl")
#TODO fix this
function solve(size, dir)
    K = 5
    M = 3
    contador = 1
    B, S, P = parse.(Int, split(size, "_"))
    failed_dir_path = "instances/" * size * "/failed_instances/"
    inst_dir_path = "instances/" * size * "/"
    sol_dir_path = "out/solutions/" * size * "/"
    plot_out_dir_path = "out/plots/" * size * "/"
    plot_sol_dir_path = "out/plots/" * size * "/solutions/"
    plot_inst_dir_path = "out/plots/" * size * "/instances/"
    if !isdir(inst_dir_path)
        mkdir(inst_dir_path)
    end

    if !isdir(sol_dir_path)
        mkdir(sol_dir_path)
    end

    if !isdir(plot_out_dir_path)
        mkdir(plot_out_dir_path)
    end

    if !isdir(plot_sol_dir_path)
        mkdir(plot_sol_dir_path)
    end
    if !isdir(plot_inst_dir_path)
        mkdir(plot_inst_dir_path)
    end
    
    entradas = readdir(dir)
    for entrada in entradas
        path = dir * entrada
        pattern = Regex("[t][_]\\d{1,3}")
        index = findfirst(pattern, path)
        almost_number = path[index]
        _, number = split(almost_number, "_")
        println(number)
        instance = read_instance(path)
        model = build_model(instance)

        X, Y, obj_val, time = optimize_model(model, number)
        file_inst_path = "inst_$number" * "_" * size * ".jld2"
        file_sol_path = "sol_$number" * "_" * size * ".jld2"
        if obj_val == 0
            @error "Instancia $number no resuelta"
            full_failed_inst_path = failed_dir_path * file_inst_path
            write_instance(instance, full_failed_inst_path)
            # como determinar si la instancia no es la que es infactible? TODO
        end
        solution = Solution(instance, X, Y, obj_val, time)
        full_sol_path = sol_dir_path * file_sol_path
        full_inst_path = inst_dir_path * file_inst_path
        plot_sol_path = plot_sol_dir_path * "sol$number" * "_" * size *  ".png"
        plot_inst_path = plot_inst_dir_path * "inst$number" * "_" * size * ".png"

        write_solution(solution, full_sol_path)
        #plot_solution(solution, plot_sol_path)
    end
end

solve(ARGS[1], ARGS[2])