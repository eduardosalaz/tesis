include("../generator/generator.jl")
include("first_model.jl")
#TODO fix this
function generate_solve(size, number)
    K = 4
    M = 3
    contador = 1
    B, S, P = parse.(Int, split(size, "_"))
    println(B, S, P)
    failed_dir_path = "insts_new/" * size * "/failed_instances/"
    inst_dir_path = "insts_new/" * size * "/"
    sol_dir_path = "out_new/solutions/" * size * "/"
    plot_out_dir_path = "out_new/plots/" * size * "/"
    plot_sol_dir_path = "out_new/plots/" * size * "/solutions/"
    plot_inst_dir_path = "out_new/plots/" * size * "/instances/"
    if !isdir(inst_dir_path)
        mkdir(inst_dir_path)
    end
    if !isdir(failed_dir_path)
        mkdir(failed_dir_path)
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
    i = number
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        parameters = generate_params(B, S, P)
        instance = Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat, parameters...)
        println("builded instance $i")
        model = build_model(instance)
        println("builded model $i")
        X, Y, obj_val, time = optimize_model(model, i)
        file_inst_path = "inst_$i" * "_010_newk" * "_" * size * ".jld2"
        file_sol_path = "sol_$i" * "_010_newk" * "_" * size * ".jld2"
        if obj_val == 0
            @error "Instancia $i no resuelta"
            full_failed_inst_path = failed_dir_path * file_inst_path
            write_instance(instance, full_failed_inst_path)
            # como determinar si la instancia no es la que es infactible? con el IIS
        else
            solution = Solution(instance, X, Y, obj_val, time)
            full_sol_path = sol_dir_path * file_sol_path
            full_inst_path = inst_dir_path * file_inst_path
            plot_sol_path = plot_sol_dir_path * "sol_$(i)" * "_010_newk" * "_" * size * ".png"
            plot_inst_path = plot_inst_dir_path * "inst_$(i)" * "_010_newk" * "_" * size * ".png"
            contador += 1
            write_instance(instance, full_inst_path)
            write_solution(solution, full_sol_path)
            plot_instance(instance, plot_inst_path)
            plot_solution(solution, plot_sol_path)
        end
    
end


generate_solve(ARGS[1], ARGS[2])
