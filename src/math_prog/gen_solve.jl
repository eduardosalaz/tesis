include("../generator/generator.jl")
include("first_model.jl")
#TODO fix this
function generate_solve()
    size = "500_120_50"
    K = 5
    M = 3
    B = 500
    S = 120
    P = 50
    failed_dir_path = "instances/experiments7/failed_instances/"
    inst_dir_path = "instances/experiments7/"
    sol_dir_path = "out/solutions/experiments7/"
    #=
    if !isdir(failed_dir_path)
        mkdir(failed_dir_path)
    end
    if !isdir(inst_dir_path)
        mkdir(inst_dir_path)
    end
    if !isdir(sol_dir_path)
        mkdir(sol_dir_path)
    end
    =#
    for i in 1:100
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        parameters = generate_params(B, S, P)
        instance = Types.Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat, parameters...)
        model = build_model(instance)
        X, Y, obj_val = optimize_model(model)
        file_inst_path = "inst_$i" * "_" * size* ".jld2"
        file_sol_path = "sol_$i" * "_" * size * ".jld2"
        if obj_val == 0
            @error "Instancia $i no resuelta"
            full_failed_inst_path = failed_dir_path * file_inst_path
            Types.write_instance(instance, full_failed_inst_path)
            continue
            # como determinar si la instancia no es la que es infactible? TODO
        end
        solution = Types.Solution(instance, X, Y, obj_val)
        full_sol_path = sol_dir_path * file_sol_path
        full_inst_path = inst_dir_path * file_inst_path
        plot_sol_dir_path = "out/plots/experiments7/solutions/"
        plot_sol_path = plot_sol_dir_path * "sol$i" * "_" * size *  ".png"
        plot_inst_dir_path = "out/plots/experiments7/instances/"
        plot_inst_path = plot_inst_dir_path * "inst$i" * "_" * size * ".png"
        Types.write_instance(instance, full_inst_path)
        Types.write_solution(solution, full_sol_path)
        Types.plot_instance(instance, plot_inst_path)
        Types.plot_solution(solution, plot_sol_path)
    end
end

generate_solve()
