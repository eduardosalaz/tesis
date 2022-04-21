include("../generator/generator.jl")
include("first_model.jl")
using .Types

function generate_solve()
    size = "150_80_65"
    K = 5
    M = 3
    B = 150
    S = 80
    P = 65
    failed_dir_path = "failed_instances/"
    inst_dir_path = "instances/experiments/"
    sol_dir_path = "out/solutions/experiments/"
    for i in 1:100
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        parameters = generate_params(B, S)
        instance = Types.Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat, parameters...)
        model = build_model(instance)
        obj_val, X, Y = optimize_model(model)
        file_inst_path = "inst_$i" * size* ".jld2"
        file_sol_path = "sol_$i" * size * ".jld2"
        if obj_val == 0
            @error "Instancia $i no resuelta"
            full_failed_inst_path = failed_dir_path * file_sol_path
            write_instance(instance, full_failed_inst_path)
            # como determinar si la instancia no es la que es infactible? TODO
        end
        solution = Types.Solution(instance, X, Y, obj_val)
        full_sol_path = sol_dir_path * file_sol_path
        full_inst_path = inst_dir_path * file_inst_path
        plot_sol_path = "out/plots/solutions/sol$i" * size *  ".png"
        plot_inst_path = "out/plots/instances/inst$i" * size * ".png"
        write_instance(instance, full_inst_path)
        write_solution(solution, full_sol_path)
        plot_instance(instance, plot_inst_path)
        plot_solution(solution, plot_sol_path)
    end
end

generate_solve()
