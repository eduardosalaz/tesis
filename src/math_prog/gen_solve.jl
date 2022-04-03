include("../generator/generator.jl")
include("first_model.jl")
using .Types

function generate()

    K = 5
    M = 3
    size = "S"
    B = 25
    S = 10
    P = 8
    for i in 1:10
        #BU_coords, S_coords = generate_coords(B, S)
        #dist_mat = generate_dist(BU_coords, S_coords, B, S)
        #parameters = generate_params(size, B, S)
        #instance = Types.Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat,parameters...)
        dir_path = "instances/"
        file_name = "inst_written$i" * size * ".jld2"
        full_path = dir_path * file_name
        #jldsave(full_path; instance)
            instancia = jldopen(full_path)
            instance = instancia["instance"]
            build_model(instance, i)
    end
end

generate()
