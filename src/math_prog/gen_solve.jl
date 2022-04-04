include("../generator/generator.jl")
include("first_model.jl")
using .Types

function generate()

    K = 5
    M = 3
    size = "S"
    B = 150
    S = 80
    P = 65
    for i in 1:100
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = generate_dist(BU_coords, S_coords, B, S)
        parameters = generate_params(size, B, S)
        instance = Types.Instance(B, S, K, M, P, BU_coords, S_coords, dist_mat,parameters...)
        try
            build_model(instance, i)
            dir_path = "instL_100/"
            file_name = "inst_new_$i" * size * ".jld2"
            full_path = dir_path * file_name
            jldsave(full_path; instance)
        catch
            @error "En la $i"
        end

    end
end

generate()
