include("constructive.jl")
include("local.jl")


function main(path)
    if !isdir(path)
        init_methods = ["relax", "pdisp", "random"]
        assign_methods = ["naive", "opp"]
        for init_method in init_methods, assign_method in assign_methods 
            solucion = constructive(path, init_method, assign_method)
            solucionLs = main(solucion)
        end
    else
        init_methods = ["relax", "pdisp", "random"]
        assign_methods = ["naive", "opp"]
        paths = readdir(path)
        for path_str in paths
            for init_method in init_methods, assign_method in assign_methods 
                solucion = constructive(path, init_method, assign_method)
                solucionLs = main(solucion)
            end
        end
    end
end

main(ARGS[1])