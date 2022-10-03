include("constructive.jl")
include("local.jl")


function main_cons(path)
    words = split(path, "_")
    id = words[4]
    if !isdir(path)
        init_methods = ["relax", "pdisp", "random"]
        assign_methods = ["naive", "opp"]
        for init_method in init_methods, assign_method in assign_methods 
            solucion = constructive(path, id, init_method, assign_method)
            solucionLs = main(solucion)
        end
    end
end

main_cons(ARGS[1])