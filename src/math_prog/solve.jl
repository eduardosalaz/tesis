include("first_model.jl")
using Types
function solve(path::String)
    instance = read_instance(path)
    model = build_model(instance)
    X, Y, obj_val = optimize_model(model)
    println(obj_val)
    solution = Solution(instance, X, Y, obj_val)
    write_solution(solution, "1solucionperra1.jld2")
    plot_solution(solution, "1plotperro1.png")
end


solve(ARGS[1])