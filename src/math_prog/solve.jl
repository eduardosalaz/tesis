include("first_model.jl")
using Types
function solve(path::String)
    instance = read_instance(path)
    model = build_model(instance)
    X, Y, obj_val = optimize_model(model)
    solution = Solution(instance, X, Y, obj_val)
    write_solution(solution, "sol_nuevo.jld2")
    plot_solution(solution, "plot_nuevo.png")
end


solve(ARGS[1])