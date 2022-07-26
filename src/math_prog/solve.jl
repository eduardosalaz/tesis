include("first_model.jl")
using Types
function solve(path::String)
    instance = read_instance(path)
    model = build_model(instance)
    X, Y, obj_val = optimize_model(model)
    solution = Solution(instance, X, Y, obj_val)
    write_solution(solution, "sol_5_cpx.jld2")
    plot_solution(solution, "plot_5_cpx.png")
end


solve(ARGS[1])