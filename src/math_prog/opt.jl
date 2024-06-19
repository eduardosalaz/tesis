include("first_model.jl")
using Gurobi
using CPLEX
using HiGHS
using DelimitedFiles

instance = read_instance(ARGS[1])
model = build_model(instance)
X, Y, obj_value, time_int = optimize_model(model, 0)
println(time_int)
solution = Solution(instance, X, Y, obj_value, time_int)
plot_solution(solution, "solucion_gurobi_S_1.png")
write_solution(solution, "solucion_gurobi_S_1.jld2")
