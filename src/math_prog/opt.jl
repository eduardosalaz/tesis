include("first_model.jl")
using CPLEX
using DelimitedFiles

instance = read_instance(ARGS[1])
model = build_model(instance)
X, Y, obj_value, time_int = optimize_model(model; solver = CPLEX::Module)
solution = Solution(instance, X, Y, obj_value, time_int)
plot_solution(solution, "solucion_cplex2.png")
write_solution(solution, "solucion_cplex2.jld2")
