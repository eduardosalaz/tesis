include("first_model.jl")
using Gurobi

instance = read_instance(ARGS[1])
model = build_model(instance)
res = optimize_model(model)