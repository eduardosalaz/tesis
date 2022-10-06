include("first_model.jl")

instance = read_instance(ARGS[1])
model = build_model(instance)
res = optimize_model(model)