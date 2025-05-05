using Types, JuMP, Gurobi, DelimitedFiles
function build_model(instance::Instance)
    B = instance.B
    S = instance.S
    D = instance.D
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    P = instance.P
    V = instance.V
    μ = instance.μ
    T = instance.T
    R = instance.R
    β = instance.β
    M = instance.M # Number of activities
    K = instance.K # Number of facility types
    
    model = Model()
    @variable(model, x[1:S, 1:B], Bin)
    @variable(model, y[1:S], Bin)
    @objective(model, Min, sum(D .* x))
    
    @constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)
    @constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y[i])
    @constraint(model, cardinality, sum(y) == P)
    @constraint(model, risk[i in 1:S], sum(x[i, j] * R[j] for j in 1:B) <= β[i])
    
    @constraint(
        model,
        tol_l[i in 1:S, m in 1:M],
        sum(x[i, j] * V[m][j] for j in 1:B) >= y[i] * ceil(Int, μ[m][i] * (1 - T[m]))
    )
    @constraint(
        model,
        tol_u[i in 1:S, m in 1:M],
        sum(x[i, j] * V[m][j] for j in 1:B) <= y[i] * floor(Int, μ[m][i] * (1 + T[m]))
    )
    
    @constraint(
        model,
        low_k[k in 1:K],
        Lk[k] <= sum(y[i] for i in Sk[k]),
    )
    @constraint(
        model,
        upp_k[k in 1:K],
        sum(y[i] for i in Sk[k]) <= Uk[k],
    )
    
    return model
end

function optimize_model(model::Model, number, log_file, results_file, time_limit=1800.0, method=1; verbose=true, solver=Gurobi)
    set_optimizer(model, solver.Optimizer)
    set_time_limit_sec(model, time_limit)
    
    if !verbose
        set_silent(model)
    end
    
    set_optimizer_attribute(model, "LogFile", log_file)
    set_optimizer_attribute(model, "Method", method)
    
    optimize!(model)
    
    tiempo = MOI.get(model, MOI.SolveTimeSec())
    time_int = trunc(Int, tiempo)
    
    if primal_status(model) != MOI.FEASIBLE_POINT
        @error "Feasible point not reached"
        obj_value = 0 
        X = [0 0; 0 0]
        Y = [0, 0]
    else
        obj_value = trunc(Int, objective_value(model))
        X = value.(model[:x])
        X = round.(Int, X)
        Y = value.(model[:y])
        Y = round.(Int, Y)
        obj_val = objective_value(model)
        bb_val = dual_objective_value(model)
        gap = relative_gap(model)
        
        writedlm(results_file, [obj_val, bb_val, gap, time_int])
        println("Solved in $time_int seconds")
    end
    
    if termination_status(model) != MOI.OPTIMAL
        @warn "Optimum not reached"
    end
    
    return X, Y, obj_value, time_int
end