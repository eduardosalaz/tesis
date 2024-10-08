using Gurobi, JuMP, JLD2
using Types
using MathOptInterface
const MOI = MathOptInterface
using DelimitedFiles

function optimize_model(model::Model, number; verbose=true, solver=Gurobi::Module)
    set_optimizer(model, solver.Optimizer)
    set_time_limit_sec(model, 7200.0) # 30 minutos
    if !verbose
        set_silent(model)
    end
    # show(model)
    set_optimizer_attribute(model, "LogFile", "log_2hs_2500_$number.txt")
    optimize!(model)
    tiempo = MOI.get(model, MOI.SolveTimeSec())
    time_int = trunc(Int, tiempo)
    if primal_status(model) != MOI.FEASIBLE_POINT
        @error "Punto factible no alcanzado"
        obj_value = 0 # TODO hacer algo mejor
        X = [0 0; 0 0]
        Y = [0, 0] # idk
    else
        obj_value = trunc(Int, objective_value(model))
        X = value.(model[:x])
        X = round.(Int, X)
        Y = value.(model[:y])
        Y = round.(Int, Y)
    end
    if termination_status(model) != MOI.OPTIMAL
        @warn "Óptimo no alcanzado"
    end
    obj_val = objective_value(model)
    bb_val = dual_objective_value(model)
    writedlm("obj_bb_$number.txt", [obj_val, bb_val])
    println(time_int)
    return X, Y, obj_value, time_int
end

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

    m = 3 # activities
    k = 5 #  of branches

    model = Model() # THIS IS WHERE THE FUN BEGINS

    @variable(model, x[1:S, 1:B], Bin)
    #@variable(model, x[1:S, 1:B], lower_bound = 0, upper_bound = 1)
    # num suc and num bu, Xᵢⱼ

    @variable(model, y[1:S], Bin)
    # Yᵢ

    @objective(model, Min, sum(D .* x))
    # Xᵢⱼ * Dᵢⱼ

    @constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)

    # ∑ᵢ∈S Xᵢⱼ = 1, ∀ j ∈ B

    @constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y[i])

    # Xᵢⱼ ≤ Yᵢ , ∀ i ∈ S, j ∈ B

    @constraint(model, cardinality, sum(y) == P)

    # ∑ i ∈ S Yᵢ = p

    @constraint(model, risk[i in 1:S], sum(x[i, j] * R[j] for j in 1:B) <= β[i])

    # ∑ j ∈ B Xᵢⱼ Rⱼ ≤ βᵢ, ∀ i ∈ S

    @constraint(
        model,
        tol_l[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) >= (y[i] * μ[M][i] * (1 - T[M])),
    )

    @constraint(
        model,
        tol_u[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) <= (y[i] * μ[M][i] * (1 + T[M])),
    )

    # Yᵢμₘⁱ(1-tᵐ) ≤ ∑i∈S Xᵢⱼ vⱼᵐ ≤ Yᵢμₘʲ(1+tᵐ) ∀ j ∈ B, m = 1 … 3

    @constraint(
        model,
        low_k[K in 1:k],
        Lk[K] <= sum(y[i] for i in Sk[K]),
    )

    @constraint(
        model,
        upp_k[K in 1:k],
        sum(y[i] for i in Sk[K]) <= Uk[K],
    )

    # lₖ ≤ ∑i ∈ Sₖ Yᵢ ≤ uₖ, k = 1 … 5
    return model
end