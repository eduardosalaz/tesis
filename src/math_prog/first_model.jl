include("../types/types.jl")
using CPLEX, JuMP, JLD2
using .Types

function optimize_model(model::Model; verbose=true, solver=CPLEX::Module)
    set_optimizer(model, solver.Optimizer)
    set_time_limit_sec(model, 600.0) # 10 mins timeout
    if !verbose
        set_silent(model)
    end
    show(model)
    optimize!(model)
    if primal_status(model) != MOI.FEASIBLE_POINT
        @error "Punto factible no alcanzado"
        obj_value = 0 # TODO hacer algo mejor
        X = [0 0 ; 0 0]
        Y = [0,0] # idk
    else
        obj_value = trunc(Int, objective_value(model))
        X = value.(model[:x])
        X = trunc.(Int, X)
        Y = value.(model[:y])
        Y = trunc.(Int, Y)
    end
    if termination_status(model) != MOI.OPTIMAL
        @warn "Óptimo no alcanzado"
    end
    println(obj_value)
    return X, Y, obj_value
end

function build_model(instance::Types.Instance)
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
    k = 5 # types of branches

    model = Model() # THIS IS WHERE THE FUN BEGINS

    @variable(model, x[1:S, 1:B], Bin)
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

    @constraint(model, risk[j in 1:B, i in 1:S], x[i, j] * R[j] <= β[i])

    # ∑ j ∈ B Xᵢⱼ Rⱼ ≤ βᵢ, ∀ i ∈ S

    @constraint(
        model,
        tol_l[i in 1:S, M in 1:m],
        y[i] * μ[M][i] * (1 - T[M]) <= sum(x[i, j] * V[m][j] for j in 1:B),
    )

    @constraint(
        model,
        tol_u[i in 1:S, M in 1:m],
        sum(x[i, j] * V[m][j] for j in 1:B) <= y[i] * μ[M][i] * (1 + T[M]),
    )

    # Yᵢμₘⁱ(1-tᵐ) ≤ ∑i∈S Xᵢⱼ vⱼᵐ ≤ Yᵢμₘʲ(1+tᵐ) ∀ j ∈ B, m = 1 … 3

    @constraint(
        model,
        low_k_types[K in 1:k],
        Lk[K] <= sum(y[i] for i in Sk[K]),
    )

    @constraint(
        model,
        upp_k_types[K in 1:k],
        sum(y[i] for i in Sk[K]) <= Uk[K],
    )

    # lₖ ≤ ∑i ∈ Sₖ Yᵢ ≤ uₖ, k = 1 … 5

    return model
end
