using Gurobi, JuMP, JLD2
using Types
using MathOptInterface
const MOI = MathOptInterface
using DelimitedFiles
using Plots
using CPLEX
include("../heuristics/constructive.jl")
include("first_model.jl")

function solve_phase1_model(instance::Instance)
    B = instance.B
    S = instance.S
    D = instance.D
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    P = instance.P
    V = instance.V
    μ = instance.μ
    R = instance.R
    β = instance.β
    m = 3 # activities
    k = 4 # branches
    
    model = Model()
    
    # X is continuous, Y is binary
    @variable(model, 0 <= x[1:S, 1:B] <= 1)
    @variable(model, y[1:S], Bin)
    
    @objective(model, Min, sum(D .* x))
    
    @constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)
    @constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y[i])
    @constraint(model, cardinality, sum(y) == P)
    @constraint(model, risk[i in 1:S], sum(x[i, j] * R[j] for j in 1:B) <= β[i])
    
    # Zero tolerance activity constraints
    @constraint(
        model,
        tol[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) == y[i] * μ[M][i]
    )
    
    @constraint(
        model,
        low_k[K in 1:k],
        Lk[K] <= sum(y[i] for i in Sk[K])
    )
    @constraint(
        model,
        upp_k[K in 1:k],
        sum(y[i] for i in Sk[K]) <= Uk[K]
    )
    
    return model
end

# Phase 2: Fix Y values and solve transportation problem
function solve_phase2_model(instance::Instance, y_fixed)
    B = instance.B
    S = instance.S
    D = instance.D
    V = instance.V
    μ = instance.μ
    R = instance.R
    β = instance.β
    m = 3 # activities
    ε = 1e-3
    
    model = Model()
    
    @variable(model, 0 <= x[1:S, 1:B] <= 1)
    
    @objective(model, Min, sum(D .* x))
    
    @constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)
    #@constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y_fixed[i])
    #@constraint(model, risk[i in 1:S], sum(x[i, j] * R[j] for j in 1:B) <= β[i])

    #=
    @constraint(
        model,
        tol[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) == y_fixed[i] * μ[M][i]
    )
        =#
    
    
    @constraint(
        model,
        tol_lb[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) >= (1-ε) * y_fixed[i] * μ[M][i]
    )
    @constraint(
        model,
        tol_ub[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) <= (1+ε) * y_fixed[i] * μ[M][i]
    )
    
    
    return model
end

function solve_phase3_model(instance::Instance, y_fixed, x_fixed)
    B = instance.B
    S = instance.S
    D = instance.D
    V = instance.V
    μ = instance.μ
    R = instance.R
    β = instance.β
    m = 3 # activities
    ε = 1e-3
    
    model = Model()
    
    @variable(model, 0 <= x[1:S, 1:B] <= 1)
    set_start_value.(model[:x], x_fixed)
    
    
    @objective(model, Min, sum(D .* x))
    
    @constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)
    #@constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y_fixed[i])
    @constraint(model, risk[i in 1:S], sum(x[i, j] * R[j] for j in 1:B) <= β[i])

    #=
    @constraint(
        model,
        tol[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) == y_fixed[i] * μ[M][i]
    )
        =#
    
    
    @constraint(
        model,
        tol_lb[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) >= (1-ε) * y_fixed[i] * μ[M][i]
    )
    @constraint(
        model,
        tol_ub[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) <= (1+ε) * y_fixed[i] * μ[M][i]
    )
    
    
    return model
end

function analyze_splits(x_sol::Matrix{Float64}, epsilon::Float64=1e-6)
    splits = 0
    S, B = size(x_sol)
    split_info = Dict()
    
    for j in 1:B
        assignments = [(i, x_sol[i,j]) for i in 1:S if abs(x_sol[i,j]) > epsilon]
        if length(assignments) > 1
            splits += 1
            println("\nBusiness unit $j is split among $(length(assignments)) successors:")
            for (i, val) in sort(assignments, by=x->x[2], rev=true)
                println("  Successor $i: $(round(val, digits=4))")
            end
            split_info[j] = assignments
        end
    end
    
    println("\nTotal number of split business units: $splits")
    return splits, split_info
end

function build_modified_model(Y, instance::Instance)
    B = instance.B
    S = instance.S
    D = instance.D
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    P = instance.P
    V = instance.V
    μ = instance.μ
    R = instance.R
    β = instance.β
    T = instance.T
    m = 3 # activities
    k = 4 # branches
    
    model = Model()
    
    # Relaxed x variable - now continuous between 0 and 1
    @variable(model, x[1:S, 1:B], Bin)
    @variable(model, y[1:S], Bin)
    #fix.(y, Y)
     # Only x variables now, relaxed to [0,1]
    # @variable(model, 0 <= x[1:S, 1:B] <= 1)
    
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
        sum(x[i, j] * V[M][j] for j in 1:B) >= y[i] * ceil(Int, μ[M][i] * (1 - T[M]))
    )
    @constraint(
        model,
        tol_u[i in 1:S, M in 1:m],
        sum(x[i, j] * V[M][j] for j in 1:B) <= y[i] * floor(Int, μ[M][i] * (1 + T[M]))
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

    return model
end

function write_full_iis(model::Model, filename::String)
    # Compute IIS
    compute_conflict!(model)
    iis_model, reference_map = copy_conflict(model)

    # Open file in write mode with UTF-8 encoding
    open(filename, "w") do f
        # Get all variables
        #println(f, "Variables:")
        #for var in all_variables(iis_model)
        #    println(f, var)
        #end

        println(f, "\nObjective:")
        println(f, objective_function(iis_model))

        println(f, "\nConstraints:")
        # Get all constraint types
        for (F, S) in list_of_constraint_types(iis_model)
            println(f, "\nConstraint type: $F-in-$S")
            # Get all constraints of this type
            for con in all_constraints(iis_model, F, S)
                println(f, con)
            end
        end
        #=
        println(f, "\nBounds:")
        for var in all_variables(iis_model)
            lb = has_lower_bound(var) ? lower_bound(var) : "-∞"
            ub = has_upper_bound(var) ? upper_bound(var) : "∞"
            println(f, "$(name(var)): $lb ≤ $(var) ≤ $ub")
        end
        =#
    end
end

function analyze_splits(x_sol::Matrix{Float64}, epsilon::Float64=1e-6)
    splits = 0
    S, B = size(x_sol)
    
    for j in 1:B
        # For each business unit, count how many non-zero assignments exist
        non_zero = sum(abs(x_sol[i,j]) > epsilon for i in 1:S)
        if non_zero > 1
            splits += 1
            println("Business unit $j is split among $(non_zero) successors")
            # Print the actual split values
            for i in 1:S
                if abs(x_sol[i,j]) > epsilon
                    println("  Successor $i: $(round(x_sol[i,j], digits=4))")
                end
            end
        end
    end
    
    println("\nTotal number of split business units: $splits")
    return splits
end
function plot_ys(Instance, Y, path::String)
    Y = Bool.(Y)
    #plot_font = "Computer Modern";
    default(framestyle=:box, grid=false, tickfontsize=7)
    
    S_coords = Instance.S_coords
    Sk = Instance.Sk
    
    # For each set Sₖ, get the coordinates and corresponding colors from Y
    S₁ = Sk[1]
    S₁_coords = S_coords[S₁, :]
    Y₁_colors = [Y[j] ? :red : :white for j in S₁]
    
    S₂ = Sk[2]
    S₂_coords = S_coords[S₂, :]
    Y₂_colors = [Y[j] ? :red : :white for j in S₂]
    
    S₃ = Sk[3]
    S₃_coords = S_coords[S₃, :]
    Y₃_colors = [Y[j] ? :red : :white for j in S₃]
    
    S₄ = Sk[4]
    S₄_coords = S_coords[S₄, :]
    Y₄_colors = [Y[j] ? :red : :white for j in S₄]
    
    # Plot each set with its corresponding colors
    Plots.scatter!(
        S₁_coords[:, 1],
        S₁_coords[:, 2],
        label = nothing,
        markershape = :hexagon,
        markercolor = Y₁_colors,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        dpi = 500
    )
    
    Plots.scatter!(
        S₂_coords[:, 1],
        S₂_coords[:, 2],
        label = nothing,
        markershape = :diamond,
        markercolor = Y₂_colors,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        dpi = 500
    )
    
    Plots.scatter!(
        S₃_coords[:, 1],
        S₃_coords[:, 2],
        label = nothing,
        markershape = :star5,
        markercolor = Y₃_colors,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        dpi = 500
    )
    
    Plots.scatter!(
        S₄_coords[:, 1],
        S₄_coords[:, 2],
        label = nothing,
        markershape = :pentagon,
        markercolor = Y₄_colors,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        dpi = 500
    )
    
    png(path)
    @debug "Wrote plot and coords"
end

function plot_ys2(Instance, Y, path::String)
    #plot_font = "Computer Modern"
    Y = Bool.(Y)
    
    S_coords = Instance.S_coords
    Sk = Instance.Sk
    
    # Create separate plots for each Sk
    
    # Plot for S₁
    p1 = plot(framestyle=:box, grid=false, tickfontsize=7)
    S₁ = Sk[1]
    S₁_coords = S_coords[S₁, :]
    Y₁_colors = [Y[j] ? :red : :white for j in S₁]
    
    scatter!(
        p1,
        S₁_coords[:, 1],
        S₁_coords[:, 2],
        label = "S₁",
        markershape = :hexagon,
        markercolor = Y₁_colors,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        dpi = 500,
        title = "S₁ Locations"
    )
    
    # Plot for S₂
    p2 = plot(framestyle=:box, grid=false, tickfontsize=7)
    S₂ = Sk[2]
    S₂_coords = S_coords[S₂, :]
    Y₂_colors = [Y[j] ? :red : :white for j in S₂]
    
    scatter!(
        p2,
        S₂_coords[:, 1],
        S₂_coords[:, 2],
        label = "S₂",
        markershape = :diamond,
        markercolor = Y₂_colors,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        dpi = 500,
        title = "S₂ Locations"
    )
    
    # Plot for S₃
    p3 = plot(framestyle=:box, grid=false, tickfontsize=7)
    S₃ = Sk[3]
    S₃_coords = S_coords[S₃, :]
    Y₃_colors = [Y[j] ? :red : :white for j in S₃]
    
    scatter!(
        p3,
        S₃_coords[:, 1],
        S₃_coords[:, 2],
        label = "S₃",
        markershape = :star5,
        markercolor = Y₃_colors,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        dpi = 500,
        title = "S₃ Locations"
    )
    
    # Plot for S₄
    p4 = plot(framestyle=:box, grid=false, tickfontsize=7)
    S₄ = Sk[4]
    S₄_coords = S_coords[S₄, :]
    Y₄_colors = [Y[j] ? :red : :white for j in S₄]
    
    scatter!(
        p4,
        S₄_coords[:, 1],
        S₄_coords[:, 2],
        label = "S₄",
        markershape = :pentagon,
        markercolor = Y₄_colors,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
        dpi = 500,
        title = "S₄ Locations"
    )
    
    # Combine all plots in a 2x2 layout
    final_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(800,800))
    
    # Save the combined plot
    savefig(final_plot, path)
    @debug "Wrote plots for each Sk"
end

function main()
    instance = read_instance(ARGS[1])
    Y, time = pdisp_2(instance)
    println(Y)
    Y_bool = zeros(Int, instance.S)
    for idx in Y
        Y_bool[idx] = 1
    end
    plot_ys2(instance, Y_bool, "pdisp_plotted_all.png")
    #println(Y_bool)
    model_original = build_model(instance)
    set_optimizer(model_original, Gurobi.Optimizer)
    set_time_limit_sec(model_original, 7200)
    write_to_file(model_original, "original_model_cplex.lp")
    optimize!(model_original)
    x_solution = value.(model_original[:x])
    y_sol = value.(model_original[:y])
    X = round.(Int, x_solution)
    Y2 = round.(Int, y_sol)
    model_transport = solve_phase2_model(instance, Y2)
    set_optimizer(model_transport, CPLEX.Optimizer)
    write_to_file(model_transport, "model_transportequality.lp")
    set_time_limit_sec(model_transport, 300)
    optimize!(model_transport)
    println(termination_status(model_transport))
    #write_full_iis(model_transport, "salida_cplex.ilp")
    x_trans = value.(model_transport[:x])
    #println(x_trans)
    #y_trans = value.(model_transport[:y])
    #X_trans = round.(Int, x_trans)
    #Y_trans = round.(Int, y_trans)
    model_transport_grb = solve_phase3_model(instance, Y2, x_trans)
    set_optimizer(model_transport_grb, Gurobi.Optimizer)
    optimize!(model_transport_grb)
    write
    #write_full_iis(model_transport, "salida.ilp")

    #analyze_splits(x_solution)
end
main()