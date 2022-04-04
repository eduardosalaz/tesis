include("../types/types.jl")
using CPLEX, JuMP, JLD2
using .Types: Instance as inst

function read_jld2(path::String)
    instancia = jldopen(path)
    instance = instancia["instance"]
    build_model(instance)
end


function read_file(path::String)
    buffer = open(path)
    content = readlines(buffer)
    close(buffer)

    dist_ph = [] # distance placeholder
    distance_matrix = Vector{Vector{Float32}} # matrix placeholder
    startl = 0 # line in which content starts
    endl = 0 # line in which it ends
    for (index, line) in enumerate(content)
        if line == "D"
            startl = index + 1
        end
        if line == "S"
            endl = index - 1
            break # computational power isn't to be wasted
        end
    end
    # DISTANCE MATRIX STARTS HERE
    distances_str = content[startl:endl]
    for distance_str in distances_str
        distance_flt = [parse(Float32, distance) for distance in split(distance_str)]
        push!(dist_ph, distance_flt)
    end

    distance_matrix = hcat(dist_ph...) # https://stackoverflow.com/a/52257481
    distance_matrix = trunc.(Int, distance_matrix) # remove later
    num_BU, num_Suc = size(distance_matrix)

    # as for now the k param is hardcoded to 5
    k = 5
    Sk = [Array{Int64}(undef, 5) for _ in 1:k] # terrible workaround
    Sk_line = endl + 2
    j = 1
    for i in Sk_line:length(content)
        if content[i] == "U"
            endl = i
            break
        end
        Sk[j] = [parse(Int, elem) for elem in split(content[i])]
        j = j + 1
    end
    Uk = Array{Int64}(undef, k) # gotta move this to a function so DRY
    Uk_line = endl + 1
    j = 1
    for i in Uk_line:length(content)
        if content[i] == "L"
            endl = i
            break
        end
        Uk[j] = parse(Int, content[i])
        j = j + 1
    end
    Lk = Array{Int64}(undef, k)
    Lk_line = endl + 1
    j = 1
    for i in Lk_line:length(content)
        if content[i] == "P"
            endl = i
            break
        end
        Lk[j] = parse(Int, content[i])
        j = j + 1
    end
    endl = endl + 1
    P = parse(Int, content[endl])
    endl = endl + 2
    R = [parse(Int, elem) for elem in split(content[endl])]
    endl = endl + 2
    m = 2 # activities
    V = [Array{Int64}(undef, 5) for _ in 1:m]
    j = 1
    for i in endl:length(content)
        if content[i] == "μ"
            endl = i
            break
        end
        V[j] = [parse(Int, value) for value in split(content[i])]
        j = j + 1
    end
    endl = endl + 1
    μ = [Array{Int64}(undef, 5) for _ in 1:m]
    j = 1
    for i in endl:length(content)
        if content[i] == "T"
            endl = i
            break
        end
        μ[j] = [parse(Int, value) for value in split(content[i])]
        j = j + 1
    end
    endl = endl + 1
    T = Array{Float64}(undef, m)
    j = 1
    for i in endl:length(content)
        if content[i] == "α"
            endl = i
            break
        end
        T[j] = parse(Float32, content[i]) # consider truncing
        j = j + 1
    end
    endl = endl + 1
    α = [parse(Int, thresh) for thresh in split(content[endl])]
    return distance_matrix, Sk, Uk, Lk, P, R, V, μ, T, α, num_BU, num_Suc
    # asserts for dims
end

function build_model(instance::inst, i::Int64)

    # D = data[1] # matrix float 32
    # D = transpose(D) # SO ROWS, COLS NOT COLS, ROWS
    # Sk = data[2] # vector of vectors of int
    # Uk = data[3] # vector of int
    # Lk = data[4] # vector of int
    # P = data[5] # int
    # R = data[6] # vector of int
    # V = data[7] # vector of vectors of int
    # μ = data[8] # vector of vectors of int
    # T = data[9] # vector of floats
    # α = data[10] # vector of int
    # num_BU = data[11] # int
    # num_Suc = data[12] # int
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

    model = Model(CPLEX.Optimizer) # THIS IS WHERE THE FUN BEGINS
    # set_silent(model)

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

    # write_to_file(model, "modelo8.lp")

    optimize!(model)
    #println(model)

    # instance = Instance(B, S, k, m, P, D, D, D, R, Sk, Uk, Lk, V, μ, T, α)


    X = value.(model[:x])

    X = trunc.(Int, X)

    #println(X)

    Y = value.(model[:y])
    Y = trunc.(Int, Y)
    #println(Y)


    obj_value = trunc(Int, objective_value(model))

    println(obj_value)
    println("------------------$i------------------")

    solucion = Types.Solution(instance, X, Y, obj_value)
    # println(value.(model[:x]))
#    println(value.(y))

    jldsave("outL_100/" * "solution_new_$i" * ".jld2"; solucion)
    return model
end


function main()
    #    path = "../instancias/instancia.txt"
    path = ARGS[1]
#    data = read_file(path)
    model = read_jld2(path)

    return model
end

# main()
