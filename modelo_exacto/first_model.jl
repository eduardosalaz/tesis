using CPLEX, JuMP


function read_file(path::String)
    buffer = open(path)
    content = readlines(buffer)
    close(buffer)
    dist_ph = [] # distance placeholder
    distance_matrix = Vector{Vector{Float32}} # matrix placeholder
    startl = 0 # line in which distance matrix starts
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
    num_BU, num_Suc = size(distance_matrix)

    # as for now the k param is hardcoded to 5
    k = 5
    Sk = [Array{Int64}(undef,5) for _ in 1:k] # terrible workaround
    Sk_line = endl+2
    j = 1
    for i in Sk_line:length(content)
        if content[i] == "U"
            endl = i
            break
        end
        Sk[j] = [parse(Int, elem) for elem in split(content[i])]
        j = j +1
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
        j = j +1
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
        j = j +1
    end
    endl = endl + 1
    P = parse(Int, content[endl])
    endl = endl + 2
    R = [parse(Int, elem) for elem in split(content[endl])]
    endl = endl + 2
    m = 3 # activities
    V = [Array{Int64}(undef,5) for _ in 1:m]
    j = 1
    for i in endl:length(content)
        if content[i] == "μ"
            endl = i
            break
        end
        V[j] = [parse(Int, value) for value in split(content[i])]
        j = j +1
    end
    endl = endl + 1
    μ = [Array{Int64}(undef,5) for _ in 1:m]
    j = 1
    for i in endl:length(content)
        if content[i] == "T"
            endl = i
            break
        end
        μ[j] = [parse(Int, value) for value in split(content[i])]
        j = j +1
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
        j = j +1
    end
    endl = endl + 1
    α = [parse(Int, thresh) for thresh in split(content[endl])]
    return distance_matrix, Sk, Uk, Lk, P, R, V, μ, T, α, num_BU, num_Suc
    # asserts for dims
end

function build_model(data)

    D = data[1] # matrix float 32
    D = transpose(D) # SO ROWS, COLS NOT COLS, ROWS
    Sk = data[2] # vector of vectors of int
    Uk = data[3] # vector of int
    Lk = data[4] # vector of int
    P = data[5] # int
    R = data[6] # vector of int
    V = data[7] # vector of vectors of int
    μ = data[8] # vector of vectors of int
    T = data[9] # vector of floats
    α = data[10] # vector of int
    num_BU = data[11] # int
    num_Suc = data[12] # int
    m = 3
    k = 5

    model = Model(CPLEX.Optimizer) # THIS IS WHERE THE FUN BEGINS

    @variable(model, x[1:num_Suc, 1:num_BU], Bin) # num suc and num bu, Xᵢⱼ
    @variable(model, y[1:num_Suc], Bin) # Yᵢ

    @objective(model, Min, sum(D .* x)) # Xᵢⱼ * Dᵢⱼ

    @constraint(model, bu_service[j in 1:num_BU], sum(x[i,j] for i in 1:num_Suc) == 1)

    # ∑ᵢ∈S Xᵢⱼ = 1, ∀ j ∈ B

    @constraint(model, use_branch[i in 1:num_Suc, j in 1:num_BU], x[i, j] <= y[i])

    # Xᵢⱼ ≤ Yᵢ , ∀ i ∈ S, j ∈ B

    @constraint(model, sum(y) == P)

    # ∑ i ∈ S Yᵢ = p

    @constraint(model, risk[i in 1:num_Suc, j in 1:num_BU], x[i,j]*R[j] <= α[i])

    # ∑ j ∈ B Xᵢⱼ Rⱼ ≤ αᵢ, ∀ i ∈ S

    #= Constraints below are infeasible

    # infeasibility row bu_service[1] : 0 == 1
    @constraint(
        model,
        tolerance_lower[i in 1:20, m in 1:3, j in 1:8],
        sum((y[j]*μ[m][j]) * (1-T[m])) <= sum((x[i,j]*V[m][i])),
    )

    @constraint(
        model,
        tolerance_upper[i in 1:20, m in 1:3, j in 1:8],
        sum((x[i,j]*V[m][i])) <= sum((y[j]*μ[m][j]) * (1+T[m])),
    )

    # Yᵢμₘʲ(1-tᵐ) ≤ ∑i∈S Xᵢⱼ vᵢᵐ ≤ Yᵢμₘʲ(1+tᵐ) ∀ j ∈ B, m = 1 … 3


    I have no idea how to do this
    #=
    @constraint(
        lower_k[k in 1:5, i in 1:8],
        Lk[k] <= length(Y)
    )
    =#
    # lₖ ≤ ∑i ∈ Sₖ Yᵢ ≤ uₖ, k = 1 … 5
    =#



    return model




end


function main()
    path = "../instancias/instancia.txt"
    #ARGS[1]
    data = read_file(path)
    model = build_model(data)
    optimize!(model)
    println(objective_value(model))
    return model
end

main()
