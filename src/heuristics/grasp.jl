include("constructive.jl")
include("ls.jl")
using Types
using Dates
using DelimitedFiles

function grasp(αₗ, αₐ, iters_no_improve, instance)
    start = now()
    bestSol = nothing
    bestWeight = 100000000
    instancia = instance
    iter_improve = 0
    B = instancia.B
    S = instancia.S
    P = instancia.P
    X = Matrix{Int64}(undef, S, B)
    Y = Vector{Int64}(undef, S)
    D = instancia.D
    Weight = 0
    iters = 0
    targets_lower, targets_upper = calculate_targets(instance)
    Xs = []
    Ys = []

    while iters < iters_no_improve
        println(iters)
        Y, time_loc = pdisp_grasp(instance, αₗ)
        X, time_alloc = oppcost_grasp(instance, Y, αₐ)
        Weight = 0
        indices = findall(x -> x == 1, X)
        for indice in indices
            Weight += D[indice]
        end
        @show Weight
        if X ∉ Xs
            println("X nueva")
            push!(Xs, X)
        else
            println("X ya existente")
        end
        if Y ∉ Ys
            println("Y nueva")
            push!(Ys, Y)
        else
            println("Y ya existente")
        end
        iters += 1
    end
end

function oppcost_grasp(instance::Types.Instance, Y, α)
    D = copy(instance.D)
    P = instance.P
    S = instance.S
    B = instance.B
    M = instance.M
    V = instance.V
    M = instance.M
    R = instance.R
    β = instance.β[1]
    X = zeros(Int64, S, B)
    start_alloc = Dates.now()
    count = 0
    not_assigned_y = findall(y -> y == 0, Y)
    #println(not_assigned_y)
    for j in not_assigned_y
        D[j, :] .= -1
    end
    oppMatrix = copy(D)
    for i in 1:B
        minimal = 0
        minimals, _ = minimums(D[:, i], 1)
        minimal = minimals[1]
        oppMatrix[:, i] .= D[:, i] .- minimal
    end
    count = 0
    todos = false
    targets = calculate_target(instance)
    values_matrix = Matrix{Int64}(undef, S, M)
    risk_vec = Vector{Int64}(undef, S)
    n = round(Int, (P / 2)) # falta tweakear
    while !todos
        values, indices = maximums4(oppMatrix, n)
        values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
        constraints = Int64[]
        # aqui calculamos inicialmente las constraints
        # vjm <= μim para cada i en 1:S, para cada m en 1:3, son S*3 constraints
        # Rj > βi para cada i en 1:S, S constraints
        for indice::CartesianIndex in indices
            #X_copy = Matrix{Int64}(undef, S, B)
            #unsafe_copyto!(X_copy, 1, X, 1, S * B)
            col = indice[2]
            row = findfirst(x -> x == 0, oppMatrix[:, col])
            #X[row, col] = 1
            constraints_v = update_constraints(M, V, R, row, col, values_matrix, risk_vec, targets, β)
            push!(constraints, constraints_v)
        end
        picked = CartesianIndex(1, 1)::CartesianIndex{2}
        if all(x -> x == 0, constraints) # si no se violan constraints, agarra el 0 de la columna del costo maximo
            bestVal = values[1]
            cutoff = trunc(Int, (bestVal + (bestVal * α)))
            rcl = [x for x in values if x <= cutoff]
            idx_rcl = [x for x in eachindex(values) if values[x] <= cutoff]
            idx_random = idx_rcl[rand(1:end)] # escoge uno al azar de la rcl
            indice = indices[idx_random]::CartesianIndex{2}
            col = indice[2]::Int64
            row = findfirst(x -> x == 0, oppMatrix[:, col])
            picked = CartesianIndex(row, col)
        else
            if all(x -> x ≠ 0, constraints) # si todas violan constraints
                original_cons, idx = findmin(constraints) # agarra el que viole menos constraints
                cutoff_cons = trunc(Int, (original_cons + (original_cons * α)))
                rcl_inner = [x for x in constraints if x <= cutoff_cons]
                idx_rcl_inner = [x for x in eachindex(constraints) if constraints[x] <= cutoff_cons]
                idx_random = idx_rcl_inner[rand(1:end)] # escoge uno al azar de la rcl
                indice = indices[idx_random] # de nuevo escoge una al azar
                col = indice[2]
                original_row = findfirst(x -> x == 0, oppMatrix[:, col])::Int64
                # queremos buscar en esa columna el siguiente valor más cercano a 0 en oppMatrix
                # al hacerlo, nos acercamos al minimo valor de la matriz de distancias
                busqueda = round(Int, (P / 2)) # a lo mejor cambiarlo despues
                _, idxs_inner::Array{Int64} = minimums(oppMatrix[:, col], busqueda)
                # el primer minimo es el 0 de nuestra localizacion optima
                # lo desechamos para darle variedad a la busqueda
                idxs_inner2::Array{Int64} = idxs_inner[2:end]
                for row in idxs_inner2
                    #X_copy = Matrix{Int64}(undef, S, B)
                    #unsafe_copyto!(X_copy, 1, X, 1, S * B)
                    try
                        #X_copy[row, col] = 1
                        constraints_v2 = update_constraints(M, V, R, row, col, values_matrix, risk_vec, targets, β)
                        if constraints_v2 < original_cons
                            original_cons = constraints_v2
                            original_row = row
                            break
                        end
                    catch
                        # no se porque hay un try catch aqui
                    end
                end
                picked = CartesianIndex(original_row, col)
            else # si hay una que no viola constraints
                for idx in eachindex(constraints)
                    if constraints[idx] == 0 # agarrala
                        indice = indices[idx]
                        col = indice[2]
                        row = findfirst(x -> x == 0, oppMatrix[:, col])
                        picked = CartesianIndex(row, col)
                        break
                    end
                end
            end
        end
        X[picked] = 1
        column = picked[2]::Int64
        oppMatrix[:, column] .= -1 # "apagamos" esa columna
        todos = true
        count += 1
        for col in eachcol(X)
            if all(x -> x == 0, col)
                todos = false
                break
            end
        end
    end
    end_alloc = Dates.now()
    delta_alloc = end_alloc - start_alloc
    alloc_millis = round(delta_alloc, Millisecond)
    return X, alloc_millis.value
end

function pdisp_simple_grasp(d, p, N, α)
    maxdist = 0
    bestpair = (0, 1)
    for i in 1:N
        for j in i+1:N
            if d[i, j] > maxdist
                maxdist = d[i, j]
                bestpair = (i, j)
            end
        end
    end
    P = Set([])
    push!(P, bestpair[1])
    push!(P, bestpair[2])

    while length(P) < p
        maxdist = 0
        vbest = 0
        minimal = Inf
        for v in 1:N
            if v in P
                continue
            end
            mindist = Inf
            for vprime in P
                if d[v, vprime] < mindist
                    mindist = d[v, vprime]
                end
            end
            if mindist > maxdist
                maxdist = mindist
                minimal = mindist
            end
        end
        maxdist = 0
        vbest = 0
        rcl = []
        for v in 1:N
            if v in P
                continue
            end
            mindist = Inf
            for vprime in P
                if d[v, vprime] < mindist
                    mindist = d[v, vprime]
                end
            end
            if mindist > maxdist
                maxdist = mindist
                if mindist <= round(Int, (minimal + minimal * α))
                    push!(rcl, v)
                end
            end
        end
        selected = false
        while !selected
            v_candidate = rand(rcl)
            if v_candidate != 0 && !(v_candidate in P)
                push!(P, v_candidate)
                selected = true
            end
        end
    end
    collection = collect(P)
    return collection
end

function pdisp_grasp(instance, αₗ)
    before_init = Dates.now()
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
    k = 5
    s_coords = instance.S_coords
    metric = Distances.Euclidean()
    d = trunc.(Int, Distances.pairwise(metric, s_coords, dims=1))
    coords_S1 = s_coords[Sk[1], :]
    coords_S2 = s_coords[Sk[2], :]
    coords_S3 = s_coords[Sk[3], :]
    coords_S4 = s_coords[Sk[4], :]
    coords_S5 = s_coords[Sk[5], :]

    d1 = trunc.(Int, Distances.pairwise(metric, coords_S1, dims=1))
    d2 = trunc.(Int, Distances.pairwise(metric, coords_S2, dims=1))
    d3 = trunc.(Int, Distances.pairwise(metric, coords_S3, dims=1))
    d4 = trunc.(Int, Distances.pairwise(metric, coords_S4, dims=1))
    d5 = trunc.(Int, Distances.pairwise(metric, coords_S5, dims=1))
    N1 = length(Sk[1])
    N2 = length(Sk[2])
    N3 = length(Sk[3])
    N4 = length(Sk[4])
    N5 = length(Sk[5])
    p1 = Lk[1]
    p2 = Lk[2]
    p3 = Lk[3]
    p4 = Lk[4]
    p5 = Lk[5]

    pdisp1 = pdisp_simple_grasp(d1, p1, N1, αₗ)
    pdisp2 = pdisp_simple_grasp(d2, p2, N2, αₗ)
    pdisp3 = pdisp_simple_grasp(d3, p3, N3, αₗ)
    pdisp4 = pdisp_simple_grasp(d4, p4, N4, αₗ)
    pdisp5 = pdisp_simple_grasp(d5, p5, N5, αₗ)

    pdisp1_fixed = Sk[1][pdisp1]
    pdisp2_fixed = Sk[2][pdisp2]
    pdisp3_fixed = Sk[3][pdisp3]
    pdisp4_fixed = Sk[4][pdisp4]
    pdisp5_fixed = Sk[5][pdisp5]

    N = S
    pdisp_ok = Set(vcat([pdisp1_fixed, pdisp2_fixed, pdisp3_fixed, pdisp4_fixed, pdisp5_fixed]...))
    if length(pdisp_ok) != P
        count = count_k(pdisp_ok, Sk)
        while length(pdisp_ok) < P
            # Find the node v that maximizes the distance to its closest neighbor in P
            maxdist = 0
            vbest = 0
            rcl = []
            best_dist = 0
            for v in 1:N
                if v in pdisp_ok
                    continue
                end
                k = node_type(v, Sk)
                if count[k] >= Uk[k]
                    continue
                end
                dist = minimum([d[v, vprime] for vprime in pdisp_ok])
                if dist > maxdist
                    maxdist = dist
                    best_dist = dist
                    #vbest = v
                end
            end
            maxdist = 0
            vbest = 0
            rcl = []
            for v in 1:N
                if v in pdisp_ok
                    continue
                end
                k = node_type(v, Sk)
                if count[k] >= Uk[k]
                    continue
                end
                dist = minimum([d[v, vprime] for vprime in pdisp_ok])
                if dist <= round(Int, best_dist + best_dist * αₗ)
                    push!(rcl, v)
                    #vbest = v
                end
            end
            vbest = rand(rcl)
            # If no such node exists, stop the algorithm
            if vbest == 0
                @error "P DISP ERROR"
                break
            end
            # Add the node vbest to the set P and update the counts
            k = node_type(vbest, Sk)
            count[k] += 1
            push!(pdisp_ok, vbest)
        end
    end
    collection = collect(pdisp_ok)
    Y_bool = zeros(Int, instance.S)
    for idx in collection
        Y_bool[idx] = 1
    end
    return Y_bool
    after_init = Dates.now()
    delta_init = after_init - before_init
    delta_init_milli = round(delta_init, Millisecond)
    collection = collect(pdisp_ok)
    return collection, delta_init_milli.value
end

function evalWeight(X, D)
    Weight = 0
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += D[indice]
    end
    return Weight
end

function main_grasp()
    file_name = ARGS[1]
    instance = read_instance(file_name)
    αₗ = 0.8
    αₐ = 0.8
    iters = 20
    grasp(αₗ, αₐ, iters, instance)
    #write_solution(solucion, "solucion_grasp_multithread2threads_nuevo_1_800.jld2")
end

main_grasp()