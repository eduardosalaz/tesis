using Types
using JuMP
using Distances
using DelimitedFiles
using Gurobi
using HiGHS
using Dates
using TimerOutputs
using MathOptInterface
using Dates
using Random
using DataStructures

function remove!(V, item)
    deleteat!(V, findall(x -> x == item, V))
end

function constructive(instance, id, init_method, assign_method; withdir=false, dir="")
    instancia = instance
    B = instancia.B
    S = instancia.S
    P = instancia.P
    X = Matrix{Int64}[]
    Y = Vector{Int64}[]
    D = instancia.D
    Weight = 0
    println("hola 4")

    before1 = Dates.now()
    time_init1 = 0
    if init_method == "relax"
        Y, time_init1 = localize_facs(instancia, init_method)
        before1 = Dates.now()
    else
        Y, time = localize_facs(instancia, init_method)
    end

    if init_method ≠ "relax"
        Y_bool = zeros(Int, instancia.S)
        for idx in Y
            Y_bool[idx] = 1
        end
    else
        Y_bool = Y
    end
    after_init = Dates.now()
    delta_init1 = after_init - before1
    if init_method == "relax"
        delta_init1 = time_init1
    end
    println("Y done with time: ", delta_init1)
    if assign_method == "naive"
        #println("NAIVE")
        X = naive_assign_bu(instancia, Y_bool)
    elseif assign_method == "opp"
        println("OPP COST")
        X = oppCostAssignment(Y_bool, instancia)
    elseif assign_method == "queue"
        println("QUEUE")
        X = oppCostQueue(Y_bool, instancia)
    end
    after1 = Dates.now()
    delta_assign = after1 - after_init
    println("X done with time: ", delta_assign)
    delta1 = after1 - before1
    secs1 = round(delta1, Second)
    time1 = secs1.value
    if init_method == "relax"
        time1 += time_init1
    end

    before = Dates.now()
    time_init = 0
    if init_method == "relax"
        Y, time_init = localize_facs(instancia, init_method)
        before = Dates.now()
    else
        Y, time = localize_facs(instancia, init_method)
    end
    if init_method ≠ "relax"
        Y_bool = zeros(Int, instancia.S)
        for idx in Y
            Y_bool[idx] = 1
        end
    else
        Y_bool = Y
    end
    after_init2 = Dates.now()
    delta_init2 = after_init2 - before
    if init_method == "relax"
        delta_init2 = time_init
    end
    println("Y done with time: ", delta_init2)
    if assign_method == "naive"
        #println("NAIVE")
        X = naive_assign_bu(instancia, Y_bool)
    elseif assign_method == "opp"
        println("OPP COST")
        X = oppCostAssignment(Y_bool, instancia)
    end
    after = Dates.now()
    delta_assign2 = after - after_init2
    println("X done with time: ", delta_assign2)
    delta = after - before
    secs = round(delta, Second)
    time = secs.value
    if init_method == "relax"
        time += time_init
    end
    println("X done")
    indices = findall(x -> x == 1, X)
    for indice in indices
        Weight += D[indice]
    end
    if time > time1
        time = time1
    end
    println(time)

    str_path = "sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * init_method * "_" * assign_method
    plot_str_path = str_path * ".png"
    solution_str_path = str_path * ".jld2"

    solution = Types.Solution(instancia, X, Y_bool, Weight, time)
    Types.plot_solution(solution, plot_str_path)
    Types.write_solution(solution, solution_str_path)
    println(isFactible(solution, true))
    return solution
end

function remove!(V, item)
    deleteat!(V, findall(x -> x == item, V))
end

function pdisp_simple(d, p, N)
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
                vbest = v
            end
        end
        if vbest != 0 && !(vbest in P)

            push!(P, vbest)
        end
    end
    collection = collect(P)
    return collection
end

function count_k(P, Sk)
    count = zeros(Int, length(Sk))
    for i in P
        k = node_type(i, Sk)
        count[k] += 1
    end
    return count
end

function node_type(i, Sk)
    for k in eachindex(Sk)
        if i in Sk[k]
            return k
        end
    end
    println("Node $i not found in Sk")
end

function localize_facs(instance, method)
    # k_type = findall(x->x==1, idx_candidate .∈ Sk) # get k type of fac
    if method == "pdisp"
        println("P-DISP")
        return pdisp_2(instance)
    elseif method == "random"
        println("RANDOM")
        return random_init(instance)
    elseif method == "relax"
        println("RELAXATION")
        return relax_init(instance)
    end
end

function pdisp_2(instance)
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

    pdisp1 = pdisp_simple(d1, p1, N1)
    pdisp2 = pdisp_simple(d2, p2, N2)
    pdisp3 = pdisp_simple(d3, p3, N3)
    pdisp4 = pdisp_simple(d4, p4, N4)
    pdisp5 = pdisp_simple(d5, p5, N5)

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
                    vbest = v
                end
            end
            # If no such node exists, stop the algorithm
            if vbest == 0
                @error "PDISP FAILED"
                println("*******************************************************************************************")
                break
            end
            # Add the node vbest to the set P and update the counts
            k = node_type(vbest, Sk)
            count[k] += 1
            push!(pdisp_ok, vbest)
        end
    end
    after_init = Dates.now()
    delta_init = after_init - before_init
    delta_init_milli = round(delta_init, Millisecond)
    collection = collect(pdisp_ok)
    return collection, delta_init_milli.value
end

function minimums(matrix::Matrix, n)::Tuple{Vector{Int64},Vector{CartesianIndex{2}}}
    type = eltype(matrix)
    vals = fill(10000000000000, n)
    arr = Array{Tuple{type,CartesianIndex}}(undef, n)
    indices = Array{Int64}(undef, n)
    @inbounds for i ∈ axes(matrix, 1), j ∈ axes(matrix, 2)
        biggest, index = findmax(vals)
        if matrix[i, j] < biggest
            arr[index] = matrix[i, j], CartesianIndex(i, j)
            vals[index] = matrix[i, j]
        end
    end
    arr = sort(arr, by=x -> x[1])
    vals = [x[1] for x in arr]
    indices = [x[2] for x in arr]
    return vals, indices
end

function minimums(vec::Vector, n)::Tuple{Vector{Int64},Vector{Int64}}
    type = eltype(vec)
    vals = fill(10000000000000, n)
    arr = Array{Tuple{type,Int64}}(undef, n)
    indices = Array{Int64}(undef, n)
    @inbounds for i ∈ eachindex(vec)
        if vec[i] > 0
            biggest, index = findmax(vals)
            if vec[i] < biggest
                arr[index] = vec[i], i
                vals[index] = vec[i]
            end
        end
    end
    arr = sort(arr, by=x -> x[1])
    vals = [x[1] for x in arr]
    indices = [x[2] for x in arr]
    return vals, indices
end
#132407

function maximums3(M, n)
    v = vec(M)
    l = length(v)
    ix = [1:l;]
    partialsortperm!(ix, v, (l-n+1):l, initialized=true)
    indices = CartesianIndices(M)[ix[(l-n+1):l]]
    return reverse!(indices)
end


function calculate_target(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    targets = Vector{Float32}(undef, M)
    for m in 1:M
        targets[m] = 1 * μ[m][1] * (1 + T[m]) # corregir esto
    end
    return targets
end

function calculate_targets_upper(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    targets = Vector{Float32}(undef, M)
    for m in 1:M
        targets[m] = 1 * μ[m][1] * (1 + T[m]) # corregir esto
    end
    return targets
end

function calculate_targets_lower(instance)
    M = instance.M
    μ = instance.μ
    T = instance.T
    targets = Vector{Float32}(undef, M)
    for m in 1:M
        targets[m] = 1 * μ[m][1] * (1 - T[m]) # corregir esto
    end
    return targets
end

function start_constraints(S, B, M, V, R, X, values_matrix, risk_vec)
    for i in 1:S
        for m in 1:M
            values_matrix[i, m] = sum(X[i, j] * V[m][j] for j in 1:B)
        end
        risk_vec[i] = sum(X[i, j] * R[j] for j in 1:B)
    end
    return values_matrix, risk_vec
end

function update_constraints(M, V, R, i, j, values_matrix, risk_vec, targets, β)
    constraints = 0
    for m in 1:M
        values_matrix[i, m] += V[m][j]
        if values_matrix[i, m] > targets[m]
            constraints += 1
        end
    end
    risk_vec[i] += R[j]
    if risk_vec[i] > β
        constraints += 1
    end
    return constraints
end

function encode_key(client_index, facility_index, max_facilities)
    return client_index * max_facilities + facility_index
end

function decode_key(key, max_facilities)
    client_index = div(key, max_facilities)
    facility_index = key % max_facilities
    return client_index, facility_index
end

function compute_assignments_and_opportunity_costs(D::Matrix{Int64}, Y::Vector{Int}, N::Int)
    _, num_clients = size(D)
    best_assignments = Dict{Int,Vector{Int64}}()
    pq = PriorityQueue{Int,Int64}(Base.Order.Reverse)

    for j in 1:num_clients
        # Use a temporary array to store the facility opportunity costs for this client
        costs = Tuple{Int64,Int}[]
        for (i, yi) in enumerate(Y)
            if yi == 1
                push!(costs, (D[i, j], i))
            end
        end
        # Sort the costs
        sort!(costs)
        # Store the top N assignments for this client
        best_assignments[j] = best_assignments[j] = [cost[2] for cost in costs[1:N]] # extrae el indice nada mas
        # Calculate the opportunity cost as the largest difference among all possible assignments
        opp_cost = costs[end][1] - costs[1][1]  # Difference between largest and smallest costs
        # Use the smallest facility index for this client to form the unique key
        enqueue!(pq, j, opp_cost)
    end
    return best_assignments, pq
end



"""
1. Calculate Opportunity Costs
    Initialize oppCostMatrix with zeros
    FOR each facility i
        FOR each client j
            Calculate opportunity cost for client j with facility i
            Insert into oppCostMatrix[i][j]

2. Compute Best N Assignments
    Initialize best_assignments as an empty dictionary
    FOR each client j
        sorted_facilities = sort facilities based on oppCostMatrix[:,j]
        best_assignments[j] = first N entries of sorted_facilities

3. Prepare the priority queue with clients (based on their highest opportunity costs)
    Initialize priority_queue with all clients, priority being the highest opp cost from oppCostMatrix

4. Assignment Loop
    Initialize assigned_clients as an empty set
    WHILE priority_queue is not empty AND length of assigned_clients < number of clients
        client = dequeue priority_queue
        IF client not in assigned_clients
            FOR each facility in best_assignments[client]
                IF can_assign(facility)
                    assign(client, facility)
                    ADD client to assigned_clients
                    UPDATE capacity of facility
                    IF facility is full
                        update_best_assignments_for_all_clients(facility, best_assignments)
                    BREAK out of loop, move to next client

    1. Assignment Loop
    assigned_clients = Set() # An empty set to keep track of clients who have been assigned
    assignment_dict = Dict() # A dictionary to store which facility each client is assigned to

    WHILE !isempty(priority_queue) AND length(assigned_clients) < number_of_clients
        client, _ = dequeue(priority_queue) # Get the client with the highest opportunity cost

        IF client not in assigned_clients
            FOR facility in best_assignments[client]
                IF can_assign(facility) # Check constraints or facility capacity
                    assignment_dict[client] = facility
                    ADD client to assigned_clients
                    DECREASE capacity of facility
                    
                    IF facility is full
                        update_best_assignments_for_all_clients(facility, best_assignments)
                    END IF

                    BREAK # Move to the next client in the priority queue
                END IF
            END FOR
        END IF
    END WHILE

    RETURN assignment_dict

5. Handle Unassigned Clients
    FOR each client in clients
        IF client not in assigned_clients
            assign client to any available facility (based on some fallback heuristic)

FUNCTION update_best_assignments_for_all_clients(facility, best_assignments)
    FOR each client in best_assignments.keys()
        IF facility in best_assignments[client]
            REMOVE facility from best_assignments[client]
            IF length of best_assignments[client] < N
                Find the next best facility not in best_assignments[client] for this client
                ADD this new facility to best_assignments[client]
"""
function oppCostQueue(Y, instance::Types.Instance)
    D = copy(instance.D)
    X = zeros(Int, size(D))
    targets = calculate_targets_lower(instance)
    β = instance.β[1]
    S = instance.S
    B = instance.B
    M = instance.M
    V = instance.V
    P = instance.P
    R = instance.R
    values_matrix, risk_vec = start_constraints(S, B, M, V, R, X, zeros(S, M), zeros(S))
    N = instance.P # podemos hacerlo en proporcion a P, no necesariamente tiene que ser P
    best_assignments, pq = compute_assignments_and_opportunity_costs(D, Y, N)
    full_centers = Set()
    assigned_bus = Set()
    unassigned_bus = Set(collect(1:B))
    while !isempty(pq) && length(assigned_bus) < B
        bu = dequeue!(pq)
        if bu ∉ assigned_bus
            for center in best_assignments[bu]
                full = false
                fulls_m = zeros(Int, M)
                if center ∉ full_centers
                    for m in 1:M
                        values_matrix[center, m] += V[m][bu]
                        if values_matrix[center, m] > targets[m]
                            fulls_m[m] = 1
                        end
                    end
                    if all(x -> x == 1, fulls_m)
                        push!(full_centers, center)
                    end
                    risk_vec[center] += R[bu]
                    if risk_vec[center] > β
                        push!(full_centers, center)
                    end
                    push!(assigned_bus, bu)
                    pop!(unassigned_bus, bu)
                    X[center, bu] = 1
                    break
                end
            end
        end
    end
    todos = true
    count = 0
    for col in eachcol(X)
        if all(x -> x == 0, col)
            count += 1
            todos = false
        end
    end
    X = handle_unassigned_clients2(X, instance, best_assignments, unassigned_bus, values_matrix, risk_vec)
    return X
end

# Function to calculate how much an assignment violates constraints
function calculate_violation(facility, client, targets_upper, V, M, R, β, values_matrix, risk_vec)
    violation = 0
    for m in 1:M
        excess = values_matrix[facility, m] + V[m][client] - targets_upper[m]
        if excess > 0
            violation += excess
        end
    end
    excess_risk = risk_vec[facility] + R[client] - β
    if excess_risk > 0
        violation += excess_risk
    end
    return violation
end



function handle_unassigned_clients2(X, instance, best_assignments, unassigned_clients, values_matrix, risk_vec)
    targets_upper = calculate_targets_upper(instance)
    S = instance.S
    M = instance.M
    V = instance.V
    β = instance.β[1]
    R = instance.R

    # Only loop over unassigned clients
    for client in unassigned_clients
        assigned = false
        # Sort facilities for this client based on their remaining risk capacity
        sorted_best_assignments = sort(best_assignments[client], by=f -> β - risk_vec[f], rev=true)
        for facility in sorted_best_assignments
            potential_assignment_valid = true
            for m in 1:M
                if values_matrix[facility, m] + V[m][client] > targets_upper[m]
                    potential_assignment_valid = false
                    break
                end
            end
            if risk_vec[facility] + R[client] > β
                potential_assignment_valid = false
            end
            if potential_assignment_valid
                X[facility, client] = 1
                assigned = true
                delete!(unassigned_clients, client)
                for m in 1:M
                    values_matrix[facility, m] += V[m][client]
                end
                risk_vec[facility] += R[client]
                break  # Assign to the first valid facility and then exit the inner loop
            end
        end
        if !assigned
            min_violation = 10000000
            least_violating_facility = -1
            for facility in sorted_best_assignments
                violation = calculate_violation(facility, client, targets_upper, V, M, R, β, values_matrix, risk_vec)
                if violation < min_violation
                    min_violation = violation
                    least_violating_facility = facility
                end
            end
            # Assign the client to the least violating facility
            X[least_violating_facility, client] = 1
            for m in 1:M
                values_matrix[least_violating_facility, m] += V[m][client]
            end
            risk_vec[least_violating_facility] += R[client]
            println("Client $client was assigned to facility $least_violating_facility with a violation of $min_violation.")
        end
    end

    println("Number of unassigned clients: $(length(unassigned_clients))")
    return X
end

#=
 \STATE \textbf{Input:} Distance matrix $d$, Decision Variable $Y$, $\lvert B \rvert$, $\lvert S \rvert$
        \STATE \textbf{Output:} Opportunity Cost Matrix
         \Procedure{OppCostMatrix}
            \FOR{$j$ \in $1$ \to $ \lvert B \rvert$}
                \STATE $minimal \gets$ \textbf{argmin}($d_[:, j]$)
                \FOR{$i$ \in $1$ \to $\lvert S \rvert$}
                    \STATE $oppMatrix[j, i] \gets D[j, i] - minimal$
                \ENDFOR
            \ENDFOR
            \STATE \textbf{return} $oppMatrix$
        \EndProcedure
=#

function oppCostAssignment(Y, instance::Types.Instance)
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
    count = 0
    not_assigned_y = findall(y -> y == 0, Y)
    #println(not_assigned_y)
    for j in not_assigned_y
        D[j, :] .= -1
    end
    oppMatrix = copy(D)
    for i in 1:B
        minimal = 0
        minimals, _ = minimums(D[:, i], 1) # cambiar por findmin en 1 solo
        minimal = minimals[1]
        oppMatrix[:, i] .= D[:, i] .- minimal
    end
    count = 0
    todos = false
    targets = calculate_target(instance)
    values_matrix = Matrix{Float32}(undef, S, M)
    risk_vec = Vector{Float32}(undef, S)
    n = round(Int, (P / 2)) # falta tweakear
    while !todos
        indices::Vector{CartesianIndex{2}} = maximums3(oppMatrix, n)
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
            indice = indices[1]::CartesianIndex{2}
            col = indice[2]::Int64
            row = findfirst(x -> x == 0, oppMatrix[:, col])
            picked = CartesianIndex(row, col)
        else
            if all(x -> x ≠ 0, constraints) # si todas violan constraints
                original_cons, idx = findmin(constraints) # agarra el que viole menos constraints
                indice = indices[idx]
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
    return X
end

function pdisp(instance)
    P = instance.P
    s_coords = instance.S_coords
    S = instance.S
    metric = Distances.Euclidean()
    s_distances = trunc.(Int, Distances.pairwise(metric, s_coords, dims=1))
    idx_1, idx_2 = Tuple(argmax(s_distances)) # two furthest facs
    S_sol = [idx_1, idx_2] # S is our solution of our p-disp problem
    T = collect(1:S) # nodes not in our solution
    # se puede usar la suma de xⱼ a S_sol
    # o usar la facility mas cercana de xⱼ
    for s in S_sol
        remove!(T, s)
    end

    while length(S_sol) < P
        distances = Tuple{Int64,Int64}[]
        idx_max_dist = 0
        for i in T
            distance = 0
            for j in S
                distance += s_distances[i, j]
            end
            tupla = (distance, i)
            push!(distances, tupla)
        end
        distances_max = [distance[1] for distance in distances]::Array{Int64}
        indexes = [distance[2] for distance in distances]::Array{Int64}
        idx_max = argmax(distances_max)::Int64
        idx_max_dist = indexes[idx_max]
        remove!(T, idx_max_dist)
        push!(S_sol, idx_max_dist)
    end
    return S_sol
end

function random_init(instance)
    # tengo que agarrar los parametros para asignar una factible
    P = instance.P
    S = instance.S
    Y = Random.shuffle(collect(1:S))
    Y = Y[1:P]
    return Y
end

function relax_init(instance)
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
    m = instance.M
    k = instance.K

    model = JuMP.Model(Gurobi.Optimizer) # THIS IS WHERE THE FUN BEGINS

    if solver_name(model) == "Gurobi"
        MOI.set(model, MOI.RawOptimizerAttribute("SolutionLimit"), 1) # retorna la primera factible encontrada, solo funciona para gurobi    
    end

    JuMP.@variable(model, x[1:S, 1:B], lower_bound = 0, upper_bound = 1)
    # num suc and num bu, Xᵢⱼ
    JuMP.@variable(model, y[1:S], Bin)
    # Yᵢ

    JuMP.@objective(model, Min, sum(D .* x))
    # Xᵢⱼ * Dᵢⱼ

    JuMP.@constraint(model, bu_service[j in 1:B], sum(x[i, j] for i in 1:S) == 1)

    # ∑ᵢ∈S Xᵢⱼ = 1, ∀ j ∈ B

    JuMP.@constraint(model, use_branch[j in 1:B, i in 1:S], x[i, j] <= y[i])

    # Xᵢⱼ ≤ Yᵢ , ∀ i ∈ S, j ∈ B

    JuMP.@constraint(model, cardinality, sum(y) == P)

    # ∑ i ∈ S Yᵢ = p

    JuMP.@constraint(model, risk[j in 1:B, i in 1:S], x[i, j] * R[j] <= β[i])

    # ∑ j ∈ B Xᵢⱼ Rⱼ ≤ βᵢ, ∀ i ∈ S

    JuMP.@constraint(
        model,
        tol_l[i in 1:S, M in 1:m],
        y[i] * μ[M][i] * (1 - T[M]) <= sum(x[i, j] * V[m][j] for j in 1:B),
    )

    JuMP.@constraint(
        model,
        tol_u[i in 1:S, M in 1:m],
        sum(x[i, j] * V[m][j] for j in 1:B) <= y[i] * μ[M][i] * (1 + T[M]),
    )

    # Yᵢμₘⁱ(1-tᵐ) ≤ ∑i∈S Xᵢⱼ vⱼᵐ ≤ Yᵢμₘʲ(1+tᵐ) ∀ j ∈ B, m = 1 … 3

    JuMP.@constraint(
        model,
        low_k[K in 1:k],
        Lk[K] <= sum(y[i] for i in Sk[K]),
    )

    JuMP.@constraint(
        model,
        upp_k[K in 1:k],
        sum(y[i] for i in Sk[K]) <= Uk[K],
    )

    JuMP.set_silent(model)
    JuMP.optimize!(model)
    tiempo = MOI.get(model, MOI.SolveTimeSec())
    time_int = round(Int, tiempo)
    Y = round.(Int, JuMP.value.(model[:y]))
    return Y, time_int
end

function isFactible(solution::Types.Solution, verbose=true)
    number_constraints_violated = 0
    instance = solution.Instance
    Sk = instance.Sk
    Lk = instance.Lk
    Uk = instance.Uk
    Y = solution.Y
    X = solution.X
    Lk = instance.Lk
    Uk = instance.Uk
    Sk = instance.Sk
    μ = instance.μ
    T = instance.T
    M = instance.M
    S = instance.S
    P = instance.P
    B = instance.B
    β = instance.β
    R = instance.R
    K = instance.K
    usables_i = Set(findall(==(1), Y))
    V = instance.V

    counts_k = []

    for y in eachindex(Y)
        if Y[y] == 1
            assignments_y = X[y, :]
            if !any(x -> x == 1, assignments_y)
                if verbose
                    println("Violando asignación de Y en: $y")
                end
                number_constraints_violated += 1
                # penalizar más aquí?
            end
        end
    end

    if sum(Y) != P
        if verbose
            println("Violando número de centros asignados ", sum(Y), " $P")
        end
        number_constraints_violated += 1
    end

    for j in 1:B
        if sum(X[i, j] for i in 1:S) ≠ 1
            if verbose
                println("Violando servicio a la BU: ", j)
                number_constraints_violated += 1
            end

        end
    end

    for j in 1:B
        assigned = findfirst(==(1), X[:, j])
        if assigned ∉ usables_i
            println("Violando $j asignada a una y no abierta $assigned")
            number_constraints_violated += 1
        end
    end

    for k_type in 1:K
        indices_k_type = Sk[k_type]
        count_k_type = 0
        for indice in indices_k_type
            if Y[indice] == 1
                count_k_type += 1
            end
        end
        push!(counts_k, count_k_type)
    end

    for k in 1:K
        if counts_k[k] < Lk[k]
            if verbose
                println("Violando Lk en $k")
            end
            number_constraints_violated += 1
        end
        if counts_k[k] > Uk[k]
            if verbose
                println("Violando Uk en $k")
            end
            number_constraints_violated += 1
        end
    end

    for i in 1:S
        for m in 1:M
            if !(Y[i] * μ[m][i] * (1 - T[m]) <= sum(X[i, j] * V[m][j] for j in 1:B))
                if verbose
                    println("violando V inferior en i: $i y m: $m")
                    println("μ: ",  Y[i] * μ[m][i] * (1 - T[m]))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
            if !(sum(X[i, j] * V[m][j] for j in 1:B) <= Y[i] * μ[m][i] * (1 + T[m]))
                if verbose
                    println("violando V superior en i: $i y m: $m")
                    println("μ: ",  [i] * μ[m][i] * (1 + T[m]))
                    println("V: ", sum(X[i, j] * V[m][j] for j in 1:B))
                end
                number_constraints_violated += 1
            end
        end
    end

    for i in 1:S
        if sum(X[i, j] * R[j] for j in 1:B) > β[i]
            if verbose
                println("violando riesgo en $i")
                println("β: ", β[i])
                println("R: ", sum(X[i, j] * R[j] for j in 1:B))
            end
            number_constraints_violated += 1
        end
    end
    # println("Restricciones violadas: $number_constraints_violated")
    # intercambiar nodo i por nodo j 
    # mover nodo de territorio i a territorio j
    # apagar un branch y prender otro branch
    # checar ILS ??
    if number_constraints_violated ≠ 0
        return false, number_constraints_violated
    else
        return true, number_constraints_violated
    end
end

function main_constructive(init_method, assign_method; path="inst", read_file=true, instance_obj=nothing, id=0)
    println("hola 2")
    instance = 1 # para traerlo al scope
    if read_file
        pattern = Regex("[t][_]\\d{1,3}")
        index = findfirst(pattern, path)
        almost_number = path[index]
        _, number = split(almost_number, "_")
        instance = Types.read_instance(path)
    else
        instance = instance_obj
        number = id
    end
    B = instance.B
    P = instance.P
    S = instance.S
    println("hola 3")

    solution = constructive(instance, number, init_method, assign_method)
    println(isFactible(solution))
    sol_path = "sol_$number" * "_$B" * "_$S" * "_$P" * "_$init_method" * "_$assign_method" * "_rev.jld2"
    write_solution(solution, sol_path)
    println("\a")
    return solution
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("hola")
    main_constructive(ARGS[2], ARGS[3]; path=ARGS[1])
end
