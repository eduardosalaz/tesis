using JuMP
using Gurobi
using LinearAlgebra
using Plots

# Create model
function build_multi_p_dispersion()
    # Data
    p = 40  # Number of points to select
    
    # Point data: x, y, type
    data = [
        (255, 75, 2),
        (318, 314, 3),
        (6, 213, 1),
        (433, 224, 1),
        (328, 377, 2),
        (90, 70, 1),
        (272, 215, 3),
        (485, 175, 2),
        (469, 250, 2),
        (395, 25, 4),
        (97, 157, 1),
        (343, 211, 2),
        (230, 375, 5),
        (39, 446, 4),
        (489, 446, 2),
        (75, 325, 3),
        (290, 121, 4),
        (157, 150, 1),
        (254, 377, 5),
        (61, 108, 2),
        (401, 368, 4),
        (293, 407, 5),
        (460, 298, 1),
        (7, 397, 5),
        (69, 450, 3),
        (56, 418, 2),
        (189, 58, 4),
        (123, 406, 2),
        (374, 359, 4),
        (250, 260, 3),
        (222, 102, 1),
        (143, 360, 1),
        (207, 130, 1),
        (299, 303, 3),
        (337, 273, 5),
        (38, 182, 5),
        (306, 18, 4),
        (265, 136, 5),
        (469, 400, 4),
        (486, 267, 1),
        (368, 110, 1),
        (225, 171, 1),
        (111, 31, 5),
        (472, 68, 1),
        (412, 496, 5),
        (280, 18, 2),
        (7, 112, 4),
        (464, 376, 1),
        (359, 463, 2),
        (167, 64, 5),
        (188, 311, 1),
        (185, 116, 5),
        (180, 255, 1),
        (130, 254, 2),
        (321, 438, 3),
        (479, 172, 2),
        (250, 347, 4),
        (188, 428, 4),
        (312, 361, 2),
        (171, 432, 4),
        (73, 130, 4),
        (199, 230, 3),
        (375, 120, 2),
        (110, 89, 2),
        (240, 51, 1),
        (456, 111, 4),
        (222, 253, 5),
        (187, 175, 5),
        (304, 437, 3),
        (289, 320, 2),
        (205, 391, 2),
        (24, 202, 4),
        (157, 437, 3),
        (437, 464, 3),
        (463, 124, 4),
        (420, 469, 4),
        (294, 381, 5),
        (274, 99, 5),
        (397, 34, 2),
        (32, 349, 5),
        (56, 355, 2),
        (391, 368, 2),
        (455, 88, 4),
        (96, 471, 1),
        (238, 246, 2),
        (279, 36, 1),
        (171, 158, 5),
        (196, 298, 4),
        (58, 381, 2),
        (144, 28, 5),
        (316, 101, 4),
        (350, 237, 4),
        (416, 53, 5),
        (176, 324, 5),
        (187, 191, 2),
        (298, 455, 3),
        (54, 188, 1),
        (383, 402, 4),
        (295, 497, 2),
        (272, 431, 1),
        (149, 157, 4),
        (486, 333, 3),
        (198, 122, 3),
        (210, 104, 3),
        (220, 484, 3),
        (399, 412, 5),
        (190, 340, 2),
        (158, 11, 3),
        (492, 188, 2),
        (114, 63, 1),
        (234, 227, 5),
        (83, 98, 1),
        (136, 367, 1),
        (400, 482, 5),
        (495, 176, 3),
        (153, 247, 1),
        (144, 490, 1),
        (181, 384, 1),
        (438, 422, 5),
        (368, 479, 4),
        (426, 399, 5),
        (72, 194, 3),
        (307, 12, 5),
        (251, 276, 2),
        (404, 318, 4),
        (473, 313, 3),
        (125, 209, 3),
        (325, 286, 3),
        (403, 491, 5),
        (374, 318, 3),
        (162, 155, 1),
        (378, 98, 3),
        (47, 145, 4),
        (364, 193, 2),
        (75, 108, 4),
        (294, 291, 1),
        (51, 79, 3),
        (402, 338, 1),
        (123, 68, 3),
        (264, 257, 1),
        (241, 264, 4),
        (24, 154, 3),
        (447, 337, 5),
        (332, 273, 4),
        (405, 274, 4),
        (153, 454, 1),
        (69, 145, 3),
        (130, 429, 5),
        (500, 297, 4),
        (282, 110, 3),
        (491, 26, 2),
        (225, 425, 5),
        (37, 302, 3),
        (103, 484, 5),
        (276, 106, 5),
        (498, 19, 2),
        (30, 387, 3),
        (310, 498, 2),
        (309, 115, 2),
        (161, 278, 4),
        (44, 349, 3),
        (94, 156, 4),
        (290, 150, 4),
        (244, 209, 5),
        (280, 60, 1),
        (496, 273, 4),
        (435, 173, 1),
        (445, 34, 5),
        (187, 392, 1),
        (446, 286, 4),
        (352, 231, 5),
        (94, 481, 2),
        (194, 465, 2),
        (345, 200, 4),
        (264, 183, 2),
        (487, 137, 2),
        (383, 499, 3),
        (497, 438, 3),
        (183, 181, 5),
        (436, 149, 1),
        (213, 463, 2),
        (243, 283, 1),
        (399, 21, 3),
        (378, 142, 1),
        (67, 55, 2),
        (164, 161, 4),
        (461, 220, 2),
        (410, 444, 5),
        (66, 63, 1),
        (323, 252, 1),
        (55, 165, 4),
        (324, 490, 1),
        (150, 456, 3),
        (354, 441, 4),
        (340, 360, 5),
        (93, 405, 1),
        (182, 99, 2),
        (261, 13, 3),
        (355, 431, 1),
        (21, 108, 5),
        (251, 350, 1),
        (244, 367, 2),
        (100, 57, 1),
        (482, 301, 4),
        (319, 409, 2),
        (343, 158, 2),
        (342, 475, 5),
        (182, 7, 4),
        (239, 327, 4),
        (102, 87, 1),
        (360, 319, 4),
        (325, 201, 2),
        (171, 398, 4),
        (158, 345, 5),
        (304, 278, 2),
        (295, 140, 1),
        (123, 129, 1),
        (294, 372, 3),
        (377, 267, 3),
        (209, 26, 4),
        (338, 78, 4),
        (3, 423, 4),
        (250, 374, 4),
        (443, 364, 3),
        (91, 230, 2),
        (386, 359, 3),
        (268, 361, 4),
        (143, 417, 3),
        (105, 292, 4),
        (482, 336, 2),
        (1, 96, 4),
        (350, 311, 3),
        (469, 14, 3),
        (221, 410, 4),
        (404, 219, 4),
        (67, 410, 4),
        (86, 480, 1),
        (286, 355, 3),
        (117, 227, 3),
        (308, 212, 1),
        (185, 57, 5),
        (19, 235, 4),
        (478, 83, 1),
        (494, 374, 1),
        (104, 159, 3),
        (455, 205, 4),
        (189, 22, 2),
        (253, 395, 2),
        (191, 296, 1),
        (69, 7, 5),
        (491, 315, 2),
        (398, 349, 2),
        (327, 57, 1),
        (210, 86, 2),
        (159, 443, 4),
        (228, 203, 1),
        (160, 288, 1),
        (202, 163, 3),
        (377, 268, 4),
        (74, 355, 4),
        (317, 452, 3),
        (369, 301, 2),
        (393, 267, 4),
        (398, 458, 4),
        (145, 404, 5),
        (40, 345, 1),
        (413, 288, 1),
        (244, 136, 5),
        (250, 465, 5),
        (132, 253, 5),
        (41, 168, 4),
        (396, 37, 5),
        (210, 79, 2),
        (500, 291, 3),
        (11, 155, 1),
        (263, 52, 3),
        (270, 358, 3),
        (439, 340, 2),
        (364, 152, 1),
        (163, 249, 5),
        (437, 241, 2),
        (354, 246, 3),
        (2, 284, 2),
        (241, 91, 3),
        (262, 412, 4),
        (225, 439, 5),
        (349, 333, 2),
        (121, 143, 4),
        (486, 483, 1),
        (443, 205, 1),
        (387, 382, 2),
        (319, 304, 3),
        (11, 149, 3),
        (387, 18, 4),
        (284, 273, 2),
        (475, 489, 5),
        (56, 272, 1),
        (180, 77, 4),
        (246, 199, 3),
        (152, 160, 1)
    ]
    n = length(data)  # Number of points
    I = 1:n           # Set of points
    
    # Calculate Euclidean distances between all pairs of points
    d = zeros(n, n)
    for i in I
        for j in I
            d[i,j] = round(sqrt((data[i][1] - data[j][1])^2 + (data[i][2] - data[j][2])^2))
        end
    end
    
    # Maximum distance
    Dmax = maximum(d)
    
    # Get unique types
    K = sort(unique([point[3] for point in data]))
    
    # Calculate type proportions and bounds
    proportion = Dict(k => count(p -> p[3] == k, data) / n for k in K)
    lb = Dict(k => floor(Int, p * proportion[k]) for k in K)
    ub = Dict(k => ceil(Int, p * proportion[k]) for k in K)
    
    println("Type proportions: ", proportion)
    println("Lower bounds: ", lb)
    println("Upper bounds: ", ub)
    
    # Model A: Multi p-Dispersion with sumsum objective
    function multi_pdp_sum()
        model = Model(Gurobi.Optimizer)
        
        # Variables
        @variable(model, y[I], Bin)  # 1 if point i is selected
        @variable(model, x[i=I, j=I; i < j], Bin)  # 1 if both i and j are selected
        @variable(model, w[K], Int)  # Number of points of type k selected
        
        # Objective: Maximize sum of distances between selected points
        @objective(model, Max, sum(d[i,j] * x[i,j] for i in I for j in I if i < j))
        
        # Constraint: Select exactly p points
        @constraint(model, sum(y) == p)
        
        # Constraints: Relationship between x and y variables
        for i in I
            for j in I
                if i < j
                    @constraint(model, x[i,j] <= y[i])
                    @constraint(model, x[i,j] <= y[j])
                    @constraint(model, x[i,j] >= y[i] + y[j] - 1)
                end
            end
        end
        
        # Constraint: Relationship between y and w variables
        for k in K
            @constraint(model, w[k] == sum(y[i] for i in I if data[i][3] == k))
        end
        
        # Constraints: Type balance
        for k in K
            @constraint(model, w[k] >= lb[k])
            @constraint(model, w[k] <= ub[k])
        end
        
        return model
    end
    
    # Model B: Multi p-Dispersion with minmin objective
    function multi_pdp_min()
        model = Model(Gurobi.Optimizer)
        
        # Variables
        @variable(model, y[I], Bin)  # 1 if point i is selected
        @variable(model, u >= 0)     # Minimum distance between any two selected points
        @variable(model, w[K], Int)  # Number of points of type k selected
        
        # Objective: Maximize the minimum distance
        @objective(model, Max, u)
        
        # Constraint: Linearization of minimum distance
        for i in I
            for j in I
                if i < j
                    @constraint(model, u <= d[i,j] + Dmax * (2 - (y[i] + y[j])))
                end
            end
        end
        
        # Constraint: Select exactly p points
        @constraint(model, sum(y) == p)
        
        # Constraint: Relationship between y and w variables
        for k in K
            @constraint(model, w[k] == sum(y[i] for i in I if data[i][3] == k))
        end
        
        # Constraints: Type balance
        for k in K
            @constraint(model, w[k] >= lb[k])
            @constraint(model, w[k] <= ub[k])
        end
        
        return model
    end
    
    # Model C: p-Dispersion with minmin objective (no type constraints)
    function pdp_min()
        model = Model(Gurobi.Optimizer)
        
        # Variables
        @variable(model, y[I], Bin)  # 1 if point i is selected
        @variable(model, u >= 0)     # Minimum distance between any two selected points
        
        # Objective: Maximize the minimum distance
        @objective(model, Max, u)
        
        # Constraint: Linearization of minimum distance
        for i in I
            for j in I
                if i < j
                    @constraint(model, u <= d[i,j] + Dmax * (2 - (y[i] + y[j])))
                end
            end
        end
        
        # Constraint: Select exactly p points
        @constraint(model, sum(y) == p)
        
        return model
    end
    
    return pdp_min, multi_pdp_min, multi_pdp_sum, data, d, lb, ub, K
end

function multi_type_pdp_heuristic(data, p; objective="sumsum", local_search=true)
    n = length(data)
    
    # Calculate distances between all points
    d = zeros(n, n)
    for i in 1:n
        for j in 1:n
            if i != j
                d[i,j] = round(sqrt((data[i][1] - data[j][1])^2 + (data[i][2] - data[j][2])^2))
            end
        end
    end
    
    # Get unique types and their proportions
    types = sort(unique([point[3] for point in data]))
    type_counts = Dict(k => count(p -> p[3] == k, data) for k in types)
    proportion = Dict(k => type_counts[k] / n for k in types)
    
    # Calculate bounds for each type
    lb = Dict(k => floor(Int, p * proportion[k]) for k in types)
    ub = Dict(k => ceil(Int, p * proportion[k]) for k in types)
    
    # Initialize solution
    selected = Int[]
    selected_by_type = Dict(k => 0 for k in types)
    
    # Phase 1: Greedy construction with type constraints
    remaining = p
    
    # First, ensure we meet lower bounds for each type
    for k in types
        # We need to select at least lb[k] points of type k
        candidates = findall(i -> data[i][3] == k && !(i in selected), 1:n)
        needed = lb[k]
        
        while needed > 0 && !isempty(candidates)
            best_idx = -1
            best_value = -Inf
            
            for idx in candidates
                # For the first point, use distance from center
                if isempty(selected)
                    center_x = sum(data[i][1] for i in 1:n) / n
                    center_y = sum(data[i][2] for i in 1:n) / n
                    value = sqrt((data[idx][1] - center_x)^2 + (data[idx][2] - center_y)^2)
                else
                    # For subsequent points, evaluate based on objective
                    if objective == "sumsum"
                        value = sum(d[idx, j] for j in selected)
                    else # minmin
                        value = minimum(d[idx, j] for j in selected)
                    end
                end
                
                if value > best_value
                    best_value = value
                    best_idx = idx
                end
            end
            
            if best_idx != -1
                push!(selected, best_idx)
                selected_by_type[k] += 1
                needed -= 1
                remaining -= 1
                filter!(x -> x != best_idx, candidates)
            else
                break
            end
        end
    end
    
    # Phase 2: Fill remaining slots greedily while respecting upper bounds
    while remaining > 0
        best_idx = -1
        best_value = -Inf
        
        for i in 1:n
            if !(i in selected) && selected_by_type[data[i][3]] < ub[data[i][3]]
                # Evaluate contribution to objective
                if isempty(selected)
                    # For the first point, use distance from center
                    center_x = sum(data[j][1] for j in 1:n) / n
                    center_y = sum(data[j][2] for j in 1:n) / n
                    value = sqrt((data[i][1] - center_x)^2 + (data[i][2] - center_y)^2)
                else
                    # For subsequent points
                    if objective == "sumsum"
                        value = sum(d[i, j] for j in selected)
                    else # minmin
                        value = minimum(d[i, j] for j in selected)
                    end
                end
                
                if value > best_value
                    best_value = value
                    best_idx = i
                end
            end
        end
        
        if best_idx == -1
            break  # No valid points left to add
        end
        
        push!(selected, best_idx)
        selected_by_type[data[best_idx][3]] += 1
        remaining -= 1
    end
    
    # Phase 3: Local search improvement
    if local_search && length(selected) >= 2
        improved = true
        iterations = 0
        max_iterations = 100
        
        while improved && iterations < max_iterations
            improved = false
            iterations += 1
            
            # Calculate current objective value
            current_obj = 0
            if objective == "sumsum"
                current_obj = sum(d[s1, s2] for s1 in selected for s2 in selected if s1 < s2)
            else # minmin
                current_obj = minimum(d[s1, s2] for s1 in selected for s2 in selected if s1 < s2)
            end
            
            # Try swapping each selected point with each non-selected point of same type
            for i in selected
                type_i = data[i][3]
                
                for j in 1:n
                    if !(j in selected) && data[j][3] == type_i
                        # Simulate swap
                        new_selected = copy(selected)
                        replace!(new_selected, i => j)
                        
                        # Calculate new objective
                        new_obj = 0
                        if objective == "sumsum"
                            new_obj = sum(d[s1, s2] for s1 in new_selected for s2 in new_selected if s1 < s2)
                        else # minmin
                            new_obj = minimum(d[s1, s2] for s1 in new_selected for s2 in new_selected if s1 < s2)
                        end
                        
                        # If improvement found, make the swap
                        if new_obj > current_obj
                            replace!(selected, i => j)
                            improved = true
                            break
                        end
                    end
                end
                
                if improved
                    break
                end
            end
        end
    end
    
    # Calculate final objective value
    obj_value = 0
    if length(selected) >= 2
        if objective == "sumsum"
            obj_value = sum(d[s1, s2] for s1 in selected for s2 in selected if s1 < s2)
        else # minmin
            obj_value = minimum(d[s1, s2] for s1 in selected for s2 in selected if s1 < s2)
        end
    end
    
    return selected, selected_by_type, obj_value
end

# NEW FUNCTION: Simple p-dispersion algorithm
function pdisp_simple(d, p, N)
    maxdist = 0
    bestpair = (0, 1)
    
    # Find the pair of points with maximum distance
    for i in 1:N
        for j in i+1:N
            if d[i, j] > maxdist
                maxdist = d[i, j]
                bestpair = (i, j)
            end
        end
    end
    
    # Initialize solution with the maximum distance pair
    P = Set([])
    push!(P, bestpair[1])
    push!(P, bestpair[2])
    
    # Iteratively add points to maximize the minimum distance
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

# NEW FUNCTION: Determine the type of a node
function node_type(i, Sk)
    for k in eachindex(Sk)
        if i in Sk[k]
            return k
        end
    end
    println("Node $i not found in Sk")
    return 0  # Return default if not found
end

# NEW FUNCTION: Count how many points of each type are in solution P
function count_k(P, Sk)
    count = zeros(Int, length(Sk))
    for i in P
        k = node_type(i, Sk)
        count[k] += 1
    end
    return count
end

# NEW FUNCTION: Enhanced p-dispersion algorithm with type constraints
function pdisp_2(data, p)
    # Extract data
    n = length(data)  # Number of points
    
    # Calculate distances
    d = zeros(n, n)
    for i in 1:n
        for j in 1:n
            d[i,j] = round(sqrt((data[i][1] - data[j][1])^2 + (data[i][2] - data[j][2])^2))
        end
    end
    
    # Get unique types and create Sk (sets of nodes by type)
    types = sort(unique([point[3] for point in data]))
    Sk = [Int[] for _ in 1:length(types)]
    
    for i in 1:n
        type_idx = findfirst(t -> t == data[i][3], types)
        push!(Sk[type_idx], i)
    end
    
    # Calculate proportions and bounds for each type
    proportion = Dict(k => length(Sk[k]) / n for k in 1:length(types))
    Lk = Dict(k => floor(Int, p * proportion[k]) for k in 1:length(types))
    Uk = Dict(k => ceil(Int, p * proportion[k]) for k in 1:length(types))
    
    # Create distance matrices for each type
    d_by_type = []
    for k in 1:length(types)
        d_k = zeros(length(Sk[k]), length(Sk[k]))
        for i in 1:length(Sk[k])
            for j in 1:length(Sk[k])
                d_k[i, j] = d[Sk[k][i], Sk[k][j]]
            end
        end
        push!(d_by_type, d_k)
    end
    
    # Solve p-dispersion for each type
    p_disp_by_type = []
    for k in 1:length(types)
        p_k = pdisp_simple(d_by_type[k], Lk[k], length(Sk[k]))
        # Convert indices back to original data indices
        p_k_fixed = [Sk[k][i] for i in p_k]
        push!(p_disp_by_type, p_k_fixed)
    end
    
    # Combine solutions
    pdisp_ok = Set(vcat(p_disp_by_type...))
    
    # Adjust solution if needed
    if length(pdisp_ok) != p
        count = count_k(pdisp_ok, Sk)
        
        while length(pdisp_ok) < p
            # Find the node v that maximizes the distance to its closest neighbor in P
            maxdist = 0
            vbest = 0
            
            for v in 1:n
                if v in pdisp_ok
                    continue
                end
                
                k = node_type(v, Sk)
                if k == 0 || count[k] >= Uk[k]
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
                println("PDISP_2 could not find a valid solution that satisfies all constraints")
                break
            end
            
            # Add the node vbest to the set P and update the counts
            k = node_type(vbest, Sk)
            if k > 0  # Ensure valid type
                count[k] += 1
                push!(pdisp_ok, vbest)
            end
        end
    end
    
    # Calculate objective value - minimum distance between any two selected points
    collection = collect(pdisp_ok)
    min_dist = Inf
    
    if length(collection) >= 2
        for i in 1:length(collection)
            for j in i+1:length(collection)
                min_dist = min(min_dist, d[collection[i], collection[j]])
            end
        end
    end
    
    return collection, min_dist
end

# Plot the solution with different shapes for types and colors for selection status
function plot_solution(data, selected_indices; filename="solution_plot.png", title="Multi-Type p-Dispersion Solution")
    # Create a new plot
    p = plot(
        size=(800, 800),
        xlabel="X",
        ylabel="Y",
        title=title,
        legend=:topright,
        framestyle=:box,
        grid=false
    )
    
    # Define marker shapes for different types
    # Use different markers for different types: :circle, :rect, :star5, :diamond, :hexagon
    shapes = Dict(
        1 => :circle,
        2 => :rect,
        3 => :star5,
        4 => :diamond,
        5 => :hexagon
    )
    
    # Define colors for selected/non-selected
    selected_color = :red
    non_selected_color = :blue
    
    # Get unique types
    types = sort(unique([point[3] for point in data]))
    
    # First plot non-selected points (blue)
    for t in types
        # Get indices of non-selected points of type t
        non_selected = [i for i in 1:length(data) if data[i][3] == t && !(i in selected_indices)]
        
        if !isempty(non_selected)
            x_coords = [data[i][1] for i in non_selected]
            y_coords = [data[i][2] for i in non_selected]
            
            scatter!(
                p, 
                x_coords, 
                y_coords, 
                markershape=shapes[t],
                markercolor=non_selected_color, 
                markerstrokecolor=:black,
                markersize=6,
                markerstrokewidth=1,
                label="Type $t (not selected)"
            )
        end
    end
    
    # Then plot selected points (red) so they appear on top
    for t in types
        # Get indices of selected points of type t
        selected = [i for i in selected_indices if data[i][3] == t]
        
        if !isempty(selected)
            x_coords = [data[i][1] for i in selected]
            y_coords = [data[i][2] for i in selected]
            
            scatter!(
                p, 
                x_coords, 
                y_coords, 
                markershape=shapes[t],
                markercolor=selected_color, 
                markerstrokecolor=:black,
                markersize=8,
                markerstrokewidth=1.5,
                label="Type $t (selected)"
            )
        end
    end
    
    # Save the plot to a file
    savefig(p, filename)
    
    return p
end

# Compare multiple solutions
function compare_multiple_solutions(data, solutions_dict; filename="multi_solution_comparison.png")
    # Create a new plot
    p = plot(
        size=(1000, 800),
        xlabel="X",
        ylabel="Y",
        title="Comparison of Multiple Solution Methods",
        legend=:topright,
        framestyle=:box,
        grid=false
    )
    
    # Define marker shapes for different types
    shapes = Dict(
        1 => :circle,
        2 => :rect,
        3 => :star5,
        4 => :diamond,
        5 => :hexagon
    )
    
    # Define colors for different solutions
    solution_colors = Dict(
        "Optimal" => :red,
        "Heuristic" => :blue,
        "Simple" => :green,
        "Type-Constrained" => :purple
    )
    
    # Get unique types
    types = sort(unique([point[3] for point in data]))
    
    # First plot non-selected points (light gray)
    non_selected = Set(1:length(data))
    for (_, indices) in solutions_dict
        non_selected = setdiff(non_selected, indices)
    end
    
    for t in types
        # Get indices of non-selected points of type t
        type_non_selected = [i for i in non_selected if data[i][3] == t]
        
        if !isempty(type_non_selected)
            x_coords = [data[i][1] for i in type_non_selected]
            y_coords = [data[i][2] for i in type_non_selected]
            
            scatter!(
                p, 
                x_coords, 
                y_coords, 
                markershape=shapes[t],
                markercolor=:lightgray, 
                markerstrokecolor=:gray,
                markersize=4,
                markerstrokewidth=0.5,
                label=t==1 ? "Non-selected" : ""  # Only show in legend once
            )
        end
    end
    
    # Plot each solution
    for (method_name, indices) in solutions_dict
        color = solution_colors[method_name]
        
        for t in types
            # Get indices of selected points of type t
            selected = [i for i in indices if data[i][3] == t]
            
            if !isempty(selected)
                x_coords = [data[i][1] for i in selected]
                y_coords = [data[i][2] for i in selected]
                
                scatter!(
                    p, 
                    x_coords, 
                    y_coords, 
                    markershape=shapes[t],
                    markercolor=color, 
                    markerstrokecolor=:black,
                    markersize=8,
                    markerstrokewidth=1.5,
                    label="Type $t ($method_name)"
                )
            end
        end
    end
    
    savefig(p, filename)
    println("Multi-solution comparison plot saved to '$filename'")
    
    return p
end

# Execute the models
function run_models()
    pdp_min_model, multi_pdp_min_model, multi_pdp_sum_model, data, d, lb, ub, K = build_multi_p_dispersion()
    p = 40  # Number of points to select
    
    # Store all solutions for comparison
    solutions = Dict{String, Vector{Int}}()
    
    # Solve Model B: Multi p-Dispersion with minmin objective
    println("\nSolving multi p-dispersion model (multi_pdp_min)...")
    model_b = multi_pdp_min_model()
    set_time_limit_sec(model_b, 2000)
    optimize!(model_b)
    
    if termination_status(model_b) == MOI.OPTIMAL || termination_status(model_b) == MOI.TIME_LIMIT
        y_val_b = value.(model_b[:y])
        y_round = round.(Int, y_val_b)
        z_min_val_b = objective_value(model_b)
        
        println("Solution status: ", termination_status(model_b))
        println("Objective value (z_min): ", z_min_val_b)
        
        selected_indices = findall(x -> x == 1, y_round.data)
        println("Selected points (optimal): ", selected_indices)
        
        # Plot the optimal solution
        plot_solution(data, selected_indices, filename="optimal_solution.png", title="Optimal Solution (MinMin Objective)")
        println("Optimal solution plot saved to 'optimal_solution.png'")
        
        # Store solution
        solutions["Optimal"] = selected_indices
    else
        println("Model B could not be solved to optimality.")
        println("Status: ", termination_status(model_b))
    end
    
    # Run the multi-type p-dispersion heuristic
    println("\nRunning multi-type heuristic solution...")
    heuristic_selected, type_distribution, heuristic_obj = multi_type_pdp_heuristic(data, p, objective="minmin")
    println("Heuristic solution selected points: ", heuristic_selected)
    println("Heuristic objective value: ", heuristic_obj)
    println("Type distribution: ", type_distribution)
    
    # Plot the heuristic solution
    plot_solution(data, heuristic_selected, filename="heuristic_solution.png", title="Heuristic Solution (MinMin Objective)")
    println("Heuristic solution plot saved to 'heuristic_solution.png'")
    
    # Store solution
    solutions["Heuristic"] = heuristic_selected
    
    # Run the simple p-dispersion algorithm (no type constraints)
    println("\nRunning simple p-dispersion algorithm...")
    n = length(data)
    simple_selected = pdisp_simple(d, p, n)
    
    # Calculate objective value for simple solution
    simple_obj = Inf
    for i in simple_selected
        for j in simple_selected
            if i < j
                simple_obj = min(simple_obj, d[i, j])
            end
        end
    end
    
    println("Simple p-dispersion solution: ", simple_selected)
    println("Simple p-dispersion objective value: ", simple_obj)
    
    # Check type distribution of simple solution
    simple_type_count = Dict(k => count(i -> data[i][3] == k, simple_selected) for k in K)
    println("Simple solution type distribution: ", simple_type_count)
    
    # Plot the simple solution
    plot_solution(data, simple_selected, filename="simple_solution.png", title="Simple p-Dispersion Solution (No Type Constraints)")
    println("Simple solution plot saved to 'simple_solution.png'")
    
    # Store solution
    solutions["Simple"] = simple_selected
    
    # Run the type-constrained p-dispersion algorithm (pdisp_2)
    println("\nRunning type-constrained p-dispersion algorithm...")
    type_constrained_selected, type_constrained_obj = pdisp_2(data, p)
    
    println("Type-constrained p-dispersion solution: ", type_constrained_selected)
    println("Type-constrained p-dispersion objective value: ", type_constrained_obj)
    
    # Check type distribution of type-constrained solution
    type_constrained_count = Dict(k => count(i -> data[i][3] == k, type_constrained_selected) for k in K)
    println("Type-constrained solution type distribution: ", type_constrained_count)
    
    # Plot the type-constrained solution
    plot_solution(data, type_constrained_selected, filename="type_constrained_solution.png", title="Type-Constrained p-Dispersion Solution")
    println("Type-constrained solution plot saved to 'type_constrained_solution.png'")
    
    # Store solution
    solutions["Type-Constrained"] = type_constrained_selected
    
    # Create a multi-solution comparison plot
    compare_multiple_solutions(data, solutions)
    
    return data, solutions
end

# Run the models with visualization
run_models()