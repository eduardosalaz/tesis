using JuMP
using Gurobi
using LinearAlgebra
using CPLEX
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
        
        @variable(model, y[I], Bin)      # 1 if point i is selected
        @variable(model, u >= 0)         # Minimum distance between any two selected points
        @variable(model, z[i=I,j=I; i<j], Bin)  # 1 if both points i and j are selected
        @variable(model, w[K], Int)      # Number of points of type k selected
        
        # Objective: Maximize the minimum distance
        @objective(model, Max, u)
        
        # Constraint: z[i,j] = 1 if and only if both i and j are selected
        for i in I
            for j in I
                if i < j
                    @constraint(model, z[i,j] <= y[i])
                    @constraint(model, z[i,j] <= y[j])
                    @constraint(model, z[i,j] >= y[i] + y[j] - 1)
                end
            end
        end
        
        # Minimum distance constraints with tighter linearization
        for i in I
            for j in I
                if i < j
                    @constraint(model, u <= d[i,j] + Dmax * (1 - z[i,j]))
                end
            end
        end
        
        # Remaining constraints...
        @constraint(model, sum(y) == p)
        
        for k in K
            @constraint(model, w[k] == sum(y[i] for i in I if data[i][3] == k))
            @constraint(model, w[k] >= lb[k])
            @constraint(model, w[k] <= ub[k])
        end
        

        
        
        return model
    end

    function multi_pdp_min_alternative()
        print("hola")
        model = Model(Gurobi.Optimizer)
        
        # Variables for point selection
        @variable(model, y[I], Bin)  # 1 if point i is selected
        @variable(model, w[K], Int)  # Number of points of type k selected
        
        # Get all unique distance values
        possible_distances = sort(unique([round(sqrt((data[i][1] - data[j][1])^2 + 
                                                   (data[i][2] - data[j][2])^2)) 
                                         for i in I for j in I if i < j]))
        
        # Binary variable for each possible distance - equals 1 if this is the minimum distance
        @variable(model, z[possible_distances], Bin)
        
        # Exactly one distance value is the minimum
        @constraint(model, sum(z[dist] for dist in possible_distances) == 1)
        
        # Link variables: if z[dist]=1, then no pair of points closer than dist can be selected
        for dist in possible_distances
            for i in I
                for j in I
                    if i < j && d[i,j] < dist
                        @constraint(model, y[i] + y[j] <= 1 + (1 - z[dist]))
                    end
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
        
        # Objective: maximize the selected minimum distance
        @objective(model, Max, sum(dist * z[dist] for dist in possible_distances))
        
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
    
    return pdp_min, multi_pdp_min, multi_pdp_min_alternative, data
end

# Execute the models
function run_models()
    pdp_min_model, multi_pdp_min_model, multi_pdp_min_alternative, data = build_multi_p_dispersion()
    #=
    # Solve Model C: p-Dispersion with minmin objective
    println("\nSolving basic p-dispersion model (pdp_min)...")
    model_c = pdp_min_model()
    #set_silent(model_c)
    set_time_limit_sec(model_c, 2000)
    optimize!(model_c)
    
    if termination_status(model_c) == MOI.OPTIMAL
        y_val_c = value.(model_c[:y])
        z_min_val_c = objective_value(model_c)
        
        println("Solution status: ", termination_status(model_c))
        println("Objective value (z_min): ", z_min_val_c)
        #println("Selected points (y): ", findall(x -> x > 0.5, y_val_c))
    else
        println("Model C could not be solved to optimality.")
        println("Status: ", termination_status(model_c))
    end
    =#
    # Solve Model B: Multi p-Dispersion with minmin objective
    println("\nSolving multi p-dispersion model (multi_pdp_min)...")
    model_b = multi_pdp_min_model()
    #set_silent(model_b)
    set_time_limit_sec(model_b, 2000)
    optimize!(model_b)
    
    if termination_status(model_b) == MOI.OPTIMAL
        y_val_b = value.(model_b[:y])
        w_val_b = value.(model_b[:w])
        z_min_val_b = objective_value(model_b)
        
        println("Solution status: ", termination_status(model_b))
        println("Objective value (z_min): ", z_min_val_b)
        #println("Selected points (y): ", findall(x -> x > 0.5, y_val_b))
        #println("Points by type (w): ", Dict(k => w_val_b[k] for k in keys(w_val_b)))
        
        # Show actual points selected
        #=
        selected_indices = findall(x -> x > 0.5, y_val_b)
        selected_points = [data[i] for i in selected_indices]
        println("Selected point coordinates and types:")
        for (idx, point) in zip(selected_indices, selected_points)
            println("P$(lpad(idx, 3, "0")): ($(point[1]), $(point[2])), type $(point[3])")
        end
        =#
    else
        println("Model B could not be solved to optimality.")
        println("Status: ", termination_status(model_b))
    end
end

# Run the models
run_models()