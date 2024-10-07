include("constructive.jl")
include("ls.jl")
using DataFrames
using CSV
using LinearAlgebra
using DelimitedFiles
function main_test()
    entradas = readdir(ARGS[1])
    arr_constructive = []
    arr_ls = []
    best_sols_gurobi_625 = Dict(
        "1" => 464274,
        "2" => 475106,
        "3" => 451907,
        "4" => 460877,
        "5" => 469190,
        "6" => 491611,
        "7" => 454550,
        "8" => 462113,
        "9" => 488971,
        "10" => 471810,
        "11" => 468195,
        "12" => 461823,
        "13" => 468480,
        "14" => 459301,
        "15" => 483950,
        "16" => 462365,
        "17" => 462745,
        "18" => 455692,
        "19" => 461881,
        "20" => 465070
    )

    best_bounds_gurobi_625 = Dict(
        "1" => 462154.445,
        "2" => 474267.150,
        "3" => 451740.158,
        "4" => 460520.424,
        "5" => 468298.942,
        "6" => 490305.188,
        "7" => 453366.406,
        "8" => 460970.188,
        "9" => 487836.913,
        "10" => 467610.276,
        "11" => 466413.292,
        "12" => 461040.141,
        "13" => 468091.441,
        "14" => 459097.840,
        "15" => 479672.397,
        "16" => 461967.277,
        "17" => 462126.537,
        "18" => 454556.375,
        "19" => 458774.946,
        "20" => 464547.884
    )

    time_gurobi_625 = Dict(
        "1" => 1800,
        "2" => 1800,
        "3" => 1800,
        "4" => 1800,
        "5" => 1800,
        "6" => 1800,
        "7" => 1800,
        "8" => 1800,
        "9" => 1800,
        "10" => 1800,
        "11" => 1800,
        "12" => 1800,
        "13" => 1800,
        "14" => 1800,
        "15" => 1800,
        "16" => 1800,
        "17" => 1800,
        "18" => 1800,
        "19" => 1800,
        "20" => 1800
    )


    best_sols_gurobi_1250 = Dict(
        "1" => 949005,
        "2" => 975979,
        "3" => 988670,
        "4" => 974186,
        "5" => 958455,
        "6" => 964088,
        "7" => 968567,
        "8" => 986315,
        "9" => 951731,
        "10" => 999473,
        "11" => 961532,
        "12" => 942864,
        "13" => 1012710,
        "14" => 987762,
        "15" => 961675,
        "16" => 990828,
        "17" => 992812,
        "18" => 980761,
        "19" => 973562,
        "20" => 972605
    )

    best_bounds_gurobi_1250 = Dict(
        "1" => 944758.818,
        "2" => 966579.707,
        "3" => 982010.310,
        "4" => 971268.401,
        "5" => 957947.388,
        "6" => 961925.879,
        "7" => 962977.434,
        "8" => 983038.375,
        "9" => 947306.650,
        "10" => 991050.288,
        "11" => 954783.337,
        "12" => 940956.264,
        "13" => 1009464.27,
        "14" => 985021.255,
        "15" => 957101.370,
        "16" => 983915.691,
        "17" => 991539.114,
        "18" => 974474.121,
        "19" => 970867.408,
        "20" => 968213.682
    )

    time_gurobi_1250 = Dict(
        "1" => 1800,
        "2" => 1800,
        "3" => 1800,
        "4" => 1800,
        "5" => 1800,
        "6" => 1800,
        "7" => 1800,
        "8" => 1800,
        "9" => 1800,
        "10" => 1800,
        "11" => 1800,
        "12" => 1800,
        "13" => 1800,
        "14" => 1800,
        "15" => 1800,
        "16" => 1800,
        "17" => 1800,
        "18" => 1800,
        "19" => 1800,
        "20" => 1800
    )

    for entrada in entradas
        path = ARGS[1] * entrada
        instancia = read_instance(path)
        instance = instancia
        pattern = Regex("[t][_]\\d{1,3}")
        index = findfirst(pattern, path)
        almost_number = path[index]
        _, number = split(almost_number, "_")
        println(number)
        id = number
        #sol_exac_path = "out\\solutions\\625_78_32\\sol_" * id * "_625_78_32.jld2"
        #sol_exac = read_solution(sol_exac_path)

        weight_exac = best_sols_gurobi_625[number]
        bound_exac = best_bounds_gurobi_625[number]
        time_exac = time_gurobi_625[number]

        B = instancia.B
        S = instancia.S
        P = instancia.P
        X = Matrix{Int64}(undef, S, B)
        Y = Vector{Int64}(undef, S)
        D = instancia.D
        Weight = 0
        

        location_methods = ["pdisp"]
        alloc_methods = ["queue"]
        for location_method in location_methods, alloc_method in alloc_methods
            time_loc = 1000
            if location_method == "relax"
                Y, time_loc = localize_facs(instancia, location_method)
                if time_loc == 0
                    time_loc = 1000
                else
                    time_loc = time_loc * 1000
                end
            else
                Y, time_loc = localize_facs(instancia, location_method)
            end
            if location_method ≠ "relax"
                Y_bool = zeros(Int, instancia.S)
                for idx in Y
                    Y_bool[idx] = 1
                end
            else
                Y_bool = Y
            end
            delta_loc_milli = time_loc
            #println("Y done with time: ", delta_loc_milli)
            before_alloc = Dates.now()
            unassigned_count = 1
            if alloc_method == "opp"
                #println("OPP COST")
                X = oppCostAssignment(Y_bool, instancia)
                unassigned_count = 0
            elseif alloc_method == "queue"
                X, unassigned_count = oppCostQueue(Y_bool, instancia)
            end
            after_alloc = Dates.now()
            delta_alloc = after_alloc - before_alloc
            delta_alloc_milli = round(delta_alloc, Millisecond)
            delta = delta_loc_milli + delta_alloc_milli.value
            println("X asignada")
            #micros = round(delta, Millisecond)
            time_cons = delta
            Weight = dot(X, D)
            if !isdir("pruebas_new")
                mkdir("pruebas_new")
            end

            for col in 1:B
                if all(==(0), X[:, col])
                    println("uuhhh all 0s in $col 1")
                end
            end

            str_path = "out\\solutions\\625_78_32\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * location_method * "_" * alloc_method * "_new_final2"
            plot_str_path = str_path * ".png"
            solution_str_path = str_path * ".jld2"
            solution = Types.Solution(instancia, X, Y_bool, Weight, time_cons)

            #Types.plot_solution(solution, plot_str_path)
            Types.write_solution(solution, solution_str_path)

            oldSol = solution
            oldTime = oldSol.Time
            targets_lower, targets_upper = calculate_targets(instance)
            targets_lower_op, targets_upper_op = calculate_targets_optimized(instance)


            factible_old, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper)
            to_remove = length(remove)
            to_add = length(add)
            really_cons_old = constraints
            weight_really_before = oldSol.Weight
            weight_before = 1000000000
            repair_delta = 0

            instance = solution.Instance
            factible_after_repair = false
            targets_lower, targets_upper = calculate_targets(instance)
            factible, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper)
            #println(remove)
            #println(add)
            repaired = oldSol
            original_weight = 10000000000000
            weight_before = 0
            total_time = 0
            before_repair = Dates.now()
            repair_algorithm = 1
            if factible
                println("Factible")
                original_weight = solution.Weight
                weight_before = original_weight
                factible_after_repair = true
            else
                println("Reparando")
                #println(isFactible(oldSol, true))
                #writedlm("aver.txt", oldSol.X)
                repaired_1 = repair_solution1(oldSol, constraints, targets_lower, targets_upper, remove, add)
                fac_repaired_1, cons = isFactible(repaired_1, false)
                if !fac_repaired_1
                    repaired_2 = repair_solution2(oldSol, constraints, targets_lower, targets_upper, remove, add)
                    fac_repaired_2, cons = isFactible(repaired_2, false)
                    if fac_repaired_2
                        repair_algorithm = 2
                        factible_after_repair = true
                        repaired = repaired_2
                    end
                else
                    repaired = repaired_1
                    factible_after_repair = true
                end
                if factible_after_repair
                    original_weight = repaired.Weight
                    weight_before = repaired.Weight
                    println("Reparada")
                end
            end
            after_repair = Dates.now()
            repair_delta_time = after_repair - before_repair
            repair_delta = round(repair_delta_time, Millisecond)
            delta_ls_value = 0
            weight_after = 10000000
            counter = 0
            counter_improve_simple = 0
            counter_improve_interchange = 0
            counter_improve_deactivate = 0
            if !factible_after_repair
                @error "INFACTIBLE"
            else
                oldSol = repaired
                oldSol2 = repaired
                improvement = true
                before_ls = Dates.now()
                improvement = true
                loop = 0
                while improvement
                    improvement = false  # Reset the flag at the start of each loop iteration
                    prev_weight = oldSol.Weight

                    # Array to keep track of individual improvements
                    improvements = Bool[]
                    # First improvement function
                    sol_moved_bu = simple_bu_improve_optimized(oldSol, targets_lower_op, targets_upper_op, :ff)
                    new_weight_moved = sol_moved_bu.Weight
                    push!(improvements, new_weight_moved < prev_weight)
                    if improvements[end]
                        prev_weight = new_weight_moved
                        #println("En el loop loop el movimiento simple mejora con un new_weight_moved")
                        oldSol = sol_moved_bu  # Update oldSol if there was an improvement
                    end
                    #println(isFactible(sol_moved_bu, true))
                    # Second improvement function
                    sol_interchanged_bu = interchange_bu_improve_optimized(oldSol, targets_lower_op, targets_upper_op, :ff)
                    new_weight_moved = sol_interchanged_bu.Weight
                    push!(improvements, new_weight_moved < prev_weight)
                    if improvements[end]
                        prev_weight = new_weight_moved
                        #println("En el loop loop el movimiento intercambio mejora con un new_weight_moved")
                        oldSol = sol_interchanged_bu  # Update oldSol if there was an improvement
                    end
                    #println(isFactible(sol_interchanged_bu, true))

                    # Third improvement function
                    sol_deactivated_center = deactivate_center_improve(oldSol, targets_lower, targets_upper)
                    new_weight_moved = sol_deactivated_center.Weight
                    push!(improvements, new_weight_moved < prev_weight)
                    if improvements[end]
                        prev_weight = new_weight_moved
                        #println("En el loop $loop el movimiento desactivar mejora con un $new_weight_moved")
                        oldSol = sol_deactivated_center  # Update oldSol if there was an improvement
                    end
                    counter += 1
                    #println(isFactible(sol_deactivated_center, true))

                    # Check for any improvements

                    improvement = any(improvements)
                    #println(isFactible(oldSol, true))
                end
                after_ls = Dates.now()
                delta_ls = after_ls - before_ls
                delta_ls_millis = round(delta_ls, Millisecond)
                delta_ls_value = delta_ls_millis.value
                weight_after = oldSol.Weight
                str_path = "out\\solutions\\625_78_32\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * location_method * "_" * alloc_method * "_ls_ff_new_final2"
                plot_str_path = str_path * ".png"
                solution_str_path = str_path * ".jld2"
                total_time = time_cons + repair_delta.value + delta_ls_value
                solution = Types.Solution(instancia, oldSol.X, oldSol.Y, weight_after, total_time)
                #Types.plot_solution(solution, plot_str_path)
                Types.write_solution(solution, solution_str_path)
            end
            diff_repaired = abs(weight_before - weight_exac) / abs(weight_before)
            gap_repaired = abs(weight_before - bound_exac) / abs(weight_before)

            diff_improved = abs(weight_after - weight_exac) / abs(weight_after)
            gap_improved = abs(weight_after - bound_exac) / abs(weight_after)
            abs_improved = weight_before - weight_after
            rel_improved = ((abs_improved) / weight_after) * 100
            #time_exac = time_exac
            beat_gurobi = false
            if weight_after < weight_exac
                beat_gurobi = true
            end

            row_cons = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli,
                "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" => weight_really_before, "TotalTimeCons" => time_cons, "Factible" => factible_old,
                "Add" => to_add, "Remove" => to_remove, "ValueOptim" => weight_exac, "TimeOptim" => time_exac, "BestBound" => bound_exac, "Unassigned" => unassigned_count)
            push!(arr_constructive, row_cons)

            if factible_after_repair
                row_ls_factible = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli,
                    "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" => weight_really_before, "TimeCons" => time_cons, "Factible" => factible_old,
                    "Add" => to_add, "Remove" => to_remove, "TimeRepair" => repair_delta.value, "ValueRepair" => weight_before, "ImproveAbsLS" => abs_improved,
                    "ImproveRelLS" => rel_improved, "GapRepair" => gap_repaired, "GapLS" => gap_improved, "DiffRepair" => diff_repaired, "DiffLS" => diff_improved, "TimeLS" => delta_ls_value, "ValueOptim" => weight_exac, "TimeOptim" => time_exac, "TimeTotal" => total_time,
                    "CyclesLS" => counter, "SimpleImprovedTimes" => counter_improve_simple, "InterchangeImprovedTimes" => counter_improve_interchange, "DeactivateImprovedTimes" => counter_improve_deactivate, "ValueLS" => weight_after,
                    "RepairAlgorithm" => repair_algorithm, "BeatGurobi" => beat_gurobi, "BestBound" => bound_exac, "Repaired" => factible_after_repair)
                push!(arr_ls, row_ls_factible)
            else
                row_ls_infactible = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli,
                    "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" => weight_really_before, "TimeCons" => time_cons, "Factible" => factible_old,
                    "Add" => to_add, "Remove" => to_remove, "TimeRepair" => repair_delta.value, "ValueRepair" => 10000000, "ImproveAbsLS" => 0, "TimeTotal" => time_cons + repair_delta.value,
                    "ImproveRelLS" => 0, "GapRepair" => 100, "GapLS" => 100, "DiffRepair" => 100, "DiffLS" => 100, "TimeLS" => 0, "ValueOptim" => weight_exac, "TimeOptim" => time_exac, "ValueLS" => weight_after,
                    "CyclesLS" => 0, "SimpleImprovedTimes" => 0, "InterchangeImprovedTimes" => 0, "DeactivateImprovedTimes" => 0, "RepairAlgorithm" => 0, "BeatGurobi" => beat_gurobi, "BestBound" => bound_exac, "Repaired" => factible_after_repair)
                push!(arr_ls, row_ls_infactible)
            end
        end
        df1 = vcat(DataFrame.(arr_constructive)...)
        df2 = vcat(DataFrame.(arr_ls)...)
        df1 = df1[:, sortperm(names(df1))]
        df2 = df2[:, sortperm(names(df2))]
        CSV.write("new2/df_cons_625_pdisp_queue_ff_new2.csv", df1)
        CSV.write("new2/df_ls_625_pdisp_queue_ff_new2.csv", df2)
    end
end

main_test()
println("\a")
