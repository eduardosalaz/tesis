include("constructive.jl")
include("ls.jl")
using DataFrames
using CSV
function main_test()
    entradas = readdir(ARGS[1])
    arr_constructive = []
    arr_ls = []
    contador = 1
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
        sol_exac_path = "out\\solutions\\625_78_32\\sol_" * id * "_625_78_32.jld2"
        sol_exac = read_solution(sol_exac_path)
        weight_exac = sol_exac.Weight
        time_exac = sol_exac.Time

        B = instancia.B
        S = instancia.S
        P = instancia.P
        X = Matrix{Int64}[]
        Y = Vector{Int64}[]
        D = instancia.D
        Weight = 0

        location_methods = ["relax", "pdisp"]
        alloc_methods = ["opp"]
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
            if alloc_method == "opp"
                #println("OPP COST")
                X = oppCostAssignment(Y_bool, instancia)
            end
            after_alloc = Dates.now()
            delta_alloc = after_alloc - before_alloc
            delta_alloc_milli = round(delta_alloc, Millisecond)
            delta = delta_loc_milli + delta_alloc_milli.value
            println("X asignada")
            #micros = round(delta, Millisecond)
            time_cons = delta
            indices = findall(x -> x == 1, X)
            for indice in indices
                Weight += D[indice]
            end
            if !isdir("pruebas")
                mkdir("pruebas")
            end

            str_path = "out\\solutions\\625_78_32\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * location_method * "_" * alloc_method
            plot_str_path = str_path * ".png"
            solution_str_path = str_path * ".jld2"
            solution = Types.Solution(instancia, X, Y_bool, Weight, time_cons)

            Types.plot_solution(solution, plot_str_path)
            Types.write_solution(solution, solution_str_path)

            oldSol = solution
            oldTime = oldSol.Time
            targets_lower, targets_upper = calculate_targets(instance)


            factible_old, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper)
            to_remove = length(remove)
            to_add = length(add)
            really_cons_old = constraints
            weight_really_before = oldSol.Weight
            weight_before = 1000000000
            valor_previo = 0
            repair_delta = 0

            instance = solution.Instance
            factible_after_repair = false
            targets_lower, targets_upper = calculate_targets(instance)
            factible, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper)
            repaired = oldSol
            original_weight = 10000000000000
            weight_before = 0
            before_repair = Dates.now()
            if factible
                println("Factible")
                original_weight = solution.Weight
                weight_before = original_weight
                factible_after_repair = true
            else
                println("Reparando")
                repaired_1 = repair_solution1(oldSol, constraints, targets_lower, targets_upper, remove, add)
                fac_repaired_1, cons = isFactible(repaired_1, false)
                if !fac_repaired_1
                    repaired_2 = repair_solution2(oldSol, constraints, targets_lower, targets_upper, remove, add)
                    fac_repaired_2, cons = isFactible(repaired_2, false)
                    if fac_repaired_2
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
                while improvement
                    counter += 1
                    improvement = false
                    prev_weight = oldSol.Weight
                    sol_moved_bu = simple_bu_improve(oldSol, targets_lower, targets_upper, :bf)
                    new_weight_moved = sol_moved_bu.Weight
                    if new_weight_moved < prev_weight
                        improvement = true
                        prev_weight = new_weight_moved
                        counter_improve_simple += 1
                    end
                    sol_interchanged_bu = interchange_bu_improve(sol_moved_bu, targets_lower, targets_upper, :bf)
                    new_weight_moved = sol_interchanged_bu.Weight
                    if new_weight_moved < prev_weight
                        improvement = true
                        prev_weight = new_weight_moved
                        counter_improve_interchange += 1
                    end
                    sol_deactivated_center = deactivate_center_improve(sol_interchanged_bu, targets_lower, targets_upper)
                    new_weight_moved = sol_deactivated_center.Weight
                    if new_weight_moved < prev_weight
                        improvement = true
                        prev_weight = new_weight_moved
                        counter_improve_deactivate += 1
                    end
                    oldSol = sol_deactivated_center
                end
                after_ls = Dates.now()
                delta_ls = after_ls - before_ls
                delta_ls_millis = round(delta_ls, Millisecond)
                delta_ls_value = delta_ls_millis.value
                weight_after = oldSol.Weight
                str_path = "out\\solutions\\625_78_32\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * location_method * "_" * alloc_method * "_lsbf"
                plot_str_path = str_path * ".png"
                solution_str_path = str_path * ".jld2"
                total_time = time_cons + repair_delta.value + delta_ls_value
                solution = Types.Solution(instancia, oldSol.X, oldSol.Y, weight_after, total_time)
                Types.plot_solution(solution, plot_str_path)
                Types.write_solution(solution, solution_str_path)
            end
            gap_repaired = (1 - (weight_exac / weight_before)) * 100
            gap_improved = (1 - (weight_exac / weight_after)) * 100
            abs_improved = weight_before - weight_after
            rel_improved = ((abs_improved) / weight_after) * 100
            
            counter = counter - 1
            row_cons = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli,
                "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "Value" => weight_really_before, "TotalTimeCons" => time_cons, "Factible" => factible_old,
                "Add" => to_add, "Remove" => to_remove, "ValueOptim" => weight_exac, "TimeOptim" => time_exac*1000)
            push!(arr_constructive, row_cons)

            if factible_after_repair
                row_ls_factible = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli,
                    "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" => weight_really_before, "TimeCons" => time_cons, "Factible" => factible_old,
                    "Add" => to_add, "Remove" => to_remove, "TimeRepair" => repair_delta.value, "ValueRepair" => weight_before, "Strategy" => "BF", "ImproveAbsLS" => abs_improved,
                    "ImproveRelLS" => rel_improved, "GapRepair" => gap_repaired, "GapLS" => gap_improved, "TimeLS" => delta_ls_value, "ValueOptim" => weight_exac, "TimeOptim" => time_exac*1000,
                    "CyclesLS" => counter, "SimpleImprovedTimes" => counter_improve_simple, "InterchangeImprovedTimes" => counter_improve_interchange, "DeactivateImprovedTimes" => counter_improve_deactivate)
                push!(arr_ls, row_ls_factible)
            else
                row_ls_infactible = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli,
                    "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" => weight_really_before, "TimeCons" => time_cons, "Factible" => factible_old,
                    "Add" => to_add, "Remove" => to_remove, "TimeRepair" => repair_delta.value, "ValueRepair" => 10000000, "Strategy" => "BF", "ImproveAbsLS" => 0,
                    "ImproveRelLS" => 0, "GapRepair" => 0, "GapLS" => 0, "TimeLS" => 0, "ValueOptim" => weight_exac, "TimeOptim" => time_exac*1000,
                    "CyclesLS" => 0, "SimpleImprovedTimes" => 0, "InterchangeImprovedTimes" => 0, "DeactivateImprovedTimes" => 0)
                push!(arr_ls, row_ls_infactible)
            end
        end
        df1 = vcat(DataFrame.(arr_constructive)...)
        df2 = vcat(DataFrame.(arr_ls)...)
        df1 = df1[:, sortperm(names(df1))]
        df2 = df2[:, sortperm(names(df2))]
        CSV.write("df_cons_625.csv", df1)
        CSV.write("df_ls_625_bf.csv", df2)
        println("\a")
    end
end

main_test()
