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
        contador += 1
        path = ARGS[1] * entrada
        println(path)
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
            println("--------------------------------")
            really_cons_old = constraints
            weight_really_before = oldSol.Weight
            weight_before = 1000000000
            valor_previo = 0
            repair_delta = 0
            if factible_old
                #println("Solucion Factible")
                weight_before = solution.Weight
                valor_previo = weight_before
                repair_delta = 0
            else
                before_repair = Dates.now()
                oldSol = repair_solution1(oldSol, constraints, targets_lower, targets_upper, remove, add)
                factible, constraints, remove, add = isFactible4(oldSol, targets_lower, targets_upper)
                if factible
                    weight_before = oldSol.Weight
                    valor_previo = weight_before
                end
                after_repair = Dates.now()
                repair_delta_time = after_repair - before_repair
                repair_delta = round(repair_delta_time, Millisecond)
            end
            println("Reparada")

            gap_repaired = (1 - (weight_exac / weight_before)) * 100

            before_simple_bf = Dates.now()
            solution_moved_bf = simple_bu_improve(oldSol, targets_lower, targets_upper, :bf)
            after_simple_bf = Dates.now()
            simple_delta_time_bf = after_simple_bf - before_simple_bf
            simple_delta_bf = round(simple_delta_time_bf, Millisecond)
            weight_move_bf = solution_moved_bf.Weight

            abs_improve_moved_bf = weight_before - weight_move_bf
            rel_improve_moved_bf = (1 - (weight_move_bf / weight_before)) * 100
            gap_moved_bf = (1 - (weight_exac / weight_move_bf)) * 100


            before_interchange_bf = Dates.now()
            solution_interchanged_bf = interchange_bu_improve(oldSol, targets_lower, targets_upper, :bf)
            after_interchange_bf = Dates.now()
            interchange_delta_time_bf = after_interchange_bf - before_interchange_bf
            interchange_delta_bf = round(interchange_delta_time_bf, Millisecond)
            weight_interchange_bf = solution_interchanged_bf.Weight

            abs_improve_interchanged_bf = weight_before - weight_interchange_bf
            rel_improve_interchanged_bf = (1 - (weight_interchange_bf / weight_before)) * 100
            gap_interchanged_bf = (1 - (weight_exac / weight_interchange_bf)) * 100


            before_deactivate_ff = Dates.now()
            solution_deactivated_ff = deactivate_center_improve(oldSol, targets_lower, targets_upper)
            after_deactivate_ff = Dates.now()
            deactivate_delta_time_ff = after_deactivate_ff - before_deactivate_ff
            deactivate_delta_ff = round(deactivate_delta_time_ff, Millisecond)
            weight_deactivate_ff = solution_deactivated_ff.Weight
          
            abs_improve_deactivated_ff = weight_before - weight_deactivate_ff
            rel_improve_deactivated_ff = (1 - (weight_deactivate_ff / weight_before)) * 100
            gap_deactivated_ff = (1 - (weight_exac / weight_deactivate_ff)) * 100


            before_simple_ff = Dates.now()
            solution_moved_ff = simple_bu_improve(oldSol, targets_lower, targets_upper, :ff)
            after_simple_ff = Dates.now()
            simple_delta_time_ff = after_simple_ff - before_simple_ff
            simple_delta_ff = round(simple_delta_time_ff, Millisecond)
            weight_move_ff = solution_moved_ff.Weight

            abs_improve_moved_ff = weight_before - weight_move_ff
            rel_improve_moved_ff = (1 - (weight_move_ff / weight_before)) * 100
            gap_moved_ff = (1 - (weight_exac / weight_move_ff)) * 100

            before_interchange_ff = Dates.now()
            solution_interchanged_ff = interchange_bu_improve(oldSol, targets_lower, targets_upper, :ff)
            after_interchange_ff = Dates.now()
            interchange_delta_time_ff = after_interchange_ff - before_interchange_ff
            interchange_delta_ff = round(interchange_delta_time_ff, Millisecond)
            weight_interchange_ff = solution_interchanged_ff.Weight
        
            abs_improve_interchanged_ff = weight_before - weight_interchange_ff
            rel_improve_interchanged_ff = (1 - (weight_interchange_ff / weight_before)) * 100
            gap_interchanged_ff = (1 - (weight_exac / weight_interchange_ff)) * 100

            
            #=
            str_path = "out\\solutions\\625_78_32\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * location_method * "_" * alloc_method * "_ls"
            plot_str_path = str_path * ".png"
            solution_str_path = str_path * ".jld2"
            solution = Types.Solution(instancia, sol_ls.X, sol_ls.Y, final_weight, total_time)

            Types.plot_solution(solution, plot_str_path)
            Types.write_solution(solution, solution_str_path)
            =#
            row_cons = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P"=> P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli, 
            "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "Value" =>weight_really_before, "TotalTimeCons"=> time_cons, "Factible" => factible_old, 
            "Add" => to_add, "Remove" => to_remove, "Optim" => weight_exac, "TimeOptim" => time_exac)
            push!(arr_constructive, row_cons)

            row_ls_moved_bf = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P"=> P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli, 
            "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" =>weight_really_before, "TotalTimeCons"=> time_cons, "Factible" => factible_old, 
            "Add" => to_add, "Remove" => to_remove, "TimeRepair" => repair_delta.value, "ValueRepair"=>weight_before, "Move"=> "Simple", "Strategy" => "BF", "AbsImproveMove" => abs_improve_moved_bf,
            "RelImproveMove" => rel_improve_moved_bf, "GapRepair" => gap_repaired, "GapMove" => gap_moved_bf, "TimeMove" => simple_delta_bf.value, "Optim" => weight_exac, "TimeOptim" => time_exac, "ValueLS" => weight_move_bf)
            
            row_ls_interchanged_bf = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P"=> P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli, 
            "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" =>weight_really_before, "TotalTimeCons"=> time_cons, "Factible" => factible_old, 
            "Add" => to_add, "Remove" => to_remove, "TimeRepair" => repair_delta.value, "ValueRepair"=>weight_before, "Move"=> "Interchange", "Strategy" => "BF", "AbsImproveMove" => abs_improve_interchanged_bf,
            "RelImproveMove" => rel_improve_interchanged_bf, "GapRepair" => gap_repaired, "GapMove" => gap_interchanged_bf, "TimeMove" => interchange_delta_bf.value, "Optim" => weight_exac, "TimeOptim" => time_exac, "ValueLS" => weight_interchange_bf)


            row_ls_moved_ff = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P"=> P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli, 
            "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" =>weight_really_before, "TotalTimeCons"=> time_cons, "Factible" => factible_old, 
            "Add" => to_add, "Remove" => to_remove, "TimeRepair" => repair_delta.value, "ValueRepair"=>weight_before, "Move"=> "Simple", "Strategy" => "FF", "AbsImproveMove" => abs_improve_moved_ff,
            "RelImproveMove" => rel_improve_moved_ff, "GapRepair" => gap_repaired, "GapMove" => gap_moved_ff, "TimeMove" => simple_delta_ff.value, "Optim" => weight_exac, "TimeOptim" => time_exac, "ValueLS" => weight_move_ff)
            
            row_ls_interchanged_ff = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P"=> P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli, 
            "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" =>weight_really_before, "TotalTimeCons"=> time_cons, "Factible" => factible_old, 
            "Add" => to_add, "Remove" => to_remove, "TimeRepair" => repair_delta.value, "ValueRepair"=>weight_before, "Move"=> "Interchange", "Strategy" => "FF", "AbsImproveMove" => abs_improve_interchanged_ff,
            "RelImproveMove" => rel_improve_interchanged_ff, "GapRepair" => gap_repaired, "GapMove" => gap_interchanged_ff, "TimeMove" => interchange_delta_ff.value, "Optim" => weight_exac, "TimeOptim" => time_exac, "ValueLS" => weight_interchange_ff)

            row_ls_deactivated_ff = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P"=> P, "Loc" => location_method, "Alloc" => alloc_method, "TimeLoc" => delta_loc_milli, 
            "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValueCons" =>weight_really_before, "TotalTimeCons"=> time_cons, "Factible" => factible_old, 
            "Add" => to_add, "Remove" => to_remove, "TimeRepair" => repair_delta.value, "ValueRepair"=>weight_before, "Move"=> "Deactivate", "Strategy"=> "FF", "AbsImproveMove" => abs_improve_deactivated_ff,
            "RelImproveMove" => rel_improve_deactivated_ff, "GapRepair" => gap_repaired, "GapMove" => gap_deactivated_ff, "TimeMove" => deactivate_delta_ff.value, "Optim" => weight_exac, "TimeOptim" => time_exac, "ValueLS" => weight_deactivate_ff)
            
            push!(arr_ls, row_ls_interchanged_ff)
            push!(arr_ls, row_ls_deactivated_ff)
            push!(arr_ls, row_ls_moved_ff)
            push!(arr_ls, row_ls_moved_bf)
            push!(arr_ls, row_ls_interchanged_bf)
            end
        end
    df1 = vcat(DataFrame.(arr_constructive)...)
    df2 = vcat(DataFrame.(arr_ls)...)
    CSV.write("df_cons_625.csv", df1)
    CSV.write("df_ls_625.csv", df2)
end

main_test()
