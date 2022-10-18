include("constructive.jl")
include("local.jl")
using DataFrames
using CSV
function main_test()
    entradas = readdir(ARGS[1])
    arr_constructive = []
    arr_ls = []
    contador = 0
    for entrada in entradas
        path = ARGS[1] * entrada
        contador += 1
        if contador > 10
            break
        end
        #println(path)
        instancia = read_instance(path)

        pattern = Regex("[t][_]\\d{1,3}")
        index = findfirst(pattern, path)
        almost_number = path[index]
        _, number = split(almost_number, "_")
        println(number)
        id = number
        sol_exac_path = "out\\solutions\\800_150_50\\sol_" * id * "_800_150_50.jld2"
        sol_exac = read_solution(sol_exac_path)
        weight_exac = sol_exac.Weight
        time_exac = sol_exac.Time * 1000

        B = instancia.B
        S = instancia.S
        P = instancia.P
        X = Matrix{Int64}[]
        Y = Vector{Int64}[]
        D = instancia.D
        Weight = 0

        init_methods = ["relax", "pdisp", "random"]
        alloc_methods = ["naive", "opp"]
        for init_method in init_methods, alloc_method in alloc_methods
            before_init = Dates.now()
            time_init = 1
            if init_method == "relax"   
                Y, time_init = localize_facs(instancia, init_method)
                if time_init == 0
                    time_init = 1
                end
                before_init = Dates.now()  
            else
                Y = localize_facs(instancia, init_method)
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
            delta_init = after_init - before_init
            delta_init_milli = round(delta_init, Millisecond)
            if init_method == "relax"
                delta_init_milli = time_init * 1000
            end
            #println("Y done with time: ", delta_init_milli)
            if alloc_method == "naive"
                #println("NAIVE")
                X = naive_assign_bu(instancia, Y_bool)
            elseif alloc_method == "opp"
                #println("OPP COST")
                X = oppCostAssignment(Y_bool, instancia)
            end
            after_alloc = Dates.now()
            delta_alloc = after_alloc - after_init
            delta_alloc_milli = round(delta_alloc, Millisecond)
            #println("X done with time: ", delta_alloc_milli)
            delta = after_alloc - before_init
            micros = round(delta, Millisecond)
            time_cons = micros.value
            if init_method == "relax"
                time_cons += (time_init * 1000)
            end
            indices = findall(x -> x == 1, X)
            for indice in indices
                Weight += D[indice]
            end    
            if !isdir("pruebas")
                mkdir("pruebas")
            end

            str_path = "out\\solutions\\800_150_50\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * init_method * "_" * alloc_method
            plot_str_path = str_path * ".png"
            solution_str_path = str_path * ".jld2"
            solution = Types.Solution(instancia, X, Y_bool, Weight, time_cons)

            Types.plot_solution(solution, plot_str_path)
            Types.write_solution(solution, solution_str_path)
            println(isFactible(solution, false))

            oldSol = solution
            oldTime = oldSol.Time
            factible_old, cons_v_og = isFactible(solution, false)
            really_cons_old = cons_v_og
            println("Constraints antes de LS: ", really_cons_old)
            weight_before = 1e12
            valor_previo = 0
            if factible_old
                #println("Solucion Factible")
                weight_before = solution.Weight
                valor_previo = weight_before
            else
                #println("Solucion no factible")
            end

           

            before_simple = Dates.now()
            solution, cons_v = simple_move_bu(solution, cons_v_og)
            println("Constraints despues de simple: ", cons_v)
            after_simple = Dates.now()
            mejora_cons_sim = cons_v_og - cons_v
            println("Mejora constraints simple: ", mejora_cons_sim)
            mejora_peso_sim = 0
            cons_v_og = cons_v
            if cons_v_og == 0
                valor_previo = solution.Weight
                if factible_old
                    mejora_peso_sim = weight_before - valor_previo
                end
            end
            delta_simple = after_simple - before_simple
            delta_simple_milli = round(delta_simple, Millisecond)
            #println("Simple terminado con tiempo: ", delta_simple_milli)
            before_inter = Dates.now()
            solution, cons_v = interchange_bus(solution, cons_v)
            println("Constraints despues de intercambio: ", cons_v)
            after_inter = Dates.now()
            mejora_cons_int = cons_v_og - cons_v
            println("Mejora constraints intercambio: ", mejora_cons_int)
            mejora_peso_int = 0
            cons_v_og = cons_v
            if cons_v_og == 0
                mejora_peso_int = valor_previo - solution.Weight
                valor_previo = solution.Weight
            end
            delta_inter = after_inter - before_inter
            delta_inter_milli = round(delta_inter, Millisecond)
            #println("Intercambio terminado con tiempo: ", delta_inter_milli)
            before_deacN = Dates.now()
            solution, cons_v = deactivateBranch(solution, cons_v)
            println("Constraints despues de desactivar N: ", cons_v)
            after_deacN = Dates.now()
            mejora_cons_deacN = cons_v_og - cons_v
            println("Mejora constraints deacN: ", mejora_cons_deacN)
            cons_v_og = cons_v
            mejora_peso_deacN = 0
            if cons_v_og == 0
                mejora_peso_deacN = valor_previo - solution.Weight
                valor_previo = solution.Weight
            end
            delta_deacN = after_deacN - before_deacN
            delta_deacN_milli = round(delta_deacN, Millisecond)
            before_deacS = Dates.now()
            solution, cons_v = deactivateBranch2(solution, cons_v)
            println("Constraints despues de desactivar S: ", cons_v)
            after_deacS = Dates.now()
            mejora_cons_deacS = cons_v_og - cons_v
            println("Mejora constraints deacS: ", mejora_cons_deacS)
            cons_v_og = cons_v
            mejora_peso_deacS = 0
            if cons_v_og == 0
                mejora_peso_deacS = valor_previo - solution.Weight
                valor_previo = solution.Weight
            end
            delta_deacS = after_deacS - before_deacS
            delta_deacS_milli = round(delta_deacS, Millisecond)
            val_delta = delta_deacS_milli.:value
            #println(val_delta)
            #println(typeof(val_delta))
            #println("Desactivar Smart terminado con tiempo: ", delta_deacS)
            delta_ls = after_deacS - before_simple
            micros_ls = round(delta_ls, Millisecond)
            time_ls = micros_ls.value
            total_time = time_ls + time_cons
            #println(total_time)
            final_weight = 1e12

            gap = 100
            if cons_v == 0
                final_weight = solution.Weight
                gap = ( 1 - (weight_exac / final_weight)) * 100
            end
            if init_method ≠ "relax"
                delta_init_milli = delta_init_milli.value
            end
            sol_ls = solution
            str_path = "out\\solutions\\800_150_50\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * init_method * "_" * alloc_method * "_ls"
            plot_str_path = str_path * ".png"
            solution_str_path = str_path * ".jld2"
            solution = Types.Solution(instancia, sol_ls.X, sol_ls.Y, final_weight, total_time)

            Types.plot_solution(solution, plot_str_path)
            Types.write_solution(solution, solution_str_path)
            if init_method ≠ "relax"
                if typeof(delta_init_milli) != Int64
                    delta_init_milli = delta_init_milli.value
                end
            end
            row_cons = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P"=> P, "Init" => init_method, "Alloc" => alloc_method, "TimeInit" => delta_init_milli, 
            "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "Value" =>weight_before, "TotalTimeCons"=> time_cons)
            push!(arr_constructive, row_cons)
            row_cons_ls = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P"=> P, "Init" => init_method, "Alloc" => alloc_method, "TimeInit" => delta_init_milli, 
            "TimeAlloc" => delta_alloc_milli.value, "Constraints" => really_cons_old, "ValuePreLS" =>weight_before, 
            "TimeSimple"=> delta_simple_milli.value, "ConsSimpleImprov" => mejora_cons_sim, "ValueSimpleImprov" => mejora_peso_sim, 
            "TimeInter"=> delta_inter_milli.value, "ConsInterImprov" => mejora_cons_int, "ValueInterImprov" => mejora_peso_int, 
            "TimeDeacN"=> delta_deacN_milli.value, "ConsDeacNImprov" => mejora_cons_deacN, "ValueDeacNImprov" => mejora_peso_deacN, 
            "TimeDeacS"=> val_delta, "ConsDeacSImprov" => mejora_cons_deacS, "ValueDeacSImprov" => mejora_peso_deacS, 
            "TotalTimeCons"=> time_cons, "TotalTimeLs" => time_ls, "TotalTime" => total_time, "ExactTime" => time_exac, 
            "ValuePostLs"=> final_weight, "Optim" => weight_exac, "Gap%"=> gap)
            push!(arr_ls, row_cons_ls)
            # println(row_cons_ls)
            end
        end
    df1 = vcat(DataFrame.(arr_constructive)...)
    df2 = vcat(DataFrame.(arr_ls)...)
    sort!(df1, :ID)
    sort!(df2, :ID)
    CSV.write("df_cons_800.csv", df1)
    CSV.write("df_ls_800.csv", df2)
end

main_test()
