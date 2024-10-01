include("constructive.jl")
include("ls.jl")
include("grasp_parallel.jl")
using DataFrames
using CSV
function main_test()
    entradas = readdir(ARGS[1])
    println(entradas)
    arr_grasp = []
    best_sols_gurobi_625 = Dict(
        "1" => 464246,
        "2" => 475106,
        "3" => 451879,
        "4" => 460877,
        "5" => 469136,
        "6" => 491677,
        "7" => 454550,
        "8" => 462146,
        "9" => 488765,
        "10" => 472705,
        "11" => 468053,
        "12" => 462072,
        "13" => 468505,
        "14" => 459301,
        "15" => 484674,
        "16" => 462365,
        "17" => 462745,
        "18" => 455692,
        "19" => 460741,
        "20" => 465070
    )
    
    best_bounds_gurobi_625 = Dict(
        "1" => 462006.983,
        "2" => 474508.581,
        "3" => 451829.403,
        "4" => 460361.988,
        "5" => 468430.222,
        "6" => 490028.106,
        "7" => 453977.334,
        "8" => 460849.226,
        "9" => 487419.813,
        "10" => 466940.063,
        "11" => 466771.252,
        "12" => 460985.561,
        "13" => 468312.438,
        "14" => 459248.658,
        "15" => 479009.77,
        "16" => 462309.765,
        "17" => 462066.303,
        "18" => 454192.089,
        "19" => 458314.077,
        "20" => 464814.975
    )
    
    time_gurobi_625 = Dict(
        "1" => 1800,
        "2" => 1800,
        "3" => 1483,
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
        "14" => 976,
        "15" => 1800,
        "16" => 1384,
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
        time_exac = time_gurobi_625[number]

        B = instancia.B
        S = instancia.S
        P = instancia.P
        X = Matrix{Int64}(undef, S, B)
        Y = Vector{Int64}(undef, S)
        D = instancia.D
        Weight = 0
        #test_αₗ = [0.1, 0.3, 0.5, 0.7, 0.9]
        test_αₗ = [0.1]
        #test_αₐ = [0.1, 0.3, 0.5, 0.7, 0.9]
        test_αₐ = [0.3]
        #test_iters = [10, 30, 50, 70, 90]
        test_iters = [70]
        before = now()
        for αₗ in test_αₗ
            for αₐ in test_αₐ
                for iters in test_iters
                    println("$αₗ, $αₐ, $iters")
                    sol, time_grasp = grasp(αₗ, αₐ,  iters, instance)
                    if sol !== nothing
                        bestWeight = sol.Weight

                        diff_grasp = abs(bestWeight - best_sols_gurobi_625[number]) / abs(bestWeight)
                        gap_grasp = abs(bestWeight - best_bounds_gurobi_625[number]) / abs(bestWeight)
                        time_exac = time_gurobi_625[number]
                        beat_gurobi = false
                        if bestWeight < best_sols_gurobi_625[number]
                            beat_gurobi = true
                        end

                        str_path = "out\\solutions\\625_78_32\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * "$αₗ" * "_" * "$αₐ" * "_" * "$iters" * "_grasp_final"
                        plot_str_path = str_path * ".png"
                        solution_str_path = str_path * ".jld2"
                        #Types.plot_solution(sol, plot_str_path)
                        Types.write_solution(sol, solution_str_path)
                        
                        # gap_repaired = (1 - (weight_exac / bestWeight)) * 100
                        row_grasp = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => "pdisp", "Alloc" => "queue", "TotalTimeGRASP" => time_grasp, "TimeOptim" => time_exac,
                        "ValueOptim" => weight_exac, "BestBound" => best_bounds_gurobi_625[number], "ValueGRASP" => bestWeight, "Gap%"=> gap_grasp, "Diff%"=> diff_grasp, "alpha_l"=>αₗ, "alpha_a"=>αₐ,
                        "BeatGurobi"=>beat_gurobi, "iters"=>iters)
                        push!(arr_grasp, row_grasp)
                    else
                        row_grasp = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => "pdisp", "Alloc" => "queue", "TotalTimeGRASP" => 0, "TimeOptim" => time_exac,
                        "ValueOptim" => weight_exac, "BestBound"=> best_bounds_gurobi_625[number], "ValueGRASP" => 0, "Gap%"=> 100, "Diff%"=>100, "alpha_l"=>αₗ, "alpha_a"=>αₐ, "iters"=>iters, "BeatGurobi"=>false)
                        push!(arr_grasp, row_grasp)
                    end
                    
                end
            end
        end
        df1 = vcat(DataFrame.(arr_grasp)...)
        df1 = df1[:, sortperm(names(df1))]
        CSV.write("df_grasp_625_78_32_70iters.csv", df1)
    end
end

main_test()
println("\a")
