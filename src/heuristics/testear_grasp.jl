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
        "10" => 999440,
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
        "1" => 944849.511,
        "2" => 967012.289,
        "3" => 982048.419,
        "4" => 971396.864,
        "5" => 957950.011,
        "6" => 961976.211,
        "7" => 963149.897,
        "8" => 983207.456,
        "9" => 947419.238,
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

    best_sols_gurobi_2500 = Dict(
        "1" => 883693,
        "2" => 893408,
        "3" => 882602,
        "4" => 913347,
        "5" => 895581,
        "6" => 896157,
        "7" => 893684,
        "8" => 895125,
        "9" => 877510,
        "10" => 892898
    )

    best_bounds_gurobi_2500 = Dict(
        "1" => 881117.835,
        "2" => 892105.537,
        "3" => 882275.835,
        "4" => 910961.688,
        "5" => 893657.944,
        "6" => 894518.511,
        "7" => 892487.885,
        "8" => 892662.594,
        "9" => 876485.604,
        "10" => 891652.122
    )

    time_gurobi_2500 = Dict(
        "1" => 7200,
        "2" => 7200,
        "3" => 7200,
        "4" => 7200,
        "5" => 7200,
        "6" => 7200,
        "7" => 7200,
        "8" => 7200,
        "9" => 7200,
        "10" => 7200
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
        weight_exac = best_sols_gurobi_1250[number]
        time_exac = time_gurobi_1250[number]
        bound = best_bounds_gurobi_1250[number]

        B = instancia.B
        S = instancia.S
        P = instancia.P
        X = Matrix{Int64}(undef, S, B)
        Y = Vector{Int64}(undef, S)
        D = instancia.D
        Weight = 0
        test_αₗ = [0.1, 0.3, 0.5, 0.7]
        #test_αₗ = [0.1]
        test_αₐ = [0.1, 0.3, 0.5, 0.7]
        #test_αₐ = [0.1]
        test_iters = [25, 50, 75, 100]
        #test_iters = [100]
        before = now()
        for αₗ in test_αₗ
            for αₐ in test_αₐ
                for iters in test_iters
                    println("$αₗ, $αₐ, $iters")
                    sol, time_grasp = grasp(αₗ, αₐ,  iters, instance)
                    if sol !== nothing
                        bestWeight = sol.Weight

                        diff_grasp = abs(bestWeight - weight_exac) / abs(bestWeight)
                        gap_grasp = abs(bestWeight - bound) / abs(bestWeight)
                        #time_exac = time_gurobi_1250[number]
                        beat_gurobi = false
                        if bestWeight < weight_exac
                            beat_gurobi = true
                        end

                        str_path = "out\\solutions\\1250_155_62\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * "$αₗ" * "_" * "$αₐ" * "_" * "$iters" * "_grasp_tuning"
                        plot_str_path = str_path * ".png"
                        solution_str_path = str_path * ".jld2"
                        #Types.plot_solution(sol, plot_str_path)
                        Types.write_solution(sol, solution_str_path)
                        
                        # gap_repaired = (1 - (weight_exac / bestWeight)) * 100
                        row_grasp = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => "pdisp", "Alloc" => "queue", "TotalTimeGRASP" => time_grasp, "TimeOptim" => time_exac,
                        "ValueOptim" => weight_exac, "BestBound" => bound, "ValueGRASP" => bestWeight, "Gap%"=> gap_grasp, "Diff%"=> diff_grasp, "alpha_l"=>αₗ, "alpha_a"=>αₐ,
                        "BeatGurobi"=>beat_gurobi, "iters"=>iters)
                        push!(arr_grasp, row_grasp)
                    else
                        row_grasp = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => "pdisp", "Alloc" => "queue", "TotalTimeGRASP" => 0, "TimeOptim" => time_exac,
                        "ValueOptim" => weight_exac, "BestBound"=> bound, "ValueGRASP" => 0, "Gap%"=> 100, "Diff%"=>100, "alpha_l"=>αₗ, "alpha_a"=>αₐ, "iters"=>iters, "BeatGurobi"=>false)
                        push!(arr_grasp, row_grasp)
                    end
                    
                end
            end
        end
        df1 = vcat(DataFrame.(arr_grasp)...)
        df1 = df1[:, sortperm(names(df1))]
        CSV.write("df_grasp_1250_allcombinations_insts89.csv", df1)
    end
end

main_test()
println("\a")
