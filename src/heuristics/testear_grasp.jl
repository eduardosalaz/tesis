include("constructive.jl")
include("ls.jl")
include("grasp_parallel.jl")
using DataFrames
using CSV
function main_test()
    entradas = readdir(ARGS[1])
    arr_grasp = []
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
        X = Matrix{Int64}(undef, S, B)
        Y = Vector{Int64}(undef, S)
        D = instancia.D
        Weight = 0
        test_αₗ = [0.1, 0.3, 0.5, 0.7, 0.9]
        #test_αₗ = [0.1, 0.3]
        test_αₐ = [0.1, 0.3, 0.5, 0.7, 0.9]
        # test_αₐ = [0.1, 0.3]
        test_iters = [10, 30, 50, 70, 90]
        before = now()
        for αₗ in test_αₗ
            for αₐ in test_αₐ
                for iters in test_iters
                    println("$αₗ, $αₐ, $iters")
                    sol, time_grasp = grasp(αₗ, αₐ,  iters, instance)
                    if sol !== nothing
                        str_path = "out\\solutions\\625_78_32\\heurs\\sol" * "_" * string(id) * "_" * "$B" * "_" * "$S" * "_" * "$P" * "_" * "$αₗ" * "_" * "$αₐ" * "_" * "$iters" * "_grasp"
                        plot_str_path = str_path * ".png"
                        solution_str_path = str_path * ".jld2"
                        Types.plot_solution(sol, plot_str_path)
                        Types.write_solution(sol, solution_str_path)
                        bestWeight = sol.Weight
                        gap_repaired = (1 - (weight_exac / bestWeight)) * 100
                        row_grasp = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => "pdisp", "Alloc" => "queue", "TotalTimeGRASP" => time_grasp, "TimeOptim" => time_exac*1000,
                        "ValueOptim" => weight_exac, "ValueGRASP" => bestWeight, "Gap%"=> gap_repaired, "alpha_l"=>αₗ, "alpha_a"=>αₐ, "iters"=>iters)
                        push!(arr_grasp, row_grasp)
                    else
                        row_grasp = Dict("ID" => parse(Int, id), "B" => B, "S" => S, "P" => P, "Loc" => "pdisp", "Alloc" => "queue", "TotalTimeGRASP" => 0, "TimeOptim" => time_exac*1000,
                        "ValueOptim" => weight_exac, "ValueGRASP" => 0, "Gap%"=> 0, "alpha_l"=>αₗ, "alpha_a"=>αₐ, "iters"=>iters)
                        push!(arr_grasp, row_grasp)
                    end
                    
                end
            end
        end
        df1 = vcat(DataFrame.(arr_grasp)...)
        df1 = df1[:, sortperm(names(df1))]
        CSV.write("df_grasp_625_all_combs.csv", df1)
    end
end

main_test()
println("\a")
