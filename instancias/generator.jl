using Plots, Distances, DelimitedFiles

function generate(num_BU, num_Suc)
    BU_coords = rand(10:500, (num_BU, 2))
    S_coords = rand(10:500, (num_Suc, 2))


    Plots.scatter(
    BU_coords[:,1],
    BU_coords[:,2],
    label = "BUs",
    markershape = :circle,
    markercolor = :blue,
)
Plots.scatter!(
    S_coords[:,1],
    S_coords[:,2],
    label = "Branches",
    markershape = :square,
    markercolor = :white,
    markersize = 6,
    markerstrokecolor = :red,
    markerstrokewidth = 2,
)

    # plot(scatter(BU_coords[:, 1], BU_coords[:, 2], lab = "BU", marker = ([:hex :d], 5, 0.8, Plots.stroke(3, :gray))))
    # ploteado = scatter!(S_coords[:, 1], S_coords[:, 2], lab = "S", c = :black, markersize = :6)
    png("plot_instance_bu$num_BU"*"_suc" *"$num_Suc" *".png")
    return BU_coords, S_coords
end

function distanceM(B, S)
    metrica = Euclidean()
    mat = zeros(8, 20)
    for i in 1:8
        for j in 1:20
            distancia = metrica(B[j, :], S[i, :])
            mat[i, j] = distancia
        end
    end

    return mat
end


function write_file(mat, B, S, num_BU, num_Suc)
    open("instance_numBu$num_BU"*"_num_Suc"*"$num_Suc"*".txt", "w") do io
        write(io, "BU\n")
        writedlm(io, B, ' ')
        write(io, "Suc\n")
        writedlm(io, S, ' ')
        write(io, "D\n")
        writedlm(io, mat, ' ')
        write(io, "S\n")
        write(io, "1 3\n") # s1
        write(io, "2 5\n") # s2
        write(io, "4\n") # s3
        write(io, "6 7\n") #s4
        write(io, "8\n") # s5

        write(io, "U\n")
        write(io, "2\n") #u1
        write(io, "2\n")
        write(io, "1\n")
        write(io, "2\n")
        write(io, "1\n")

        write(io, "L\n")
        write(io, "1\n") # l1  rango de 1 ≤ 2 = 1 o 2
        write(io, "1\n") # l2 rango de 1 ≤ 2 = 1 o 2
        write(io, "1\n") # l3 rango de 1 ≤ 1 = 1
        write(io, "1\n") # l4 rango de 1 ≤ 2 = 1 o 2
        write(io, "0\n") # l5 rango de 0 ≤ 1 = 0 o 1

        write(io, "P\n")
        write(io, "8\n") # solo 8 centros

    end
end




function main()
    num_BU = parse(Int, ARGS[1])
    num_Suc = parse(Int, ARGS[2])
    B, S = generate(num_BU, num_Suc)
    M = distanceM(B, S)
    write_file(M,B,S, num_BU, num_Suc)
    @info "Done"
end

main()
