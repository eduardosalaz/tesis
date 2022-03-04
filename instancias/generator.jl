using Plots, Distances, DelimitedFiles

function generate()
    BU_coords = rand(10:500, (20, 2))
    S_coords = rand(10:500, (8, 2))
    plot(scatter(BU_coords[:, 1], BU_coords[:, 2], lab = "BU", marker = ([:hex :d], 5, 0.8, Plots.stroke(3, :gray))))
    ploteado = scatter!(S_coords[:, 1], S_coords[:, 2], lab = "S", c = :black, markersize = :6)
    png("plot_instance.png")
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


function write_file(mat, B, S)
    open("instancia1.txt", "w") do io
        write(io, "BU\n")
        writedlm(io, B, ' ')
        write(io, "S\n")
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
        write(io, "6\n") # solo 6 centros

    end
end




function main()
    B, S = generate()
    M = distanceM(B, S)
    write_file(M,B,S)
    @info "Done"
end
