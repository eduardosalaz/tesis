using Plots, Distances, DelimitedFiles

function generate()
    BU_coords = rand(10:500, (20,2))
    S_coords = rand(10:500,(8,2))
    plot(scatter(BU_coords[:,1], BU_coords[:,2], lab="BU", marker = ([:hex :d], 5, 0.8, Plots.stroke(3, :gray))))
    ploteado = scatter!(S_coords[:,1], S_coords[:,2], lab="S", c=:black, markersize=:6)
    png("ploteado.png")
    return BU_coords, S_coords
end

function distanceM(B,S)
    metrica = Euclidean()
    mat = zeros(8,20)
    for i in 1:8
        for j in 1:20
	   distancia = metrica(B[j,:], S[i,:])
	   mat[i,j] = distancia
	end
    end
    open("distance_mat.txt", "w") do io
        writedlm(io, mat)
    end
    return mat
end

B, S = generate()
matriz = distanceM(B,S)
