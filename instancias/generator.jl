using Plots, Distances

function generate()
    BU_coords = rand(10:500, (20,2))
    S_coords = rand(10:500,(8,2))
    plot(scatter(BU_coords[:,1], BU_coords[:,2], lab="BU", marker = ([:hex :d], 5, 0.8, Plots.stroke(3, :gray))))
    ploteado = scatter!(S_coords[:,1], S_coords[:,2], lab="S", c=:black, markersize=:6)
    png("ploteado.png")
    return BU_coords, S_coords
end

generate()
