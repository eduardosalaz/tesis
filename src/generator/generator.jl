include("../types/types.jl")
using Distances, JLD2, Plots, .Types

function write_jld2(size::String)
    B, S = 0
    if size == "S"
        B = 15
        S = 5
        BU_coords, S_coords = generate_coords(B, S)
        dist_mat = distanceM(BU_coords, S_coords, B, S)
        parameters = generate_params(size)



    elseif size == "M"

    elseif size == "L"

    elseif size == "XL"

    else
        @error "Size not recognized"
        exit(1)
    end

end


function generate_coords(B, S)
    BU_coords = rand(10:500, (B, 2))
    S_coords = rand(10:500, (S, 2))

    Plots.scatter(
        BU_coords[:, 1],
        BU_coords[:, 2],
        label = "BUs",
        markershape = :circle,
        markercolor = :blue,
    )
    Plots.scatter!(
        S_coords[:, 1],
        S_coords[:, 2],
        label = "Branches",
        markershape = :square,
        markercolor = :white,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    png("../out/plot_instance_bu$B" * "_suc" * "$S" * ".png")
    return BU_coords, S_coords
end

function distanceM(BU_coords, S_coords, B, S)

    metrica = Euclidean()
    mat = zeros(S, B)
    for i in 1:S
        for j in 1:B
            distancia = metrica(BU_coords[j, :], S_coords[i, :])
            mat[i, j] = distancia
        end
    end

    return trunc.(Int, mat)
end

function generate_params(size::String)
    M = 3
    K = 5
    params = nothing
    if size == "S"


    elseif size == "M"

    elseif size == "L"

    elseif size == "XL"

    end


end


function write_file(mat, B, S, num_BU, num_Suc)
    open("instance_numBu$num_BU" * "_num_Suc" * "$num_Suc" * ".txt", "w") do io
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
    size = ARGS[1]
    write_jld2(size)
    @info "Done"
end

main()
