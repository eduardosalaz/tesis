using JLD2

function plot_instance(Instance::Instance, path::String)
    Plots.scatter(
        Instance.BU_coords[:, 1],
        Instance.BU_coords[:, 2],
        label = "BUs",
        markershape = :circle,
        markercolor = :blue,
    )
    Plots.scatter!(
        Instance.S_coords[:, 1],
        Instance.S_coords[:, 2],
        label = "Branches",
        markershape = :square,
        markercolor = :white,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    sizeB = Instance.B
    sizeS = Instance.S
    png(path)
    @debug "Wrote plot and coords"
end

function plot_solution(Solution::Solution, path::String)
    Instance = Solution.Instance
    BU_coords = Instance.BU_coords
    S_coords = Instance.S_coords
    Y = Solution.Y
    X = Solution.X
    S = Instance.S
    B = Instance.B
    is_used = [(Y[j] ? :red : :white) for j in 1:S]
    Sk = Instance.Sk
    S₁ = Sk[1]
    S₁_coords = S_coords[S₁, :]
    S₂ = Sk[2]
    S₂_coords = S_coords[S₂, :]
    S₃ = Sk[3]
    S₃_coords = S_coords[S₃, :]
    S₄ = Sk[4]
    S₄_coords = S_coords[S₄, :]
    S₅ = Sk[5]
    S₅_coords = S_coords[S₅, :]
    Plots.scatter(
        BU_coords[:, 1],
        BU_coords[:, 2],
        label = "BUs",
        markershape = :circle,
        markercolor = :blue,
    )
    Plots.scatter!(
        S₁_coords[:, 1],
        S₁_coords[:, 2],
        label = "Branches K₁",
        markershape = :hexagon,
        markercolor = is_used,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    Plots.scatter!(
        S₂_coords[:, 1],
        S₂_coords[:, 2],
        label = "Branches K₂",
        markershape = :diamond,
        markercolor = is_used,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    Plots.scatter!(
        S₃_coords[:, 1],
        S₃_coords[:, 2],
        label = "Branches K₃",
        markershape = :heptagon,
        markercolor = is_used,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    Plots.scatter!(
        S₄_coords[:, 1],
        S₄_coords[:, 2],
        label = "Branches K₄",
        markershape = :pentagon,
        markercolor = is_used,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    Plots.scatter!(
        S₅_coords[:, 1],
        S₅_coords[:, 2],
        label = "Branches K₅",
        markershape = :star4,
        markercolor = is_used,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )

    for i in 1:B
        for j in 1:S
            if X[i, j] == 1
                Plots.plot!(
                    [BU_coords[i,1], S_coords[j,1]],
                    [BU_coords[i,2], S_coords[j,2]],
                    color = :black,
                    label = nothing,
                )
            end
        end
    end
    png(path)
    @debug "Wrote plot"

end

function read_instance(path::String)
    instancia = jldopen(path)
    instance = instancia["instance"]
    return instance
end

function write_instance(instance::Instance, path::String)
    jldsave(path; instance)
end

function read_solution(path::String)
    sol = jldopen(path)
    solution = sol["solution"]
    return solution
end

function write_solution(solution::Solution, path::String)
    jldsave(path; solution)
end
