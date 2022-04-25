module Types
using JLD2, Plots
include("utils.jl")
struct Instance
    B::Int64
    S::Int64
    K::Int64
    M::Int64
    P::Int64
    BU_coords::Matrix{Int64}
    S_coords::Matrix{Int64}
    D::Matrix{Int64}
    Sk::Vector{Vector{Int64}}
    Lk::Vector{Int64}
    Uk::Vector{Int64}
    V::Vector{Vector{Int64}}
    μ::Vector{Vector{Int64}}
    T::Vector{Float64}
    R::Vector{Int64}
    β::Vector{Int64}
end

struct Solution
    Instance::Instance
    X::Matrix{Int64}
    Y::Vector{Int64}
    Weight::Int64
end

function plot_instance(Instance, path::String)
    Plots.scatter(
        Instance.BU_coords[:, 1],
        Instance.BU_coords[:, 2],
        label = nothing,
        markershape = :circle,
        markercolor = :blue,
    )
    Plots.scatter!(
        Instance.S_coords[:, 1],
        Instance.S_coords[:, 2],
        label = nothing,
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

function plot_solution(Solution, path::String)
    Instance = Solution.Instance
    BU_coords = Instance.BU_coords
    S_coords = Instance.S_coords
    Y = Solution.Y
    # Y = iszero.(Y)
    Y2 = BitArray(Y)
    X = Solution.X
    S = Instance.S
    B = Instance.B
    Sk = Instance.Sk
    S₁ = Sk[1]
    S₁_coords = S_coords[S₁, :]
    Y₁_used = [(Y2[j] ? :red : :white) for j in S₁]
    S₂ = Sk[2]
    S₂_coords = S_coords[S₂, :]
    Y₂_used = [(Y2[j] ? :red : :white) for j in S₂]
    S₃ = Sk[3]
    S₃_coords = S_coords[S₃, :]
    Y₃_used = [(Y2[j] ? :red : :white) for j in S₃]
    S₄ = Sk[4]
    S₄_coords = S_coords[S₄, :]
    Y₄_used = [(Y2[j] ? :red : :white) for j in S₄]
    S₅ = Sk[5]
    S₅_coords = S_coords[S₅, :]
    Y₅_used = [(Y2[j] ? :red : :white) for j in S₅]
    Plots.scatter(
        BU_coords[:, 1],
        BU_coords[:, 2],
        markershape = :circle,
        markercolor = :blue,
        label = nothing
    )
    Plots.scatter!(
        S₁_coords[:, 1],
        S₁_coords[:, 2],
        label = nothing,
        markershape = :hexagon,
        markercolor = Y₁_used,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    Plots.scatter!(
        S₂_coords[:, 1],
        S₂_coords[:, 2],
        label = nothing,
        markershape = :diamond,
        markercolor = Y₂_used,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    Plots.scatter!(
        S₃_coords[:, 1],
        S₃_coords[:, 2],
        label = nothing,
        markershape = :star5,
        markercolor = Y₃_used,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    Plots.scatter!(
        S₄_coords[:, 1],
        S₄_coords[:, 2],
        label = nothing,
        markershape = :pentagon,
        markercolor = Y₄_used ,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )
    Plots.scatter!(
        S₅_coords[:, 1],
        S₅_coords[:, 2],
        label = nothing,
        markershape = :star4,
        markercolor = Y₅_used ,
        markersize = 6,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )

    for i in 1:B
        for j in 1:S
            if X[j,i] == 1
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

function write_instance(instance, path::String)
    jldsave(path; instance)
end

function read_solution(path::String)
    sol = jldopen(path)
    solution = sol["solution"]
    return solution
end

function write_solution(solution, path::String)
    jldsave(path; solution)
end

end
