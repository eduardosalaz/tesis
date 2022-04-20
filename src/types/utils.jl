using JLD2

function plot_instance(Instance::Instance)
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
    png("out/plots/plot_instance_bu$size" * "_suc" * "$sizeS" * ".png")
    #@info "Wrote plot and coords"
end

function plot_solution(Solution::Solution)
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
