using DelimitedFiles, Distances

function remove!(V, item)
    deleteat!(V, findall(x->x==item, V))
end
function main()
    m = readdlm("src\\heuristics\\data.txt")
    coords = m[:, 2:3]
    metric = Euclidean()
    s_distances = pairwise(metric, coords, dims=1)
    S = 100
    P = 10
    idx_1, idx_2 = Tuple(argmax(s_distances))
    S_sol = [idx_1, idx_2] # S is our solution of our p-disp problem
    T = collect(1:S) # nodes not in our solution
    for s in S_sol
        remove!(T, s)
    end
    while length(S_sol) < P
        
    end
    println(S_sol)
end
main()




