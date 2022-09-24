function maximums1(matrix, n)
    type = eltype(matrix)
    vals = zeros(type, n)
    indices = Array{CartesianIndex}(undef, n)
    for i ∈ axes(matrix, 1), j ∈ axes(matrix, 2)
        smallest, index = findmin(vals)
        if matrix[i, j] > smallest
            vals[index] = matrix[i, j]
            indices[index] = CartesianIndex(i, j)
        end
    end
    combined = hcat(vals, indices)
    sorted = sortslices(combined, dims=1, by=x -> (x[1], x[2]), rev=true)
    vals = sorted[:,1]
    indices = sorted[:,2]
    return vals, indices
end

function maximums2(matrix, n)
    type = eltype(matrix)
    vals = zeros(type, n)
    arr = Array{Tuple{type, CartesianIndex}}(undef, n)
    for i ∈ axes(matrix, 1), j ∈ axes(matrix, 2)
        smallest, index = findmin(vals)
        if matrix[i, j] > smallest
            arr[index] = matrix[i,j], CartesianIndex(i, j)
            vals[index] = matrix[i, j]
        end
    end
    arr = sort(arr, by = x->x[1], rev=true)
    vals = [x[1] for x in arr]
    indices = [x[2] for x in arr]
    return vals, indices
end

for i in 1:10
    D = rand(10:10000, 10,10)
    @time maximums1(D, 10) 
    @time maximums2(D, 10)
    println()
end