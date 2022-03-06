using CPLEX, JuMP


function read_file(path::String)
    buffer = open(path)
    content = readlines(buffer)
    close(buffer)
    dist_ph = [] # distance placeholder
    distance_matrix = Vector{Vector{Float32}} # matrix placeholder
    startl = 0 # line in which distance matrix starts
    endl = 0 # line in which it ends
    for (index, line) in enumerate(content)
        if line == "D"
            startl = index + 1
        end
        if line == "S"
            endl = index - 1
            break # computational power isn't to be wasted
        end
    end
    # DISTANCE MATRIX STARTS HERE
    distances_str = content[startl:endl]
    for distance_str in distances_str
        distance_flt = [parse(Float32, distance) for distance in split(distance_str)]
        push!(dist_ph, distance_flt)
    end

    distance_matrix = hcat(dist_ph...) # https://stackoverflow.com/a/52257481
    num_BU, num_Suc = size(distance_matrix)

    # as for now the k param is hardcoded to 5
    k = 5
    Sk = [Set{Int64}(1) for _ in 1:k] # terrible workaround
    Sk_line = endl+2
    j = 1
    for i in Sk_line:length(content)
        if content[i] == "U"
            endl = i
            break
        end
        Sk[j] = Set([parse(Int, elem) for elem in split(content[i])])
        j = j +1
    end
    Uk = Array{Int64}(undef, k) # gotta move this to a function so DRY
    Uk_line = endl + 1
    j = 1
    for i in Uk_line:length(content)
        if content[i] == "L"
            endl = i
            break
        end
        Uk[j] = parse(Int, content[i])
        j = j +1
    end
    Lk = Array{Int64}(undef, k)
    Lk_line = endl + 1
    j = 1
    for i in Lk_line:length(content)
        if content[i] == "P"
            endl = i
            break
        end
        Lk[j] = parse(Int, content[i])
        j = j +1
    end
    endl = endl + 1
    P = parse(Int, content[endl])
    endl = endl + 2
    R = [parse(Int, elem) for elem in split(content[endl])]
    endl = endl + 2
    m = 3 # activities
    V = [Array{Int64}(undef,5) for _ in 1:m]
    j = 1
    for i in endl:length(content)
        if content[i] == "μ"
            endl = i
            break
        end
        V[j] = [parse(Int, value) for value in split(content[i])]
        j = j +1
    end
    endl = endl + 1
    μ = [Array{Int64}(undef,5) for _ in 1:m]
    j = 1
    for i in endl:length(content)
        if content[i] == "T"
            endl = i
            break
        end
        μ[j] = [parse(Int, value) for value in split(content[i])]
        j = j +1
    end
    endl = endl + 1
    T = Array{Float64}(undef, m)
    j = 1
    for i in endl:length(content)
        if content[i] == "α"
            endl = i
            break
        end
        T[j] = parse(Float32, content[i]) # consider truncing
        j = j +1
    end
    endl = endl + 1
    α = [parse(Int, thresh) for thresh in split(content[endl])]
    return distance_matrix, Sk, Uk, Lk, P, R, V, μ, T, α
    # falta agregar asserts para comprobar dimensiones
end

function model(data)

    D = data[1] # matrix float 32
    Sk = data[2] # vector any
    Uk = data[3] #
    Lk = data[4]
    P = data[5]
    R = data[6]
    V = data[7]
    μ = data[8]
    T = data[9]
    α = data[10]


end


function main()
    path = "../instancias/instancia.txt"
    #ARGS[1]
    data = read_file(path)
    model(data)
end

main()
