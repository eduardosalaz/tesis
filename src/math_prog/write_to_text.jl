using Main
using Types
function write_instance_txt(instance::Instance, filepath::String)
    open(filepath, "w") do file
        # Write scalar parameters with labels
        println(file, "PARAMS")
        println(file, "B ", instance.B, " S ", instance.S, " K ", instance.K, " M ", instance.M, " P ", instance.P)

        # Write BU coordinates
        println(file, "BU_COORDS")
        for i in 1:instance.B
            println(file, instance.BU_coords[i, 1], " ", instance.BU_coords[i, 2])
        end

        # Write S coordinates
        println(file, "S_COORDS")
        for i in 1:instance.S
            println(file, instance.S_coords[i, 1], " ", instance.S_coords[i, 2])
        end

        # Write distance matrix
        println(file, "DISTANCES")
        for i in 1:instance.S
            line = join([string(instance.D[i, j]) for j in 1:instance.B], " ")
            println(file, line)
        end

        # Write Sk lists
        println(file, "SK")
        for k in 1:instance.K
            println(file, "TYPE", k, " ", join(instance.Sk[k], " "))
        end

        # Write bounds
        println(file, "BOUNDS")
        println(file, "LK ", join(instance.Lk, " "))
        println(file, "UK ", join(instance.Uk, " "))

        # Write activity metrics
        println(file, "ACTIVITIES")
        for m in 1:instance.M
            println(file, "V", m, " ", join(instance.V[m], " "))
        end

        # Write capacity metrics
        println(file, "CAPACITIES")
        for m in 1:instance.M
            println(file, "MU", m, " ", join(instance.μ[m], " "))
        end

        # Write remaining values
        println(file, "TOLERANCE")
        println(file, join(instance.T, " "))
        println(file, "RISK")
        println(file, join(instance.R, " "))
        println(file, "RISK_CAP")
        println(file, join(instance.β, " "))
    end
    println("Instance written to ", filepath)
end

function main(ARGS)
    instance = read_instance(ARGS[1])
    write_instance_txt(instance, ARGS[2])
end
@main