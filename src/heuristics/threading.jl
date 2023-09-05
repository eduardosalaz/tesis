using .Threads
using Dates

function test_load()
    vals = rand(100000000)
    sort!(vals)
    vals = vals .+ 1.24214212423431 .- 0.12232431253
    return (sum(vals) * 2344.2435343433)
end

function test_load2(valor)
    valor_int = round(Int, valor)
    val1 = round(Int, valor_int / 100000000000)
    val2 = factorial(big(val1))
    val3 = val2 + rand((1,2,3,4,5))
    return val3
end

function grasp1(iters)
    start = now()
    bestVal = 1
    for i in 1:iters
        println("iter: ", i)
        resultado = test_load()
        println("resultado1: ", resultado, " threadid: $(threadid())")
        resultado2 = test_load2(resultado)
        println("resultado2: ", resultado2, " threadid: $(threadid())")
        if resultado2 > bestVal
            println("Mejor resultado que el incumbente")
            bestVal = resultado2
        end
    end
    finish = now()
    println(bestVal)
    println("Runtime de: ", finish - start)
end

function grasp2(iters)
    start = now()
    bestVal = 1
    lockVar = ReentrantLock()
    Threads.@threads for i in 1:iters
        println("iter: ", i)
        resultado = test_load()
        println("resultado1: ", resultado, " threadid: $(threadid())")
        val1 = test_load2(resultado)
        println("resultado2: ", val1, " threadid: $(threadid())")
        lock(lockVar)
        try
            if val1 > bestVal
                println("Actualizando bestVal con: $val1, previous: $bestVal")
                bestVal = val1
            end
        finally
            unlock(lockVar)
        end
    end
    finish = now()
    println(bestVal)
    println("Runtime de: ", finish - start)
end

function main()
    println(Threads.nthreads())
    grasp1(10)
    grasp2(10)
end

main()