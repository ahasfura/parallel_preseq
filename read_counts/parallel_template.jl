using DistributedArrays

@everywhere function calculation(num)
    return num*num
end

println("Worker IDs: ", workers())
println()

indices = [10, 20, 30, 40, 50, 60]
j = distribute(indices)

for i in workers()
    r = @spawn localpart(j)
    @printf("%d\t%s\n", r.where, fetch(r))
end

