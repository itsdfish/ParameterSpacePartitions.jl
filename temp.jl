cd(@__DIR__)
using Pkg
Pkg.activate("")
using Revise, Distributions, ParameterSpacePartitions

n_dims = 2
bounds = fill((0, 5), n_dims)

radius = .2

model(parms) = parms 

p_fun(data) = @. Int(floor(data))

sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

init_parms = map(_ -> sample(bounds), 1:10)

options = Options(;
    radius,
    bounds,
    n_iters = 1000,
    parallel = false,
    init_parms
)

results = find_partitions(model, p_fun, options)

df = DataFrame(results)
