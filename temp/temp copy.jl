cd(@__DIR__)
using Pkg
Pkg.activate("")
using Revise, Distributions, ParameterSpacePartitions
using LinearAlgebra, Random, DataFrames, StatsPlots

function p_fun(data, points)
    distances = map(p -> norm(data .- p), points)
    _,idx = findmin(distances)
    return idx
end

# dimensions of the hypbercue
n_dims = 2
# number of partitions
n_part = 50
# number of starting points
n_start = 1
#

# partition boundaries
points = [rand(2) for i in 1:n_part]
bounds = fill((0, 1), n_dims)

model(parms) = parms 

sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

init_parms = map(_ -> sample(bounds), 1:n_start)

options = Options(;
    radius = .1,
    bounds,
    n_iters = 500,
    parallel = false,
    init_parms
)

results = find_partitions(
    model, 
    x -> p_fun(x, points), 
    options
)

df = DataFrame(results)
transform!(
    df, 
    :parms => identity => [:p1, :p2]
)



groups = groupby(df, :pattern)

@df df scatter(:p1, :p2, group=:pattern, leg=false)