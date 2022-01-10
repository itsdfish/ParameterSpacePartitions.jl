cd(@__DIR__)
using Pkg
Pkg.activate("")
using Revise, Distributions, ParameterSpacePartitions
using ParameterSpacePartitions.TestModels
using LinearAlgebra, Random, DataFrames, StatsPlots

# dimensions of the hypbercue
n_dims = 10
# number of partitions
n_obj = 100
# distribution of radii
r_dist = () -> rand(Uniform(.2, .25))
# number of starting points
n_start = 1

# partition boundaries
hyperspheres = set_locations(r_dist, n_dims, n_obj)
bounds = fill((0, 1), n_dims)

sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

init_parms = map(_ -> sample(bounds), 1:n_start)

options = Options(;
    radius = .4,
    bounds,
    n_iters = 10_000,
    parallel = false,
    init_parms
)

results = find_partitions(
    model, 
    x -> p_fun(x, hyperspheres), 
    options
)

df = DataFrame(results)


# transform!(
#     df, 
#     :parms => identity => [:p1, :p2]
# )



# groups = groupby(df, :pattern)

# @df df scatter(:p1, :p2, group=:pattern, leg=false)