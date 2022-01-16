cd(@__DIR__)
using Pkg
Pkg.activate("..")
using ParameterSpacePartitions
using ParameterSpacePartitions.TestModels
using Revise, Random, Distributions
using DataFrames, LinearAlgebra
#Random.seed!(2001)

# dimensions of the hypbercue
n_dims = 2
# number of partitions
n_part = 50
# number of starting points
n_start = 1


# partition boundaries
polytopes = [Polytope(rand(n_dims)) for i in 1:n_part]
bounds = fill((0, 1), n_dims)

sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

init_parms = map(_ -> sample(bounds), 1:n_start)

options = Options(;
    radius = .1,
    bounds,
    n_iters = 500,
    parallel = false,
    init_parms,
    λ = .2
)

results = find_partitions(
    model, 
    p_fun, 
    options,
    polytopes
)

df = DataFrame(results)
transform!(
    df, 
    :parms => identity => [:p1, :p2]
)

groups = groupby(df, :pattern)

groups = groupby(df, :chain_id)
mean_accept = combine(groups, :acceptance=>mean=>:mean)


groups = groupby(df, :pattern)

@df df scatter(:p1, :p2, group=:pattern, leg=false)