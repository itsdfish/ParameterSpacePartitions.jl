cd(@__DIR__)
using Pkg
Pkg.activate("..")
using Revise, Distributions, ParameterSpacePartitions
using ParameterSpacePartitions.TestModels
using LinearAlgebra, Random, DataFrames

#Random.seed!(54545)
# dimensions of the hypbercue
n_dims = 7
# number of partitions
n_obj = 50
# distribution of radii
r_dist = () -> rand(Uniform(.15, .25))
# number of starting points
n_start = 10

# partition boundaries
hyperspheres = set_locations(r_dist, n_dims, n_obj)
bounds = fill((0, 1), n_dims)

sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

init_parms = map(_ -> sample(bounds), 1:n_start)

options = Options(;
    radius = .25,
    bounds,
    n_iters = 5000,
    parallel = false,
    init_parms,
    λ = .2,
    t_rate = .40,
    adapt_radius! = no_adaption!,
    trace_on = false
)

results = find_partitions(
    model, 
    p_fun, 
    options,
    hyperspheres;
)

df = DataFrame(results)

parm_names = Symbol.("p",1:n_dims)
transform!(
    df, 
    :parms => identity => parm_names
)

groups = groupby(df, :chain_id)
mean_accept = combine(groups, :acceptance=>mean=>:mean)
