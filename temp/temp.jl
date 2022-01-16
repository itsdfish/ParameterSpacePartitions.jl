# https://cookierobotics.com/007/
cd(@__DIR__)
using Pkg
Pkg.activate("..")
using Revise, Distributions, ParameterSpacePartitions
using ParameterSpacePartitions.TestModels
using Random, DataFrames, StatsBase
#Random.seed!(554)

# dimensions of the hypbercue
n_dims = 3
# partitions per dimension
n_part = 2
# number of starting points
n_start = 1

# partition boundaries
bounds = fill((0, 1), n_dims)
p_bounds = range(0, 1, n_part + 1)
hypercube = HyperCube(p_bounds)

sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

init_parms = map(_ -> sample(bounds), 1:n_start)

options = Options(;
    radius = .40,
    bounds,
    n_iters = 500,
    parallel = false,
    init_parms,
    t_rate = .3,
    λ = .2
)

results = find_partitions(
    model, 
    p_fun,
    options,
    hypercube
)

df = DataFrame(results)

transform!(
    df, 
    :parms => identity => [:p1, :p2, :p3]
)

transform!(df, :pattern => denserank => :pattern_id)

groups = groupby(df, :chain_id)
mean_accept = combine(groups, :acceptance=>mean=>:mean)

groups = groupby(df, :pattern)
summary = combine(groups, 
    :p1 => minimum, :p1 => maximum, 
    :p2 => minimum, :p2 => maximum,
    :p3 => minimum, :p3 => maximum
) |> sort

using GLMakie, ColorSchemes

# transform pattern into integer id
transform!(df, :pattern => denserank => :pattern_id)

fig, ax, scat = scatter(
    df.p1,
    df.p2,
    df.p3, 
    color = df.pattern_id, 
    axis = (;type=Axis3),
    markersize = 5000,
    colormap = ColorSchemes.seaborn_deep6.colors,
    grid = true
)