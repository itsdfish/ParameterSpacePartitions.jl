# ParameterSpacePartitions.jl

Parameter space partitioning is a model assessment method for mapping regions of the parameter space to qualitative data patterns. Parameter space partitioning uses MCMC sampling to explore the parameter space. Each chain searches a region of the parameter space for a specific pattern. As the space is sampled, a new chain is created for each newly discovered pattern. On each iteration, a proposal is sampled uniformly from the surface of a hyperspere with the same number of dimensions as the parameter space. 

# Installation

You can install ParameterSpacePartitions.jl from the REPL by switching to the package mode using `]` and then typing:

```julia 
add ParameterSpacePartitions
```

# Quick Example

In this quick example, we use PSP to find the regions of a 2D polytope. The polytope formed by partitioning the space into 5 regions based on 5 centroids. A point in the space belongs to the region with the closest centroid. This involves defining a `Polytope` type, a `model` method, which, in this case, simply returns the proposed point, and a `pfun` function which returns the pattern associated with the point. The pattern is simply the id of the region associated with the closest centroid. As expected, the PSP algorithm finds all of the regions.

```@setup quick_example
using ParameterSpacePartitions
using Random
using LinearAlgebra
using DataFrames
using Distributions
using StatsPlots

struct Polytope
    location::Vector{Float64}
end

function p_fun(location, points::Vector{Polytope}, args...; kwargs...)
    distances = map(p -> norm(location .- p.location), points)
    _,idx = findmin(distances)
    return idx
end

model(parms, args...; kwargs...) = parms

# dimensions of the hypbercue
n_dims = 2
# number of partitions
n_part = 5
# number of starting points (only supports 1 currently)
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
    init_parms
)

df = find_partitions(
    model, 
    p_fun, 
    options,
    polytopes
)

n_patterns = length(unique(df.pattern))
```

```@example quick_example 
using ParameterSpacePartitions
using Random
using LinearAlgebra
using DataFrames
using Distributions
using StatsPlots

struct Polytope
    location::Vector{Float64}
end

function p_fun(location, points::Vector{Polytope}, args...; kwargs...)
    distances = map(p -> norm(location .- p.location), points)
    _,idx = findmin(distances)
    return idx
end

model(parms, args...; kwargs...) = parms

# dimensions of the hypbercue
n_dims = 2
# number of partitions
n_part = 5
# number of starting points (only supports 1 currently)
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
    init_parms
)

df = find_partitions(
    model, 
    p_fun, 
    options,
    polytopes
)

n_patterns = length(unique(df.pattern))
```
It is also possible to see the five regions by plotting the accepted samples as a scatter plot. 
```@example quick_example 
df_accepted = filter(x -> x.acceptance, df)
@df df_accepted scatter(:p1, :p2, group = :pattern, leg=false, grid=false)
```

# References

Pitt, M. A., Kim, W., Navarro, D. J., & Myung, J. I. (2006). Global model analysis by parameter space partitioning. Psychological Review, 113(1), 57.
