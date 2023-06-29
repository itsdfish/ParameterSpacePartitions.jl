```@setup example1 
using Distributions
using ParameterSpacePartitions
using ParameterSpacePartitions.TestModels
using Random
using DataFrames
Random.seed!(544)

model(parms, args...; kwargs...) = parms 

function p_fun(location, hypercube::HyperCube, args...; kwargs...)
    p_bounds = hypercube.p_bounds
    nb = length(p_bounds)
    nd = length(location)
    vals = fill(-100, nd)
    for j in 1:nd
        for i in 1:(nb-1) 
            if (location[j] ≥ p_bounds[i]) && (location[j] ≤ p_bounds[i+1])
                vals[j] = i 
                continue
            end
        end
    end
    return vals
end

# dimensions of the hypbercue
n_dims = 3
# partitions per dimension
n_part = 2

# partition boundaries
bounds = fill((0, 1), n_dims)
p_bounds = range(0, 1, n_part + 1)
hypercube = HyperCube(p_bounds)


# number of starting points
n_start = 1
# sample function
sample(bounds) = map(b -> rand(Uniform(b...)), bounds)
# initial starting points
init_parms = map(_ -> sample(bounds), 1:n_start)
parm_names = [:p1, :p2, :p3]

options = Options(;
    radius = .10,
    bounds,
    n_iters = 1000,
    init_parms
)

results = find_partitions(
    model, 
    p_fun, 
    options,
    hypercube,
)

df = DataFrame(results)

groups = groupby(df, :pattern)
summary = combine(groups, 
    :p1 => minimum, :p1 => maximum, 
    :p2 => minimum, :p2 => maximum,
    :p3 => minimum, :p3 => maximum
) |> sort


groups = groupby(df, :pattern)

df_volume = combine(
    groups,
    x -> estimate_volume(
        model,
        p_fun, 
        x,  
        bounds,
        hypercube; 
        parm_names,
    )
)

df_volume.volume = df_volume.x1 / sum(df_volume.x1)
```

# Example

In this simple example, parameter space partitioning is applied to a cube with regions of equal volume.
## Load Dependencies

The first step is to the load dependences into your session.

```@example example1
using Distributions
using ParameterSpacePartitions
using ParameterSpacePartitions.TestModels
using Random
using DataFrames
Random.seed!(544)
```

## Model Function 

Next, we define a model that accepts a vector of parameters and returns data predictions. The function `model` is imported above from `TestModels`, but is presented below for illustration. In this simple example, the model returns the parameter inputs. In more substantive model, the returned value will be the model predictions.

```@example example1
model(parms, args...; kwargs...) = parms 
```
## Pattern Function

A second function categorizes the predicted data into a qualitative pattern. At minimum, the pattern function must recieve a data input. In this example, the pattern function `p_fun` is imported from `TestModels`. The function is presented below for illustration. `p_fun` requires the location (i.e. parameters) and `HyperCube` object, which contains  partition boundaries (which is the same for each dimension). `p_fun` categorizes `location` on each dimension and returns a vector representing a pattern. For example, the pattern `[1,2,1]` indicates the `location` is in the partition for which the x axis is 1, the y axis is 2, and the z axis is 1. 

```@example example1
function p_fun(location, hypercube::HyperCube, args...; kwargs...)
    p_bounds = hypercube.p_bounds
    nb = length(p_bounds)
    nd = length(location)
    vals = fill(-100, nd)
    for j in 1:nd
        for i in 1:(nb-1) 
            if (location[j] ≥ p_bounds[i]) && (location[j] ≤ p_bounds[i+1])
                vals[j] = i 
                continue
            end
        end
    end
    return vals
end
```

## Model Configuration

Below, we will set `n_dims` and `n_part` to create a cube with a total of 2^3 = 8 regions.

```@example example1
# dimensions of the hypbercue
n_dims = 3
# partitions per dimension
n_part = 2
```

# Boundaries

In the code below, `bounds` contains the upper and lower boundaries for each dimension. In addition, we must also define a `HyperCube` object containing the boundaries for each partition. The boundaries are stored in a variable called `p_bounds`.  In this example, the partitions are equally spaced along each dimension of the unit cube. Alternatively, we could use unequal spacing to increase the difficulty of exploring the parameter space.

```@example example1
# partition boundaries
bounds = fill((0, 1), n_dims)
p_bounds = range(0, 1, n_part + 1)
hypercube = HyperCube(p_bounds)
```

## Starting Points
The sampling algorim requires a starting point to begin exploring the parameter space. The starting points must be wrapped in a vector. The starting points are sampled uniformly within the unit cube, using `bounds` to ensure the starting point is within allowable ranges. Although one starting point is sufficient for the present example, seeding the sampling algorithm with multiple starting points can improve performance. 

```@example example1
# number of starting points
n_start = 1
# sample function
sample(bounds) = map(b -> rand(Uniform(b...)), bounds)
# initial starting points
init_parms = map(_ -> sample(bounds), 1:n_start)
parm_names = [:p1,:p2,:p3]
```

## Option Configuration

We can configure options to define and improve the performance of the algorithm. The search radius is an important configuration. The challenge is to balance the tradeoff between exploration and exploitation. If the radius is too small, it will stay in one region (or a sub-space of a region), and will fail to find new regions. By contrast, if the radius is too large, many regions will be found, but will not be well defined. Pitt et al. noted that an acceptance rate of 20% may work well in many cases, but this is a heuristic rather than a hard requirement. The options also stores the bounds and initial parameters. Threading can be enabled by setting `parallel=true`. Some exploration revealed that threading becomes advantageous once the evaluation time reaches 1 millisecond or longer. Otherwise, threading overhead will reduce the speed of the algorithm. 

```@example example1
options = Options(;
    radius = .10,
    bounds,
    n_iters = 1000,
    init_parms
)
```
It is also possible to pass a custom adaption function via the keyword `adapt_radius!`. By default, the adaption function adjusts the radius to achieve a 40% acceptance rate. Additional information for configuring the default adaption function can be found via the help feature:

```julia
? adapt!
```

## Find Partitions 

Now that the requisite functions and options have been specified, we can now explore the parameter space.
The function `find_partitions` accepts the `model` function, the pattern function `p_fun`, the options object, and additional arguments and keyword arguments for `p_fun`.

```@example example1
results = find_partitions(
    model, 
    p_fun, 
    options,
    hypercube,
)
first(results)
```

`results` is a vector of `Results` objects containing the following information:

- `iter`: the iteration of the algorithm
- `chain_id`: an index of the chain, i.e. 2 is the second chain
- `parms`: a vector of parameters
- `pattern`: the target pattern of the chain
- `acceptance`: a vector indicating whether a proposal was accepted. If accepted, `parms` is the accepted proposal. If not accepted, `parms` is the same as the previous iteration.

## Results

The next code block finds the minimum and maximum of each partition.

```@example example1
df = DataFrame(results)
groups = groupby(df, :pattern)
summary = combine(groups, 
    :p1 => minimum, :p1 => maximum, 
    :p2 => minimum, :p2 => maximum,
    :p3 => minimum, :p3 => maximum
) |> sort
```
## Volume Estimation

The function `estimate_volume` approximates the volume of each region using an ellipsoid based
on the covariance of sampled points in the region. As expected, the volume percentage estimates are close to 1/8 = .125.

```@example example1
groups = groupby(df, :pattern)

df_volume = combine(
    groups,
    x -> estimate_volume(
        model,
        p_fun, 
        x,  
        bounds,
        hypercube; 
        parm_names,
    )
)

df_volume.volume = df_volume.x1 / sum(df_volume.x1)
```