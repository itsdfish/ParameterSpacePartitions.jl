# ParameterSpacePartitions

Parameter space partitioning is a model assessment method for mapping regions of the parameter space to qualitative data patterns. Parameter space partitioning uses MCMC sampling to explore the parameter space. Each chain searches a region of the parameter space for a specific pattern. As the space is sampled, a new chain is created for each newly discovered pattern. On each iteration, a proposal is sampled uniformly from the surface of a hyperspere with the same number of dimensions as the parameter space. 

# Example

In this simple example, parameter space partitioning is applied to a cube with regions of equal volume.
## Load Dependencies

The first step is to the load dependences into your session.

```julia
using Distributions, ParameterSpacePartitions
using ParameterSpacePartitions.TestModels
using Random, DataFrames
Random.seed!(544)
```

## Model Function 

Next, we define a model that accepts a vector of parameters and returns data predictions. The function `model` is imported above from `TestModels`, but is presented below for illustration. In this simple example, the model returns the parameter inputs. In more substantive model, the returned value will be the model predictions.

```julia 
model(parms, args...; kwargs...) = parms 
```
## Pattern Function

A second function categorizes the predicted data into a qualitative pattern. At minimum, the pattern function must recieve a data input. In this example, the pattern function `p_fun` is imported from `TestModels`. The function is presented below for illustration. `p_fun` requires the location (i.e. parameters) and `HyperCube` object, which contains  partition boundaries (which is the same for each dimension). `p_fun` categorizes `location` on each dimension and returns a vector, such as `[1,2,1]`, which indicates the data is in the partition that is 1 in the x axis, 2 on the y axis and 1 on the z axis. 

```julia
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

```julia
# dimensions of the hypbercue
n_dims = 3
# partitions per dimension
n_part = 2
```

# Boundaries

In the code below, `bounds` contains the upper and lower boundaries for each dimension. In addition, we must also define a `HyperCube` object containing the boundaries for each partition. The boundaries are stored in a variable called `p_bounds`.  In this example, the partitions are equally spaced along each dimension of the unit cube. Alternatively, we could use unequal spacing to increase the difficulty of exploring the parameter space.

```julia
# partition boundaries
bounds = fill((0, 1), n_dims)
p_bounds = range(0, 1, n_part + 1)
hypercube = HyperCube(p_bounds)
```

## Starting Points
The sampling algorim requires a starting point to begin exploring the parameter space. The starting points must be wrapped in a vector. The starting points are sampled uniformly within the unit cube, using `bounds` to ensure the starting point is within allowable ranges. Although one starting point is sufficient for the present example, seeding the sampling algorithm with multiple starting points can improve performance. 

```julia 
# number of starting points
n_start = 1
# sample function
sample(bounds) = map(b -> rand(Uniform(b...)), bounds)
# initial starting points
init_parms = map(_ -> sample(bounds), 1:n_start)
```

## Option Configuration

We can configure options to define and improve the performance of the algorithm. The search radius is an important configuration. The challenge is to balance the tradeoff between exploration and exploitation. If the radius is too small, it will stay in one region (or a sub-space of a region), and will fail to find new regions. By contrast, if the radius is too large, many regions will be found, but will not be well defined. Pitt et al. noted that an acceptance rate of 20% may work well in many cases, but this is a heuristic rather than a hard requirement. The options also stores the bounds and initial parameters. Threading can be enabled by setting `parallel=true`. Some exploration revealed that threading becomes advantageous once the evaluation time reaches 1 millisecond or longer. Otherwise, threading overhead will reduce the speed of the algorithm. 

```julia
options = Options(;
    radius = .40,
    bounds,
    n_iters = 500,
    parallel = false,
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

```julia
results = find_partitions(
    model, 
    p_fun, 
    options,
    hypercube,
)
```

`results` is a vector of `Results` objects containing the following information:

- `iter`: the iteration of the algorithm
- `chain_id`: an index of the chain, i.e. 2 is the second chain
- `parms`: a vector of parameters
- `pattern`: the target pattern of the chain
- `acceptance`: a vector indicating whether a proposal was accepted. If accepted, `parms` is the accepted proposal. If not accepted, `parms` is the same as the previous iteration.

## Results

To facilitate the analysis, we will convert the results to a `DataFrame` and destructure the parameter vector into individual columns---one per parameter. 

```julia
df = DataFrame(results)
transform!(
    df, 
    :parms => identity => [:p1, :p2, :p3]
)
```

The next code block finds the minimum and maximum of each partition.

```julia
groups = groupby(df, :pattern)
summary = combine(groups, 
    :p1 => minimum, :p1 => maximum, 
    :p2 => minimum, :p2 => maximum,
    :p3 => minimum, :p3 => maximum
) |> sort
```
As shown below, the algorithm found all 64 partitions. In addition, the size of the partition is approximately 1/2 = .50, which is what we expect. 

```julia
8×7 DataFrame
 Row │ pattern    p1_minimum  p1_maximum  p2_minimum  p2_maximum  p3_minimum  p3_maximum 
     │ Array…     Float64     Float64     Float64     Float64     Float64     Float64    
─────┼───────────────────────────────────────────────────────────────────────────────────
   1 │ [1, 1, 1]    0.0         0.497357    0.0         0.498463    0.0         0.498483
   2 │ [1, 1, 2]    0.0         0.499921    0.0         0.497982    0.507963    1.0
   3 │ [1, 2, 1]    0.0         0.488883    0.507123    1.0         0.0         0.496436
   4 │ [1, 2, 2]    0.0         0.495029    0.500056    1.0         0.520285    1.0
   5 │ [2, 1, 1]    0.503028    1.0         0.0         0.49685     0.0         0.496929
   6 │ [2, 1, 2]    0.500573    1.0         0.0         0.48709     0.502865    1.0
   7 │ [2, 2, 1]    0.505885    1.0         0.500411    1.0         0.0         0.497511
   8 │ [2, 2, 2]    0.504082    1.0         0.500305    1.0         0.503377    1.0
  ```

## Visualization

The following code shows how to visualize the results. 

```julia
using GLMakie, ColorSchemes, StatsBase

# transform pattern into integer id
transform!(df, :pattern => denserank => :pattern_id)

scatter(
    df.p1,
    df.p2,
    df.p3, 
    color = df.pattern_id, 
    axis = (;type=Axis3),
    markersize = 5000,
    colormap = ColorSchemes.seaborn_deep6.colors,
    grid = true
)
```

<img src="resources/figure.png" alt="" width="500" height="300">

# References

Pitt, M. A., Kim, W., Navarro, D. J., & Myung, J. I. (2006). Global model analysis by parameter space partitioning. Psychological Review, 113(1), 57.
