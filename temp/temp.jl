cd(@__DIR__)
using Pkg
Pkg.activate("")
using Revise, Distributions, ParameterSpacePartitions
using Random, DataFrames, StatsBase
Random.seed!(554)
model(parms) = parms 

function p_fun(data, p_bounds)
    nb = length(p_bounds)
    nd = length(data)
    vals = fill(-100, nd)
    for j in 1:nd
        for i in 1:(nb-1) 
            if (data[j] ≥ p_bounds[i]) && (data[j] ≤ p_bounds[i+1])
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
# number of starting points
n_start = 1

# partition boundaries
p_bounds = range(0, 1, n_part + 1)
bounds = fill((0, 1), n_dims)

sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

init_parms = map(_ -> sample(bounds), 1:n_start)

options = Options(;
    radius = .20,
    bounds,
    n_iters = 500,
    parallel = false,
    init_parms
)

results = find_partitions(
    model, 
    x -> p_fun(x, p_bounds), 
    options
)

df = DataFrame(results)
transform!(
    df, 
    :parms => identity => [:p1, :p2, :p3]
)

transform!(df, :pattern => denserank => :pattern_id)

groups = groupby(df, :pattern)
summary = combine(groups, 
    :p1 => minimum, :p1 => maximum, 
    :p2 => minimum, :p2 => maximum,
    :p3 => minimum, :p3 => maximum
) |> sort


using GLMakie

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



using GLMakie 

points = Observable(Point3f[(df.p1[1], df.p2[1], df.p3[1])])
using ColorSchemes
fig, ax, scat = scatter(
    points, 
    color = df.pattern_id, 
    axis = (;type=Axis3),
    markersize = 5000,
    colormap = ColorSchemes.seaborn_deep6.colors,
    grid = true,
)


limits!(ax, 0, 1, 0, 1, 0, 1)
frames = 2:1000

record(fig, "append_animation.mp4", frames;
        framerate = 100) do frame
    new_point = Point3f(df.p1[frame], df.p2[frame], df.p3[frame])
    points[] = push!(points[], new_point)
end