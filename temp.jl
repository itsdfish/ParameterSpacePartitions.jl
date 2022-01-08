cd(@__DIR__)
using Pkg
Pkg.activate("")
using Revise, Distributions, ParameterSpacePartitions
using Random

n_dims = 3
bounds = fill((0, 1), n_dims)

radius = .2

model(parms) = parms 

# partition boundaries
p_bounds = range(0, 1, 3)

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

sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

init_parms = map(_ -> sample(bounds), 1:1)

options = Options(;
    radius,
    bounds,
    n_iters = 100,
    parallel = false,
    init_parms
)

results = find_partitions(
    model, 
    x -> p_fun(x, p_bounds), 
    options
)

df = DataFrame(results)
transform!(df, :parms => identity => [:p1, :p2,:p3])

@df df scatter(:p1, :p2, :p3, group=:pattern)

groups = groupby(df, :pattern)
combine(groups, 
    :p1 =>minimum, :p1=>maximum, 
    :p2=>minimum, :p2=>maximum,
    :p3=>minimum, :p3=>maximum
) |> sort
