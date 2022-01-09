"""
    find_partitions(model, p_fun, options)

Performs parameter space partitioning.

# Arguments

- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `options`: a set of options for configuring the algorithm
"""
function find_partitions(model, p_fun, options)
    # evaluate data pattern for each proposal
    patterns = options.p_eval(options.init_parms, model, p_fun, options)
    # initialize
    chains = initialize(options.init_parms, patterns, options)
    # add result type
    results = map(c -> Results(c..., 0), enumerate(chains))
    all_patterns = unique(patterns)
    for iter in 1:options.n_iters
        # generate proposal for each chain
        proposals = map(c -> generate_proposal(c, options), chains)
        # evaluate data pattern for each proposal
        patterns = options.p_eval(proposals, model, p_fun, options)
        # accept or reject proposals 
        update_position!.(chains, proposals, patterns)
        # record accepted results and update chains
        update_results!(results, chains, iter)
        # add new chain if new pattern found
        process_new_patterns!(all_patterns, patterns, proposals, chains, options)
        # adjust the radius of each chain 
        options.adapt_radius!.(chains, options)
    end
    return results
end

"""
    update_results!(results, chains, iter)

Adds chain location (parameters), iteration, chain id, and data pattern to a `Results` vector.

# Arguments 

- `results`: a vector of `Results` objects
- `chains`: a vector of chains used to explore the parameter space 
- `iter`: current iteration of the algorithm
"""
function update_results!(results, chains, iter)
    for (i,c) in enumerate(chains)
        push!(results, Results(i, c, iter))
    end
    return nothing
end

"""
"""
function t_eval_patterns(proposals, model, p_fun, options)
        return tmap(p -> p_fun(model(p)), proposals)
end

function eval_patterns(proposals, model, p_fun, options)
    return map(p -> p_fun(model(p)), proposals)
end

"""
    generate_proposal(chain::Chain, options)

Generates a proposal adding a random sample from the surface of a hypersphere
to the current location of the `chain`.

# Arguments

- `chain`: a chain for exploring the parameter space
- `options`: an `Options` object holding the configuration options for the algorithm
"""
function generate_proposal(chain::Chain, options)
    Δ = random_position(chain)
    new_parms = chain.parms + Δ
    adjust_parms!(new_parms, options)
    return new_parms 
end

function generate_proposal(parms, options)
    Δ = random_position(options.radius, length(parms))
    new_parms = parms + Δ
    adjust_parms!(new_parms, options)
    return new_parms 
end

adjust_parms!(new_parms, options::Options) = adjust_parms!(new_parms, options.bounds)

function adjust_parms!(parms, bounds)
    map!((p,b) -> adjust_parm(p, b), parms, parms, bounds)
    return nothing
end

adjust_parm(p, b) = min(max(p, b[1]), b[2])

function random_position(chain) 
    return random_position(chain.radius, chain.n_dims)
end

"""
    random_position(radius, n)

A sample from a uniform distribution over the surface of a hypersphere 

# Arguments

- `radius`: radius of the hypersphere 
- `n`: the number of dimentions of the hypersphere
"""
function random_position(radius, n)
    x = rand(Normal(0, 1), n)
    return x / norm(x) * radius
end

function initialize(init_parms, patterns, options)
    return Chain.(init_parms, patterns, options.radius)
end

function update_position!(chain, proposal, pattern)
    if pattern == chain.pattern 
        push!(chain.acceptance, true)
        chain.parms = proposal
    else
        push!(chain.acceptance, false)
    end
    return nothing
end

function process_new_patterns!(all_patterns, patterns, parms, chains, options)
    for p in 1:length(patterns) 
        if !chains[p].acceptance[end] && is_new(all_patterns, patterns[p])
            push!(all_patterns, patterns[p])
            push!(chains, Chain(parms[p], patterns[p], options.radius))
        end
    end
    return nothing
end

is_new(all_patterns, pattern) = pattern ∉ all_patterns

function no_adaption(chain, options)
    return nothing
end