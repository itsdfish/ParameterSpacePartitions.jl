"""
    find_partitions(model, p_fun, options, args...; kwargs...)

Performs parameter space partitioning.

# Arguments

- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `options`: a set of options for configuring the algorithm
- `args...`: arguments passed to `model` and `p_fun`

# Keywords 

- `kwargs...`: keyword arguments passed to `model` and `p_fun`
"""
function find_partitions(model, p_fun, options, args...; kwargs...)
    _model = x -> model(x, args...; kwargs...)
    _p_fun = x -> p_fun(x, args...; kwargs...)
    # evaluate data pattern for each proposal
    patterns = options.p_eval(options.init_parms, _model, _p_fun)
    # initialize
    chains = initialize(options.init_parms, patterns, options)
    # add result type
    results = map(c -> Results(c..., 0), enumerate(chains))
    # list of all observed patterns
    all_patterns = unique(patterns)
    for iter in 1:options.n_iters
        # generate proposal for each chain
        proposals = map(c -> generate_proposal(c, options), chains)
        # evaluate data pattern for each proposal
        patterns = options.p_eval(proposals, _model, _p_fun, options, chains)
        # accept or reject proposals 
        update_position!.(chains, proposals, patterns, options)
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
    for i in 1:length(chains)
        push!(results, Results(i, chains[i], iter))
    end
    return nothing
end

"""
    t_eval_patterns(proposals, model, p_fun, options)

Uses threading to generate patterns associated with a vector of proposals

# Arguments

- `proposals`: a vector of proposals to evaluate
- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `options`: a set of options for configuring the algorithm

"""
function t_eval_patterns(proposals, model, p_fun, options, patterns)
    return tmap(p -> p_fun(model(p)), proposals)
end

"""
    eval_patterns(proposals, model, p_fun, options)

Generates patterns associated with a vector of proposals

# Arguments

- `proposals`: a vector of proposals to evaluate
- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `options`: a set of options for configuring the algorithm

"""
function eval_patterns(proposals, model, p_fun, options, chains)
    (;bounds) = options
    return map((c,p) -> eval_pattern(p_fun, model, bounds, c, p), chains, proposals)
end

function eval_patterns(proposals, model, p_fun)
    return map(p -> p_fun(model(p)), proposals)
end

function eval_pattern(p_fun, model, bounds, chain, parms)
    return in_bounds(parms, bounds) ? p_fun(model(parms)) : chain.pattern
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
    (;x_range,n_dims) = options
    Δ = random_position(chain) * rand() ^(1 / n_dims)
    new_parms = chain.parms + Δ .* options.x_range
    return new_parms 
end

function in_bounds(parms::AbstractArray, bounds)
    for i in 1:length(bounds)
        !in_bounds(parms[i], bounds[i]) ? (return false) : nothing
    end
    return true
end

in_bounds(p::Number, b) = p ≥ b[1] && p ≤ b[2] 

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
    # option for unique
    return Chain.(init_parms, patterns, options.radius)
end

"""
    update_position!(chain, proposal, pattern)

Updates the position of the chain if proposal is accepted 

# Arguments

- `chain`: a chain object for exploring the parameter space 
- `proposal`: a proposed set of parameters for next location 
- `pattern`: a data pattern associated with the proposal
"""
function update_position!(chain, proposal, pattern, bounds)
    if in_bounds(proposal, bounds) && pattern == chain.pattern 
        push!(chain.acceptance, true)
        chain.parms = proposal
    else
        push!(chain.acceptance, false)
    end
    return nothing
end

function update_position!(chain, proposal, pattern, options::Options) 
    return update_position!(chain, proposal, pattern, options.bounds)
end

function process_new_patterns!(all_patterns, patterns, parms, chains, options)
    for p in 1:length(patterns) 
        if !chains[p].acceptance[end] && is_new(all_patterns, patterns[p]) && 
            in_bounds(parms[p], options.bounds)
            push!(all_patterns, patterns[p])
            push!(chains, Chain(parms[p], patterns[p], options.radius))
        end
    end
    return nothing
end

is_new(all_patterns, pattern) = pattern ∉ all_patterns

function no_adaption!(chain, options; kwargs...)
    return nothing
end

"""
    adapt!(
        chain, 
        options; 
        t_rate = .25, 
        λ = .025, 
        trace_on = false,
        max_past = 300, 
        kwargs...
    )

Iteratively adapts the radius to achieve a target acceptance rate. The radius is adjusted according
to the following factor `c`:

```julia
c = exp(λ * d_rate)
```
where `λ` is the adaption rate, and `d_rate` is the difference between the acceptance rate and target 
acceptance rate.

# Arguments

- `chain`: a chain for exploring the parameter space
- `options`: a set of options for configuring the algorithm

# Keyword Arguments

- `t_rate = .25`: target acceptance rate 
- `λ = .025`: adaption rate 
- `trace_on = false`: prints adaption information if true
- `max_past = 300`: maximum past acceptance values considered in adaption
- `kwargs...`: keyword arguments that are not processed
"""
function adapt!(
        chain, 
        options; 
        t_rate = .25, 
        λ = .025, 
        trace_on = false,
        max_past = 300, 
        kwargs...
    )
    n_trials = length(chain.acceptance)
    # evaluate up to max_past previous acceptance values
    start_idx = max(n_trials - max_past, 1)
    # acceptance rate
    a_rate = mean(@view chain.acceptance[start_idx:end])
    # difference between acceptance and target rate
    d_rate = a_rate - t_rate
    # adaption factor 
    c = exp(λ * d_rate)
    # adapt radius  
    chain.radius = min(1, chain.radius * c)
    # print trace
    trace_on ? print_adapt(chain, d_rate, c) : nothing
    return nothing 
end

function print_adapt(chain, d_rate, c)
    println(
        "pattern: ", chain.pattern, 
        " radius: ", round(chain.radius, digits=3), 
        " d rate ", round(d_rate, digits = 3),
        " adjustment factor ", round(c, digits = 3),
    )
    println(" ")
    return nothing
end