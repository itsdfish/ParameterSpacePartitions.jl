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

    temp_chains = map(
        p -> find_regions(
            model, 
            p_fun, 
            options, 
            p,
            args...;
            kwargs...
        ),
        options.init_parms
    )

    chains = vcat(temp_chains...)
    n_init = length(options.init_parms)
    n_init > 1 ? make_unique!(chains, options) : nothing
    n_init > 1 ? make_unique!(chains, options) : nothing
    patterns = map(c -> c.pattern, chains)
    all_patterns = unique(patterns)
    for iter in 1:options.add_iters
        # generate proposal for each chain
        proposals = map(c -> propose(c, options), chains)
        # evaluate data pattern for each proposal
        patterns = options.p_eval(proposals, _model, _p_fun, options, chains)
        # accept or reject proposals 
        update_position!.(chains, proposals, patterns, options)
        # add new chain if new pattern found
        process_new_patterns!(all_patterns, patterns, proposals, chains, options)
        # adjust the radius of each chain 
        options.adapt_radius!.(chains, options)
    end
    n_init > 1 ? make_unique!(chains, options) : nothing
    n_init > 1 ? make_unique!(chains, options) : nothing
    return to_df(chains, options)
end

function find_regions(model, p_fun, options, init_parm, args...; kwargs...)
    _model = x -> model(x, args...; kwargs...)
    _p_fun = x -> p_fun(x, args...; kwargs...)

    all_patterns = options.p_eval([init_parm], _model, _p_fun)
    chains = initialize([init_parm], all_patterns, options)
    #options.n_ == 0 ? (return chains) : nothing
    complete_chains = Vector{eltype(chains)}()
    
    while !isempty(chains)
        # generate proposal for each chain
        proposals = map(c -> propose(c, options), chains)
        # evaluate data pattern for each proposal
        patterns = options.p_eval(proposals, _model, _p_fun, options, chains)
        # accept or reject proposals 
        update_position!.(chains, proposals, patterns, options)
        # add new chain if new pattern found
        process_new_patterns!(all_patterns, patterns, proposals, chains, options)
        # adjust the radius of each chain 
        options.adapt_radius!.(chains, options)
        remove_complete!(complete_chains, chains, options)
    end
    return complete_chains 
end

function remove_complete!(complete_chains, chains, options)
    remove = fill(false, length(chains))
    for (i,c) in enumerate(chains) 
        if is_complete(c, options)
            push!(complete_chains, c)
            remove[i] = true            
        end
    end
    deleteat!(chains, remove)
    return nothing  
end

function is_complete(chain, options)
    if (chain.level == 2) && (chain.n_attempt == options.n_iters)
        return true
    end
    return false 
end 

"""
    t_eval_patterns(proposals, model, p_fun, options, chains)

Uses threading to generate patterns associated with a vector of proposals

# Arguments

- `proposals`: a vector of proposals to evaluate
- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `options`: a set of options for configuring the algorithm
- `chains`: a vector of `Chain` objects
"""
function t_eval_patterns(proposals, model, p_fun, options, chains)
    (;bounds) = options
    return tmap((c,p) -> eval_pattern(p_fun, model, bounds, c, p), chains, proposals)
end

function t_eval_patterns(proposals, model, p_fun)
    return tmap(p -> p_fun(model(p)), proposals)
end

"""
    eval_patterns(proposals, model, p_fun, options, chains)

Generates patterns associated with a vector of proposals

# Arguments

- `proposals`: a vector of proposals to evaluate
- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `options`: a set of options for configuring the algorithm
- `chains`: a vector of `Chain` objects

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
    propose(chain::Chain, options)

Generates a proposal adding a random sample from the surface of a hypersphere
to the current location of the `chain`.

# Arguments

- `chain`: a chain for exploring the parameter space
- `options`: an `Options` object holding the configuration options for the algorithm
"""
function propose(chain::Chain, options)
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
    options.last_id += 1
    last_id = options.last_id
    n_start = length(init_parms)
    ids = last_id:(last_id + n_start - 1)
    options.last_id += (n_start - 1)
    return Chain.(ids, init_parms, patterns, options.radius)
end

"""
    update_position!(chain, proposal, pattern, bounds)

Updates the position of the chain if proposal is accepted 

# Arguments

- `chain`: a chain object for exploring the parameter space 
- `proposal`: a proposed set of parameters for next location 
- `pattern`: a data pattern associated with the proposal
- `bounds`: a vector of tuples for lower and upper bounds of each parameter
"""
function update_position!(chain, proposal, pattern, bounds)
    chain.n_attempt += 1
    if in_bounds(proposal, bounds) && pattern == chain.pattern 
        push!(chain.acceptance, true)
        chain.n_accept += 1
        chain.parms = proposal
        push!(chain.radii, chain.radius)
        push!(chain.all_parms, chain.parms)
    else
        push!(chain.all_parms, chain.parms)
        push!(chain.radii, chain.radius)
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
            options.last_id += 1
            push!(chains, Chain(options.last_id, parms[p], patterns[p], options.radius))
        end
    end
    return nothing
end

is_new(all_patterns, pattern) = pattern ∉ all_patterns

function no_adaption!(chain, options; kwargs...)
    return nothing
end

"""
    exp_adapt!(
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
- `λ = .05`: adaption rate 
- `trace_on = false`: prints adaption information if true
- `max_past = 300`: maximum past acceptance values considered in adaption
- `kwargs...`: keyword arguments that are not processed
"""
function exp_adapt!(
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

"""
    adapt!(
        chain, 
        options; 
        t_rate = .20, 
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

- `t_rate = .20`: target acceptance rate 
- `kwargs...`: keyword arguments that are not processed
"""
function adapt!(
    chain, 
    options; 
    t_rate = .20, 
    kwargs...
)
    chain.level == 2 ? (return nothing) : nothing
    n_dims = length(options.bounds)
    n_attempt = chain.n_attempt  
    λ = chain.λ
    if chain.level == 0
        n = ceil(100 * 1.2^n_dims)
        if mod(n_attempt, n) == 0
            a_rate = chain.n_accept / n
            chain.n_accept = 0
            if a_rate < (t_rate - .08)
                if λ > 0
                    λ -= .5
                    chain.level = 1
                    chain.n_attempt = 0
                else
                    λ -= 1.0
                end
            elseif (a_rate ≥ (t_rate - .08)) && (a_rate < (t_rate + .16))
                chain.level = 1
                chain.n_attempt = 0
            elseif a_rate ≥ (t_rate + .16)
                if  λ < 0 
                    λ += .5
                    chain.level = 1
                    chain.n_attempt = 0
                else
                    λ += 1
                end 
            end
        end
    elseif chain.level == 1
        n = ceil(200 * 1.2^n_dims)
        v = n_attempt / n
        if mod(n_attempt, n) == 0
            a_rate = chain.n_accept / n
            chain.n_accept = 0
            if a_rate < (t_rate - .05)
                λ = λ - .25 / ceil(v / 2)
                if v == 4
                    chain.level = 2
                    chain.n_attempt = 0
                end
            elseif (a_rate ≥ (t_rate - .05)) && (a_rate < (t_rate - .01))
                λ -= .125
                chain.level = 2
                chain.n_attempt = 0
            elseif (a_rate ≥ (t_rate - .01)) && (a_rate < (t_rate + .04))
                chain.level = 2
                chain.n_attempt = 0
            elseif (a_rate ≥ (t_rate + .04)) && (a_rate < (t_rate + .10))
                λ += .125
                chain.level = 2
                chain.n_attempt = 0
            elseif a_rate > (t_rate + .10)
                λ = λ + .25 / ceil(v / 2)
                if v == 4
                    chain.level = 2
                    chain.n_attempt = 0
                end
            end
        end
    end
    chain.radius = options.radius * 2^λ
    chain.λ = λ
    return nothing 
end

function accept_rate(chain, n_trials::Int, n::Int)
    m = div(n_trials, n)
    lb = (m - 1) * n + 1
    ub = m * n 
    return mean(@view chain.accept_rate[lb:ub])
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

function to_df(chains, options)
    df = DataFrame(chains)
    col_names = [:chain_id,:pattern]
    return combine(
        groupby(df, col_names), 
        :all_parms => only => options.parm_names,
        :acceptance => only => :acceptance,
        :radii => only => :radius
    )
end