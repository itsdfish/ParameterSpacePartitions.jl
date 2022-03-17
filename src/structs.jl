"""
    Options(;
        radius = .1, 
        bounds, 
        n_iters, 
        init_parms,
        parallel = false,
        adapt_radius! = adapt!,
        kwargs...
    )

An object that holds configuration options for the MCMC sampler. 

# Fields 

- `parallel=false`: runs model on threads if true. A speed up is observed if the evaluation 
time of the function is 1 ms or greater. 
- `radius = .10`: the initial radius length for each chain 
- `bounds`: a vector of tuples representing (lowerbound, upperbound) for each dimension in 
the parameter space
- `x_range`: the range of allowable values for each parameter
- `n_iters`: number of iterations to perform 
- `p_eval`: the function that evalues the model and pattern functions
- `adapt_radius!=adapt!`: a function in the form of `func(chain, options; kwargs...)` that adapts 
the radius. 
- `init_parms`: a vector of starting points, such as [[.3,.4],[.3,.5]] in the 2 dimensional case.
- `parm_names`: a vector of symbols corresponding to parameter names. The default is [:p1,:p2,..:pn] 
"""
@concrete struct Options
    parallel
    radius
    bounds
    x_range
    n_iters
    p_eval
    adapt_radius!
    init_parms
    n_dims
    parm_names
end

function Options(;
        radius = .10, 
        bounds, 
        n_iters, 
        init_parms,
        parallel = false,
        adapt_radius! = adapt!,
        parm_names = nothing,
        kwargs...
    )

    p_eval = parallel ? t_eval_patterns : eval_patterns
    _adapt_radius! = (x,y) -> adapt_radius!(x, y; kwargs...)
    x_range = map(x -> x[2] - x[1], bounds)
    n_dims = length(bounds)
    _parm_names = isnothing(parm_names) ? Symbol.("p", 1:n_dims) : parm_names

    return Options(
        parallel,
        radius, 
        bounds,
        x_range, 
        n_iters, 
        p_eval, 
        _adapt_radius!,
        init_parms,
        n_dims,
        _parm_names
    )
end

Broadcast.broadcastable(x::Options) = Ref(x)

"""
    Chain(id, parms, pattern, radius)   

An MCMC chain object.

# Arguments

- `id`: chain index
- `parms`: a vector of parameters, i.e the current state of the chain
- `pattern`: the target pattern of the chain 
- `radius`: the radius for the jump distribution

# Fields

- `parms`: a vector of parameters, i.e the current state of the chain
- `n_dims`: the number of dimensions 
- `pattern`: the target pattern of the chain 
- `radius`: the radius for the jump distribution
- `acceptance`: a Boolean vector indicating whether a proposal was accepted
- `init_parms`: a vector of starting points, such as [[.3,.4],[.3,.5]] in the 2 dimensional case. 
"""
@concrete mutable struct Chain
    chain_id 
    parms
    n_dims
    pattern
    radius
    acceptance
    all_parms
    radii 
end

function Chain(id, parms, pattern, radius) 
    return Chain(
        id,
        parms,
        length(parms),
        pattern, 
        radius,
        [true],
        [parms],
        [radius]
    )
end