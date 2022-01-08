@concrete struct Options
    parallel
    radius
    bounds
    n_iters
    p_eval
    init_parms
end

function Options(;
    radius, 
    bounds, 
    n_iters, 
    init_parms,
    parallel = false
)
    p_eval = parallel ? t_eval_patterns : eval_patterns
    
    return Options(
        parallel,
         radius, 
         bounds, 
         n_iters, 
         p_eval, 
         init_parms
    )
end

"""
    Chain(;parms, pattern, radius) 

# Keywords
"""
@concrete mutable struct Chain 
    parms
    n_dims
    pattern
    radius
    acceptance 
end

function Chain(parms, pattern, radius) 
    return Chain(
        parms,
        length(parms),
        pattern, 
        radius,
        Bool[]
    )
end

@concrete struct Results
    iter
    chain_id
    parms 
    pattern
end

function Results(chain_id, chain, iter)
    return Results(
        iter, 
        chain_id, 
        chain.parms, 
        chain.pattern
    )
end
