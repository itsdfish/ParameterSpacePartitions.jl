using DataFrames

"""
    estimate_volume(
        model, 
        p_fun, 
        df::SubDataFrame, 
        bounds,  
        args...;
        n_sim = 10_000,
        parm_names,
        show_progress = false,
        meter = nothing,
        kwargs...
    )

Estimate volume of region with an eillipsoid and hit or miss bias adjustment. 

# Arguments

- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `df::SubDataFrame`: a SubDataFrame containing MCMC samples in columns identities by `parm_names`
    and column containing data patterns named `pattern`
- `bounds`: a vector of tuples representing (lowerbound, upperbound) for each dimension in 
- `args...`: arguments passed to `model` and `p_fun`

# Keywords 

- `n_sim=10_000`: the number of samples for hit or miss bias adjustment 
- `parm_names`: a vector of symbols representing parameter columns in `df`
- `show_progress`: shows progress meter if true 
- `meter=nothing`: progress meter
- `kwargs...`: additional keyword arguments passed to `model` or `p_fun`
"""
function estimate_volume(
    model, 
    p_fun, 
    df::SubDataFrame, 
    bounds,  
    args...;
    n_sim = 10_000,
    parm_names,
    show_progress = false,
    meter = nothing,
    kwargs...
    )
    
    points = Array(df[:,parm_names])
    target = df.pattern[1]

    df_volume = estimate_volume(
        model, 
        p_fun, 
        points, 
        target, 
        bounds, 
        args...;
        parm_names,
        n_sim,
        kwargs...
    )

    show_progress ? next!(meter) : nothing

    return df_volume
end

"""
    estimate_volume(
        model, 
        p_fun, 
        groups::GroupedDataFrame, 
        bounds,  
        args...;
        n_sim = 10_000,
        show_progress = true,
        parm_names,
        kwargs...
    )
    
    

Estimate volume of region with an eillipsoid and hit or miss bias adjustment. 

# Arguments

- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `df::SubDataFrame`: a SubDataFrame containing MCMC samples in columns identities by `parm_names`
    and column containing data patterns named `pattern`
- `bounds`: a vector of tuples representing (lowerbound, upperbound) for each dimension in 
- `args...`: arguments passed to `model` and `p_fun`

# Keywords 

- `n_sim=10_000`: the number of samples for hit or miss bias adjustment 
- `parm_names`: a vector of symbols representing parameter columns in `df`
- `show_progress`: shows progress meter if true 
- `meter=nothing`: progress meter
- `kwargs...`: additional keyword arguments passed to `model` or `p_fun`
"""
function estimate_volume(
    model, 
    p_fun, 
    groups::GroupedDataFrame, 
    bounds,  
    args...;
    n_sim = 10_000,
    show_progress = false,
    parm_names,
    kwargs...
    )
    
    n_groups = length(groups)
    meter = Progress(n_groups)

    return combine(
        groups,
        x -> estimate_volume(
            model,
            p_fun, 
            x,  
            bounds,
            args...;
            parm_names,
            meter,
            show_progress,
            n_sim,
            kwargs...
        )
    )
end