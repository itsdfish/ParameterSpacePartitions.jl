using DataFrames

function estimate_volume(
    model, 
    p_fun, 
    df::SubDataFrame, 
    bounds,  
    args...;
    n_sim = 10_000,
    parm_names,
    kwargs...
    )
    
    points = Array(df[:,parm_names])
    target = df.pattern[1]

    return estimate_volume(
        model, 
        p_fun, 
        points, 
        target, 
        bounds,  
        args...;
        n_sim,
        kwargs...
    )
end