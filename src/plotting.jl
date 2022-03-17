using StatsPlots 

"""
    psp_slices(
        model, 
        p_fun, 
        parm_names, 
        bounds, 
        v_parm; 
        v_range,
        sample,
        name_fun = x -> x,
        opt_set = (),  
        kwargs...
    )


Returns a grid of scatter plots which represent slices of the parameter space at values of
a third value `v_parm`. The expected signiture for `model` is 

    function model(parms::ComponentVector; p1=.3, p2=.2, kwargs...)
        main_model(; p1, p2, parms..., kwargs...)
    end

    function main_model(;p1, p2, p3, p4)
        return main_model(p1, p1, p2, p4)
    end

    function main_model(p1, p1, p2, p4)
        ...
    end

# Arguments 

- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern
- `parm_names`: a vector of symbols corresponding to parameter names  
- `bounds`: a vector of tuples representing (lowerbound, upperbound) for each dimension in 
the parameter space

# Keywords 

- `v_parm`: the name (symbol) of the parameter which is varied to generate subplots
- `v_range`: the range of values over which `v_parm` is varied 
- `name_fun x -> x`: a function that processes/formats the patterns. The default function is the 
identiy function 
- `sample`: a function in the form of `sample(bounds, parm_names)` which samples an initial point
- `opt_set=()`: a `NamedTuple` of settings to be based to `Options`
- `kwargs...`: optional keyword arguments for `model` and `p_fun`

"""
function psp_slices(
    model, 
    p_fun, 
    parm_names, 
    bounds, 
    v_parm; 
    v_range,
    sample,
    name_fun = x->x,
    opt_set = (),
    plot_options = (),  
    kwargs...
    )
    plots = map(
        x -> psp_slice(
            model, 
            p_fun, 
            parm_names, 
            bounds, 
            v_parm, 
            x;
            sample,
            name_fun,
            opt_set, 
            plot_options,
            kwargs...
        ),
        v_range
    )
    return plot(plots...)
end

function psp_slice(
    model, 
    p_fun, 
    parm_names, 
    bounds, 
    v_parm, 
    v_val;
    sample,
    name_fun = x->x,
    opt_set = (), 
    plot_options = (),
    kwargs...
    ) 

    init_parms = map(_ -> sample(bounds, parm_names), 1:1)

    options = Options(;
        radius = .10,
        bounds,
        n_iters = 5_000,
        init_parms,
        parm_names,
        opt_set...
    )

    df = find_partitions(
        model, 
        p_fun,
        options;
        kwargs...,
        v_parm => v_val
    )

    transform!(df, :pattern => name_fun => :pattern_var)

    p = @df df scatter(
        cols(parm_names[1]), 
        cols(parm_names[2]),
        group = :pattern_var,
        grid = false,
        legend = :outertopright,
        xaxis = string(parm_names[1]),
        yaxis = string(parm_names[2]),
        title = string(v_parm, " = ", round(v_val, digits=3)),
        titlefontsize = 10,
        palette = :tab10,
        plot_options...
    )
    return p
end
