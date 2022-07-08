
volume_hypersphere(n, r=1) = π^(n / 2) / gamma(1 + n / 2) * r^n

"""
    log_vol_hypersphere(n, r=1)

Computes the log volume of a hypersphere

# Arguments

- `n`: the number of dimensions

# Keywords

- `r=1`: the radius
"""
log_vol_hypersphere(n; r=1) = (n / 2) * log(π) - loggamma(1 + n / 2) + n * log(r)

"""
    log_vol_ellipsiod(n, cov_mat; r = 1)

Computes the log volume of a ellipsoid

# Arguments
- `n`: the number of dimensions
- `cov_mat`: covariance matrix

# Keywords

- `r=1`: the radius
"""
function log_vol_ellipsiod(n, cov_mat; r = 1)
    vol = log_vol_hypersphere(n; r)
    eig_val,_ = eigen(cov_mat)
    return vol + .5 * sum(log(n + 2) .+ log.(eig_val))
end

"""
    volume_ellipsoid(n, cov_mat; r=1)   

Computes the volume of an ellipsoid

# Arguments
- `n`: the number of dimensions
- `cov_mat`: covariance matrix

# Keywords

- `r=1`: the radius
"""
function volume_ellipsoid(n, cov_mat; r = 1)
    return exp(log_vol_ellipsiod(n, cov_mat; r))
end

volume_ellipsoid(cov_mat; r=1) = volume_ellipsoid(size(cov_mat, 2), cov_mat; r)

log_vol_ellipsiod(cov_mat; r=1) = log_vol_ellipsiod(size(cov_mat, 2), cov_mat; r)

"""
    sample_ellipsoid(μ, n_dims, cov_mat)

Samples a point uniformly from the interior of an ellipsoid 

# Arguments

- `μ`: mean of ellipsoid on each dimension  
- `n_dims`: the number of dimensions 
- `cov_mat`: covariance matrix 
"""
function sample_ellipsoid(μ, n_dims, cov_mat)
    x = random_position(rand(), n_dims)
    return μ .+ sqrt((n_dims + 2) * cov_mat) * x
end

"""
    sample_ellipsoid_surface(μ, n_dims, cov_mat; r=1)

Samples a point uniformly from the surface of an ellipsoid 

# Arguments

- `μ`: mean of ellipsoid on each dimension  
- `n_dims`: the number of dimensions 
- `cov_mat`: covariance matrix

# Keywords

- `r=1`: radius of underlying unit hypersphere 
"""
function sample_ellipsoid_surface(μ, n_dims, cov_mat; r=1)
    x = random_position(r, n_dims)
    return μ .+ sqrt((n_dims + 2) * cov_mat) * x
end

"""
    bias_correction(
        model,
        p_fun, 
        points, 
        target, 
        bounds, 
        args...; 
        parm_names = Symbol.("p", 1:size(points,2)),
        n_sim=10_000, 
        kwargs...
    )

Uses the hit or miss technique to adjust volume estimates based on ellipsoids. 

# Arguments 

- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `points`: a p x n matrix of sampled points where p is the number of parameters and n is the number of samples
- `bounds`: a vector of tuples representing (lowerbound, upperbound) for each dimension in 
- `args...`: arguments passed to `model` and `p_fun`

# Keywords

- `parm_names`: parameter names with default `p1`, `p2,`... `pn`
- `n_sim=10_000`: number of simulations for hit or miss adjustment 
- `kwargs...`: optional keyword arguments for `model` and `p_fun`
"""
function bias_correction(
    model,
    p_fun, 
    points, 
    target, 
    bounds, 
    args...; 
    parm_names = Symbol.("p", 1:size(points,2)),
    n_sim = 10_000, 
    kwargs...
    )
    
    _μ = mean(points, dims=1)[:]
    μ = ComponentArray(_μ, make_axis(parm_names))
    cov_mat = cov(points)
    n = length(μ)

    _model = x -> model(x, args...; kwargs...)
    _p_fun = x -> p_fun(x, args...; kwargs...)

    hits = 0
    for i in 1:n_sim
        parms = sample_ellipsoid(μ, n, cov_mat)
        if in_bounds(parms, bounds)   
            pattern = _p_fun(_model(parms))
            hits += hit_or_miss(pattern, target, parms)
        end
    end
    return hits / n_sim
end

function hit_or_miss(target, pattern, parms)
    if target == pattern
        return 1
    end
    return 0
end

"""
    estimate_volume(
        model, 
        p_fun, 
        points, 
        target, 
        bounds,  
        args...;
        n_sim = 10_000,
        kwargs...
    )

Estimate volume of region with an eillipsoid and hit or miss bias adjustment. 

# Arguments
- `model`: a model function that returns predictions given a vector of parameters 
- `p_fun`: a function that that evaluates the qualitative data pattern 
- `points`: a p x n matrix of sampled points where p is the number of parameters and n is the number of samples
- `target`: the target pattern associated with the `points`
- `bounds`: a vector of tuples representing (lowerbound, upperbound) for each dimension in 
- `args...`: arguments passed to `model` and `p_fun`

# Keywords 
- `n_sim=10_000`: the number of samples for hit or miss bias adjustment 
- `kwargs...`: additional keyword arguments passed to `model` or `p_fun`
"""
function estimate_volume(
    model, 
    p_fun, 
    points, 
    target, 
    bounds,  
    args...;
    n_sim = 10_000,
    parm_names = Symbol.("p", 1:size(points,2)),
    kwargs...
    )

    cf = bias_correction(
        model,
        p_fun, 
        points, 
        target, 
        bounds, 
        args...;
        n_sim,
        parm_names,
        kwargs...
    )
    
    cov_mat = cov(points)
    volume = volume_ellipsoid(cov_mat)
    volume *= cf
    return volume
end

make_axis(symbols) = Axis(NamedTuple(symbols .=> eachindex(symbols)))