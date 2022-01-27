function volume_sim(config)
    (;n_shapes,n_sim,n_dims) = config
    (;n_dims,n_part,n_cells, n_start) = config
    # partition boundaries
    bounds = fill((0,1), n_dims)
    indices = set_indices(n_dims, n_part, n_cells)
    p_bounds = range(0, 1, length = n_part + 1)

    odd_shape = OddShape(indices, p_bounds, n_cells)

    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

    init_parms = map(_ -> sample(bounds), 1:n_start)

    options = Options(;
        radius = .10,
        bounds,
        n_iters = 20_000,
        init_parms
    )

    results = find_partitions(
        model, 
        p_fun,
        options,
        odd_shape
    )

    df = DataFrame(results)

    parm_names = Symbol.("p",1:n_dims)
    transform!(
        df, 
        :parms => identity => parm_names
    )

    transform!(df, :pattern => denserank => :pattern_id)

    target_pattern = 1
    df_p = filter(x-> x.pattern == target_pattern, df)
    points = Array(df_p[:,parm_names])
    v1 = estimate_volume(
        model,
        p_fun, 
        points, 
        target_pattern, 
        bounds, 
        odd_shape;
        n_sim,
    )

    target_pattern = 2
    df_p = filter(x-> x.pattern == target_pattern, df)
    points = Array(df_p[:,parm_names])
    v2 = estimate_volume(
        model,
        p_fun, 
        points, 
        target_pattern, 
        bounds, 
        odd_shape;
        n_sim
    )
    return v1 / v2
end