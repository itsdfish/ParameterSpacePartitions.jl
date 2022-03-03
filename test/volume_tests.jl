@safetestset "2D Ellipsoid Volume" begin
    using Test, QHull, Random, Distributions, LinearAlgebra
    using ParameterSpacePartitions: volume_ellipsoid, sample_ellipsoid_surface
    Random.seed!(633)
    n = 2
    μ = fill(0.0, n)
    for _ in 1:5
        x = randn(n, n)
        cov_mat = x' * x
        fun = sample_ellipsoid_surface
        x = map(_ -> fun(μ, n, cov_mat), 1:10_000)
        points = reduce(vcat, transpose.(x))
        
        v1 = volume_ellipsoid(n, cov_mat)
        
        hull = chull(points)
        v2 = hull.volume 
        
        @test v1 ≈ v2 rtol = .03
    end
end

@safetestset "3D Ellipsoid Volume" begin
    using Test, QHull, Random, Distributions, LinearAlgebra
    using ParameterSpacePartitions: volume_ellipsoid, sample_ellipsoid_surface
    Random.seed!(584)
    n = 3
    μ = fill(0.0, n)
    for _ in 1:5
        x = randn(n, n)
        cov_mat = x' * x
        fun = sample_ellipsoid_surface
        x = map(_ -> fun(μ, n, cov_mat), 1:10_000)
        points = reduce(vcat, transpose.(x))
        
        v1 = volume_ellipsoid(n, cov_mat)
        
        hull = chull(points)
        v2 = hull.volume 
        
        @test v1 ≈ v2 rtol = .03
    end
end

@safetestset "4D Ellipsoid Volume" begin
    using Test, QHull, Random, Distributions, LinearAlgebra
    using ParameterSpacePartitions: volume_ellipsoid, sample_ellipsoid_surface
    Random.seed!(28805)
    n = 4
    μ = fill(0.0, n)
    
    x = randn(n, n)
    cov_mat = x' * x
    fun = sample_ellipsoid_surface
    x = map(_ -> fun(μ, n, cov_mat), 1:10_000)
    points = reduce(vcat, transpose.(x))
    
    v1 = volume_ellipsoid(n, cov_mat)
    
    hull = chull(points)
    v2 = hull.volume 
    
    @test v1 ≈ v2 rtol = .03
end

@safetestset "3D Odd Shape Volume" begin

    using Test, Distributions, ParameterSpacePartitions
    using ParameterSpacePartitions.TestModels
    using Random, DataFrames, StatsBase
    include("volume_functions.jl")

    Random.seed!(1495)

    c = (
        # number of shapes
        n_shapes = 2,
        # number of simulations for hit or miss 
        n_sim = 10_000,
        # number of dimensions 
        n_dims = 3,
        # number of partitions per dimension
        n_part = 5,
        # number of cells for each shapes
        n_cells = [10,20],
        # number of starting points for the algorithm
        n_start = 1,
    )

    ratios = map(_ -> volume_sim(c), 1:10)

    tv1 = c.n_cells[1] * (1 / c.n_part)^c.n_dims
    tv2 = c.n_cells[2] * (1 / c.n_part)^c.n_dims
    true_ratio = tv1 / tv2
    ratio = mean(ratios)
    @test ratio ≈ true_ratio rtol = .2
end

@safetestset "5D Odd Shape Volume" begin

    using Test, Distributions, ParameterSpacePartitions
    using ParameterSpacePartitions.TestModels
    using Random, DataFrames, StatsBase
    include("volume_functions.jl")

    Random.seed!(82544)

    c = (
        # number of shapes
        n_shapes = 2,
        # number of simulations for hit or miss 
        n_sim = 10_000,
        # number of dimensions 
        n_dims = 5,
        # number of partitions per dimension
        n_part = 5,
        # number of cells for each shapes
        n_cells = [10,20],
        # number of starting points for the algorithm
        n_start = 1,
    )

    ratios = map(_ -> volume_sim(c), 1:10)

    tv1 = c.n_cells[1] * (1 / c.n_part)^c.n_dims
    tv2 = c.n_cells[2] * (1 / c.n_part)^c.n_dims
    true_ratio = tv1 / tv2
    ratio = mean(ratios)
    @test ratio ≈ true_ratio rtol = .5
end

@safetestset "Volume 5D Polytope" begin 
 
    using Test, ParameterSpacePartitions
    using ParameterSpacePartitions.TestModels
    using Random, DataFrames, Distributions
    Random.seed!(8492)

    # dimensions of the hypbercue
    n_dims = 5
    # number of partitions
    n_part = 5
    # number of starting points
    n_start = 1

    #points = [round.(rand(n_dims), digits=3) for _ in 1:n_part]

    points = [ 
        [0.659, 0.547, 0.842, 0.349, 0.667],
        [0.121, 0.758, 0.141, 0.4, 0.955],
        [0.232, 0.696, 0.999, 0.655, 0.257],
        [0.312, 0.057, 0.977, 0.025, 0.154],
        [0.191, 0.782, 0.051, 0.697, 0.653],
    ]

    # partition boundaries
    polytopes = Polytope.(points)
    bounds = fill((0, 1), n_dims)

    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

    init_parms = map(_ -> sample(bounds), 1:n_start)

    options = Options(;
        radius = .1,
        bounds,
        n_iters = 10_000,
        parallel = false,
        init_parms,
        λ = .025,
        t_rate = .3,
    )

    results = find_partitions(
        model, 
        p_fun, 
        options,
        polytopes
    )

    parm_names = Symbol.("p", 1:n_dims)
    df = DataFrame(results)
    transform!(
        df, 
        :parms => identity => parm_names
    )

    groups = groupby(df, :chain_id)
    mean_accept = combine(groups, :acceptance=>mean=>:mean)

    groups = groupby(df, :pattern)

    df_volume = combine(
        groups,
        x -> estimate_volume(
            model,
            p_fun, 
            x,  
            bounds,
            polytopes; 
            parm_names,
        )
    )

    df_volume.volume = df_volume.x1 / sum(df_volume.x1)
    std(df_volume.volume)
    # based on 10M SMC samples
    true_volume = [
        0.40331
        0.10566
        0.17438
        0.0848
        0.23185
    ]
    true_df = DataFrame(pattern = 1:n_part, volume = true_volume)
    abs_error = abs.(true_volume .- df_volume.volume)
    max_abs_error = maximum(abs_error)
    mean_abs_error = mean(abs_error)
    sort!(df_volume, :volume)
    sort!(true_df, :volume)
    r = cor(df_volume.pattern, true_df.pattern)
    @test max_abs_error < .03
    @test mean_abs_error < .015
    @test r > .8
end