using SafeTestsets

@safetestset "Cube" begin
    using ParameterSpacePartitions
    using ParameterSpacePartitions.TestModels
    using Test, Random, Distributions
    using DataFrames

    Random.seed!(541)
    # dimensions of the hypbercue
    n_dims = 3
    # partitions per dimension
    n_part = 4
    # number of starting points
    n_start = 1

    # partition boundaries
    p_bounds = range(0, 1, n_part + 1)
    hypercube = HyperCube(p_bounds)
    bounds = fill((0, 1), n_dims)

    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

    init_parms = map(_ -> sample(bounds), 1:n_start)

    options = Options(;
        radius = .18,
        bounds,
        n_iters = 100,
        parallel = false,
        init_parms
    )

    results = find_partitions(
        model, 
        p_fun, 
        options,
        hypercube
    )

    df = DataFrame(results)
    transform!(
        df, 
        :parms => identity => [:p1, :p2, :p3]
    )

    groups = groupby(df, :pattern)
    combine(groups, 
        :p1 => minimum, :p1 => maximum, 
        :p2 => minimum, :p2 => maximum,
        :p3 => minimum, :p3 => maximum
    ) |> sort

    @test length(groups) == n_part^n_dims
end


@safetestset "5D Hypercube" begin
    using ParameterSpacePartitions
    using ParameterSpacePartitions.TestModels
    using Test, Random, Distributions
    using DataFrames

    Random.seed!(103)
    # partitions per dimension
    n_part = 4
    # dimensions of hypbercue
    n_dims = 5
    # number of starting points
    n_start = 1

    # bounds of the hypercube
    bounds = fill((0, 1), n_dims)
    # partition boundaries
    p_bounds = range(0, 1, n_part + 1)
    hypercube = HyperCube(p_bounds)
 
    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

    init_parms = map(_ -> sample(bounds), 1:n_start)

    options = Options(;
        radius = .15,
        bounds,
        n_iters = 100,
        parallel = false,
        init_parms
    )

    results = find_partitions(
        model, 
        p_fun, 
        options,
        hypercube
    )

    df = DataFrame(results)
    transform!(
        df, 
        :parms => identity => [:p1,:p2,:p3,:p4,:p5]
    )

    groups = groupby(df, :pattern)
    combine(groups, 
        :p1 => minimum, :p1 => maximum, 
        :p2 => minimum, :p2 => maximum,
        :p3 => minimum, :p3 => maximum,
        :p4 => minimum, :p4 => maximum,
        :p5 => minimum, :p5 => maximum
    ) |> sort

    @test length(groups) == n_part^n_dims
end

@safetestset "5D Polytope" begin
    using ParameterSpacePartitions
    using ParameterSpacePartitions.TestModels
    using Test, Random, Distributions
    using DataFrames, LinearAlgebra
    Random.seed!(778)

    # dimensions of the hypbercue
    n_dims = 5
    # number of partitions
    n_part = 50
    # number of starting points
    n_start = 1
    

    # partition boundaries
    polytopes = [Polytope(rand(n_dims)) for i in 1:n_part]
    bounds = fill((0, 1), n_dims)

    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

    init_parms = map(_ -> sample(bounds), 1:n_start)

    options = Options(;
        radius = .1,
        bounds,
        n_iters = 500,
        parallel = false,
        init_parms
    )

    results = find_partitions(
        model, 
        p_fun, 
        options,
        polytopes
    )

    df = DataFrame(results)
    transform!(
        df, 
        :parms => identity => [:p1, :p2, :p3, :p4, :p5]
    )

    groups = groupby(df, :pattern)
    @test length(groups) == n_part
end

@safetestset "HyperSpheres" begin
    using Test, ParameterSpacePartitions
    using ParameterSpacePartitions.TestModels
    using LinearAlgebra, Random, DataFrames, Distributions
    Random.seed!(383)
    
    # dimensions of the hypbercue
    n_dims = 7
    # number of hyperspheres
    n_obj = 50
    # distribution of radii
    r_dist = () -> rand(Uniform(.2, .25))
    # number of starting points
    n_start = 10
    
    # partition boundaries
    hyperspheres = set_locations(r_dist, n_dims, n_obj)
    bounds = fill((0, 1), n_dims)
    
    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)
    
    init_parms = map(_ -> sample(bounds), 1:n_start)
    
    options = Options(;
        radius = .15,
        bounds,
        n_iters = 5000,
        init_parms,
        λ = .05
    )
    
    results = find_partitions(
        model, 
        p_fun, 
        options,
        hyperspheres
    )
    
    df = DataFrame(results)
    groups = groupby(df, :pattern)

    @test length(groups) == (n_obj + 1)
end

@safetestset "update_position!" begin
    using ParameterSpacePartitions
    using Test
    import ParameterSpacePartitions: update_position!, Chain

    parms = [.4]
    proposal = [.3]
    pattern = [4,3]
    chain = Chain(parms, pattern, .3)
    update_position!(chain, proposal, pattern)
    @test chain.parms == [.3]
    @test chain.pattern == [4,3]

    parms = [.4]
    proposal = [.3]
    pattern = [4,3]
    chain = Chain(parms, pattern, .3)
    update_position!(chain, proposal, [4,5])
    @test chain.parms == [.4]
    @test chain.pattern == [4,3]
end

@safetestset "Adapt Radius 5D HyperCube" begin
    using ParameterSpacePartitions
    using ParameterSpacePartitions.TestModels
    using Test, Random, Distributions
    using DataFrames

    Random.seed!(587)
    # partitions per dimension
    n_part = 4
    # dimensions of hypbercue
    n_dims = 5
    # number of starting points
    n_start = 1

    # bounds of the hypercube
    bounds = fill((0, 1), n_dims)
    # partition boundaries
    p_bounds = range(0, 1, n_part + 1)
    hypercube = HyperCube(p_bounds)
 
    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

    init_parms = map(_ -> sample(bounds), 1:n_start)

    t_rate = .4

    options = Options(;
        radius = .40,
        bounds,
        n_iters = 500,
        parallel = false,
        init_parms,
        t_rate,
        λ = .0
    )

    results = find_partitions(
        model, 
        p_fun, 
        options,
        hypercube
    )

    df = DataFrame(results)
    transform!(
        df, 
        :parms => identity => [:p1,:p2,:p3,:p4,:p5]
    )

    groups = groupby(df, :chain_id)
    mean_accept = combine(groups, :acceptance=>mean=>:mean)
    # show that initial radius results in poor acceptance
    @test mean(mean_accept.mean) < .10
    
    options = Options(;
        radius = .40,
        bounds,
        n_iters = 500,
        parallel = false,
        init_parms,
        t_rate,
        λ = .2
    )

    results = find_partitions(
        model, 
        p_fun, 
        options,
        hypercube
    )

    df = DataFrame(results)
    transform!(
        df, 
        :parms => identity => [:p1,:p2,:p3,:p4,:p5]
    )

    groups = groupby(df, :chain_id)
    mean_accept = combine(groups, :acceptance=>mean=>:mean)
    @test mean(mean_accept.mean) ≈ t_rate atol = .03
    @test std(mean_accept.mean) < .05
end

@safetestset "Adapt Radius 5D Polytope" begin
    using ParameterSpacePartitions
    using ParameterSpacePartitions.TestModels
    using Test, Random, Distributions
    using DataFrames, LinearAlgebra
    Random.seed!(2002)

    # dimensions of the hypbercue
    n_dims = 5
    # number of partitions
    n_part = 50
    # number of starting points
    n_start = 1
    

    # partition boundaries
    polytopes = [Polytope(rand(n_dims)) for i in 1:n_part]
    bounds = fill((0, 1), n_dims)

    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

    init_parms = map(_ -> sample(bounds), 1:n_start)

    options = Options(;
        radius = 1.6,
        bounds,
        n_iters = 500,
        parallel = false,
        init_parms,
        adapt_radius! = no_adaption!
    )

    results = find_partitions(
        model, 
        p_fun, 
        options,
        polytopes
    )

    df = DataFrame(results)
    transform!(
        df, 
        :parms => identity => [:p1, :p2, :p3, :p4, :p5]
    )

    groups = groupby(df, :pattern)
    @test length(groups) == n_part

    groups = groupby(df, :chain_id)
    mean_accept = combine(groups, :acceptance=>mean=>:mean)
    # show that initial radius results in poor acceptance
    @test mean(mean_accept.mean) < .10

    t_rate = .4 

    options = Options(;
        radius = 1.6,
        bounds,
        n_iters = 1000,
        parallel = false,
        init_parms,
        λ = 0.2,
        t_rate
    )

    results = find_partitions(
        model, 
        p_fun, 
        options,
        polytopes
    )

    df = DataFrame(results)
    transform!(
        df, 
        :parms => identity => [:p1, :p2, :p3, :p4, :p5]
    )

    groups = groupby(df, :pattern)
    @test length(groups) == n_part

    groups = groupby(df, :chain_id)
    mean_accept = combine(groups, :acceptance=>mean=>:mean)
    # show that initial radius results in poor acceptance
    @test mean(mean_accept.mean) ≈ t_rate atol = .04
    @test std(mean_accept.mean) < .05
end

@safetestset "Process New Patterns" begin
    using ParameterSpacePartitions, Test
    import ParameterSpacePartitions: process_new_patterns!, Chain, is_new

    all_patterns = [[1,2],[1,3]]
    patterns = [[1,2],[1,4]]

    chains = [Chain([.3,.3], [1,2], .2),Chain([.3,.3], [1,3], .2)]
    chains[1].acceptance[1] = false
    chains[2].acceptance[1] = false

    parms = [[.3,.4],[.3,.5]]
    options = Options(;
        radius = .2, 
        bounds = [(0,1),(0,1)], 
        n_iters = 100, 
        init_parms = [[.2,.3]]
    )

    @test !is_new(all_patterns, patterns[1])
    @test is_new(all_patterns, patterns[2])

    process_new_patterns!(all_patterns, patterns, parms, chains, options)

    @test length(chains) == 3
    @test chains[3].pattern == [1,4]
end

@safetestset "update_position!" begin
    using ParameterSpacePartitions, Test
    import ParameterSpacePartitions: update_position!, Chain

    chain = Chain([.3,.3], [1,2], .2)
    proposal = [.2,.4]
    pattern = [1,2]
    update_position!(chain, proposal, pattern)

    @test chain.parms == proposal
    @test length(chain.acceptance) == 2
    @test chain.acceptance[2] == true

    chain = Chain([.3,.3], [1,2], .2)
    proposal = [.2,.4]
    pattern = [1,3]
    update_position!(chain, proposal, pattern)

    @test chain.parms == [.3,.3]
    @test length(chain.acceptance) == 2
    @test chain.acceptance[2] == false
end

@safetestset "adjust_parms!" begin
    using ParameterSpacePartitions, Test
    import ParameterSpacePartitions: adjust_parms!, Chain

    bounds = [(0,1),(0,1)]
    parms = [.5,.5]
    adjust_parms!(parms, bounds)
    @test parms == [.5,.5]

    bounds = [(0,1),(0,1)]
    parms = [.5,-1]
    adjust_parms!(parms, bounds)
    @test parms == [.5,.0]

    bounds = [(0,1),(0,1)]
    parms = [1.2,.5]
    adjust_parms!(parms, bounds)
    @test parms == [1.0,.50] 
end

@safetestset "3D volume" begin
    using Test, QHull, Random, Distributions, LinearAlgebra
    using ParameterSpacePartitions: volume_ellipsoid, sample_ellipsoid_surface
    Random.seed!(584)
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

@safetestset "4D volume" begin
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