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
    Random.seed!(2001)

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
        radius = 1.8,
        bounds,
        n_iters = 500,
        parallel = false,
        init_parms,
        λ = 0.0
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
        radius = 1.8,
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
    @test mean(mean_accept.mean) ≈ t_rate atol = .03
    @test std(mean_accept.mean) < .05
end