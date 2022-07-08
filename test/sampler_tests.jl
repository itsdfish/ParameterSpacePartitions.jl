@safetestset "Process New Patterns" begin
    using ParameterSpacePartitions, Test
    import ParameterSpacePartitions: process_new_patterns!, Chain, is_new

    all_patterns = [[1,2],[1,3]]
    patterns = [[1,2],[1,4]]

    chains = [Chain(1, [.3,.3], [1,2], .2),Chain(2, [.3,.3], [1,3], .2)]
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

@safetestset "in_bounds" begin
    using ParameterSpacePartitions, Test
    import ParameterSpacePartitions: in_bounds, Chain

    bounds = [(0,1),(0,1)]
    parms = [.5,.5]
    in_bounds(parms, bounds)
    @test in_bounds(parms, bounds)

    bounds = [(0,1),(0,1)]
    parms = [.5,-1]
    @test !in_bounds(parms, bounds)

    bounds = [(0,1),(0,1)]
    parms = [1.2,.5]
    @test !in_bounds(parms, bounds) 
end

@safetestset "update_position!" begin
    using ParameterSpacePartitions
    using Test
    import ParameterSpacePartitions: update_position!, Chain

    parms = [.4]
    proposal = [.3]
    pattern = [4,3]
    bounds = [(0,1)]
    chain = Chain(1, parms, pattern, .3)
    update_position!(chain, proposal, pattern, bounds)
    @test chain.parms == [.3]
    @test chain.pattern == [4,3]
    @test chain.acceptance[2] == true

    parms = [.4]
    proposal = [.3]
    pattern = [4,3]
    chain = Chain(1, parms, pattern, .3)
    update_position!(chain, proposal, [4,5], bounds)
    @test chain.parms == [.4]
    @test chain.pattern == [4,3]
    @test chain.acceptance[2] == false


    parms = [.4]
    proposal = [-.3]
    pattern = [4,3]
    chain = Chain(1, parms, pattern, .3)
    update_position!(chain, proposal, pattern, bounds)
    @test chain.parms == [.4]
    @test chain.pattern == [4,3]
    @test chain.acceptance[2] == false

end

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

    df = find_partitions(
        model, 
        p_fun, 
        options,
        hypercube
    )

    groups = groupby(df, :pattern)
    length(groups)
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

    df = find_partitions(
        model, 
        p_fun, 
        options,
        hypercube
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

    df = find_partitions(
        model, 
        p_fun, 
        options,
        polytopes
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
    n_start = 1
    
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
    )
    
    df = find_partitions(
        model, 
        p_fun, 
        options,
        hyperspheres
    )
    
    groups = groupby(df, :pattern)

    @test length(groups) == (n_obj + 1)
end