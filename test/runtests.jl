using SafeTestsets

@safetestset "Cube" begin
    using ParameterSpacePartitions
    using Test, Random, Distributions
    using DataFrames
    include("functions.jl")

    Random.seed!(541)
    # dimensions of the hypbercue
    n_dims = 3
    # partitions per dimension
    n_part = 4
    # number of starting points
    n_start = 1

    # partition boundaries
    p_bounds = range(0, 1, n_part + 1)
    bounds = fill((0, 1), n_dims)

    model(parms) = parms 

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
        x -> p_fun(x, p_bounds), 
        options
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
    using Test, Random, Distributions
    using DataFrames
    include("functions.jl")

    Random.seed!(103)
    # partitions per dimension
    n_part = 4
    # dimensions of hypbercue
    n_dims = 5
    # number of starting points
    n_start = 1

    bounds = fill((0, 1), n_dims)
    # partition boundaries
    p_bounds = range(0, 1, n_part + 1)

    model(parms) = parms 
 
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
        x -> p_fun(x, p_bounds), 
        options
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
    using Test, Random, Distributions
    using DataFrames, LinearAlgebra
    Random.seed!(778)

    function p_fun(data, points)
        distances = map(p -> norm(data .- p), points)
        _,idx = findmin(distances)
        return idx
    end

    # dimensions of the hypbercue
    n_dims = 5
    # number of partitions
    n_part = 50
    # number of starting points
    n_start = 1
    #

    # partition boundaries
    points = [rand(n_dims) for i in 1:n_part]
    bounds = fill((0, 1), n_dims)

    model(parms) = parms 

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
        x -> p_fun(x, points), 
        options
    )

    df = DataFrame(results)
    transform!(
        df, 
        :parms => identity => [:p1, :p2, :p3, :p4, :p5]
    )

    groups = groupby(df, :pattern)
    @test length(groups) == n_part
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
