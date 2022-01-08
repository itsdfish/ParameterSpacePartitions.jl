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

    # partition boundaries
    p_bounds = range(0, 1, n_part + 1)
    bounds = fill((0, 1), n_dims)

    model(parms) = parms 

    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

    init_parms = map(_ -> sample(bounds), 1:1)

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

    Random.seed!(99)
    # partitions per dimension
    n_part = 4
    # dimensions of hypbercue
    n_dims = 5

    bounds = fill((0, 1), n_dims)
    # partition boundaries
    p_bounds = range(0, 1, n_part + 1)

    model(parms) = parms 
 
    sample(bounds) = map(b -> rand(Uniform(b...)), bounds)

    init_parms = map(_ -> sample(bounds), 1:1)

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