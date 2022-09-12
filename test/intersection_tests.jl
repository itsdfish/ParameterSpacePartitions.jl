@testset verbose = true "Intersection Test" begin
    @safetestset "Non-intersecting hyperellipsoids" begin
        # Let d be the distance between centroids of ellispoids e1 and e2. The test should always return false
        # if d > m_1 / 2 + m_2 /2, where m_1 is the major axis of e1 and m_2 is the major axis of e2. 
        # returns axes of ellipsoid with scaling factor of 2
        using ParameterSpacePartitions, Test, LinearAlgebra
        using Random
        include("intersection_utilities.jl")
        Random.seed!(8744)
    
        for _ in 1:1
            n_dims = rand(2:10)
            μ₁ = fill(0, n_dims)
            Σ₁ = rand_cov_mat(n_dims)
            Σ₂ = rand_cov_mat(n_dims)
            r₁ = get_axes(Σ₁) / 2
            r₂ = get_axes(Σ₂) / 2
            r_max  = maximum(r₁) + maximum(r₂)
            x = randn(n_dims)
            μ₂ = (x / sqrt(x' * x)) * r_max * 1.01
            @test intersects(μ₁, μ₂, Σ₁, Σ₂) == false
        end 
    end 

    # @safetestset "Intersecting hyperellipsoids" begin
    #     using ParameterSpacePartitions, Test, LinearAlgebra
    #     using Random
    #     include("intersection_utilities.jl")
    #     Random.seed!(354)
    
    #     for _ in 1:100
    #         n_dims = rand(2:10)
    #         μ₁ = fill(0, n_dims)
    #         Σ₁ = rand_cov_mat(n_dims)
    #         Σ₂ = rand_cov_mat(n_dims)
    #         r₁ = get_axes(Σ₁) / 2
    #         r₂ = get_axes(Σ₂) / 2
    #         r_min  = min(minimum(r₁), minimum(r₂))
    #         x = randn(n_dims)
    #         # should always intersect if less than twice minimum axis
    #         μ₂ = (x / sqrt(x' * x)) * 2 * r_min * .99
    #         @test intersects(μ₁, μ₂, Σ₁, Σ₂)
    #     end 
    # end 

    # @safetestset "hyperspheres" begin
    #     using ParameterSpacePartitions, Test, LinearAlgebra
    #     using Random
    #     include("intersection_utilities.jl")
    #     Random.seed!(6521)
    
    #     for _ in 1:100
    #         n_dims = rand(2:10)
    #         μ₁ = fill(0, n_dims)
    #         Σ₁ = Array(I(n_dims) * 1.0)
    #         Σ₂ = Array(I(n_dims) * 1.0)
    #         r = 2 * 2
    #         x = randn(n_dims)
    #         μ₂ = (x / sqrt(x' * x)) * r * .99
    #         @test intersects(μ₁, μ₂, Σ₁, Σ₂)
    #     end 

    #     for _ in 1:100
    #         n_dims = rand(2:10)
    #         μ₁ = fill(0, n_dims)
    #         Σ₁ = Array(I(n_dims) * 1.0)
    #         Σ₂ = Array(I(n_dims) * 1.0)
    #         # radius of hypersphere 
    #         r = 2 
    #         x = randn(n_dims)
    #         # distance is twice the radius plus 1%
    #         μ₂ = (x / sqrt(x' * x)) * 2 * r * 1.01
    #         @test intersects(μ₁, μ₂, Σ₁, Σ₂) == false
    #     end
    # end 

    # @safetestset "remove redundant chains" begin
    #     using ParameterSpacePartitions, Test, LinearAlgebra
    #     using ParameterSpacePartitions: Chain
    #     using ParameterSpacePartitions: intersects, remove_redundant_chains!
    #     using ParameterSpacePartitions: get_group_indices, group_by_pattern
    #     using ParameterSpacePartitions: merge_chains!
    #     using Statistics
        
    #     chains = [
    #         Chain(1, [.1,.3], 1, 1),
    #         Chain(2, [.1,.3,], 1, 1),
    #         Chain(3, [.1,.3,], 1, 1),
    #         Chain(4, [.1,.3,], 2, 1),
    #         Chain(5, [.1,.3,], 2, 1),
        
    #     ]
    #     points = [rand(2) for _ in 1:100]
    #     push!(chains[1].all_parms, points...)
    #     push!(chains[2].all_parms, points...)
    #     points = [rand(2) .+ 2 for _ in 1:100]
    #     push!(chains[3].all_parms, points...)
    #     push!(chains[4].all_parms, points...)
    #     push!(chains[5].all_parms, points...)
        
    #     c_indices = group_by_pattern(chains)
    #     @test c_indices[1] == [1,2,3]
    #     @test c_indices[2] == [4,5]

    #     indices = map(c -> get_group_indices(chains, c), c_indices)
    #     indices = vcat(indices...)
    #     remove_redundant_chains!(chains, indices)

    #     @test length(chains) == 3
    #     @test chains[1].chain_id == 1
    #     @test chains[2].chain_id == 3
    #     @test chains[3].chain_id == 4
    # end 
end

