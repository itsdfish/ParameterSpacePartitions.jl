"""
    intersects(μ1, μ2, cov1, cov2, c=2)

Tests whether two hyperellipsoids intersect. The test returns true if 
the hyperellipsoids intersection and false otherwise. 

# Arguments

- `μ1`: centroid of ellipsoid 1 
- `μ2`: centroid of ellipsoid 2
- `cov1`: covariance matrix of ellipsoid 1
- `cov2`: covariance matrix of ellipsoid 2
- `c=2`: ellipse sclar
"""
function intersects(μ1, μ2, cov1, cov2, c=2)
    cov1 .*= c^2
    cov2 .*= c^2
    cho1 = cholesky(cov1)
    inv1 = inv(cho1.U)
     
    _Q2b = inv1' * cov2 * inv1
    Q2b = Symmetric(_Q2b)
    if !issymmetric(_Q2b)
        ϵ = sum(abs.(_Q2b .- _Q2b')) / length(_Q2b)
        if ϵ > 1e-10
            @warn "Q2b is not symmetric: ϵ = $(ϵ).   Q2b = $(Q2b)"
        end
    end
    c2b = inv1' * (μ2 - μ1)
    c2c = (c2b' * inv(cholesky(Q2b).U))'
    v2c = -c2c / sqrt(c2c' * c2c)

    test_point = (v2c' * cholesky(Q2b).U)' .+ c2b
    if test_point' * test_point < 1
        return true 
    elseif all(sign.(test_point) ≠ sign.(c2b))
        return true 
    else
        return false
    end
end


function intersects(chain1, chain2, c=2)
    mat1 = to_matrix(chain1)
    mat2 = to_matrix(chain2)
    μ1 = mean(mat1, dims=1)[:]
    μ2 = mean(mat2, dims=1)[:]
    cov1 = cov(mat1)
    cov2 = cov(mat2)
    return intersects(μ1, μ2, cov1, cov2)
end

to_matrix(x) = mapreduce(permutedims, vcat, x.all_parms)

function get_group_indices(chains, chain_indices)
    # group chain indices according to region
    indices = Vector{Vector{Int}}()
    # first index group will have c = 1
    push!(indices, [chain_indices[1]])
    n_groups = length(indices)
    n_chains = length(chain_indices)
    # group index 
    g = 1
    # loop through each chain
    for i ∈ 2:n_chains
        # loop through each index group
        c = chain_indices[i]
        while g ≤ n_groups
            # if chain c matches region index in group g, 
            # push c into index group g
            # (a stand-in for intersection test)
            if intersects(chains[indices[g][1]], chains[c])
                push!(indices[g], c)
                break
            end
            g += 1
        end
        # add new index group
        if g > n_groups
            # add new group
            push!(indices, [c])
            # increment number of groups 
            n_groups += 1
        end
        # reset group counter 
        g = 1
    end
    return indices
end

function remove_redundant_chains!(chains, indices)
    k_indices = Int[]
    for i in indices
        push!(k_indices, i[1])
    end
    r_indices = setdiff(vcat(indices...), k_indices)
    deleteat!(chains, r_indices)
    return nothing
end

function group_by_pattern(chains::T) where {T}
    patterns = map(c -> c.pattern, chains)
    u_patterns = unique(patterns)
    n_patterns = length(u_patterns)
    chain_indices = [Vector{Int}() for _ in 1:n_patterns]
    for c in 1:length(chains) 
        g = findfirst(p -> chains[c].pattern == p, u_patterns)
        push!(chain_indices[g], c)
    end
    return chain_indices
end

function make_unique!(chains)
    chain_indices = group_by_pattern(chains)
    all_indices = Vector{Vector{Int}}()
    for c in chain_indices
        temp = get_group_indices(chains, c)
        push!(all_indices, temp...)
    end
        
    remove_redundant_chains!(chains, all_indices)
    return nothing
end