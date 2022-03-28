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

    if !isposdef(Q2b)
        println("cov1")
        cov1 ./= c^2
        println(cov1)
        println("cov2")
        cov2 ./= c^2
        println(cov2)
    end

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

"""
    intersects(chain1, chain2, c=2)

Tests whether two hyperellipsoids intersect. The test returns true if 
the hyperellipsoids intersection and false otherwise. 

# Arguments

- `chain1`: chain object 
- `chain2`: chain object
- `c=2`: ellipse sclar
"""
function intersects(chain1, chain2, c=2)
    mat1 = to_matrix(chain1)
    mat2 = to_matrix(chain2)
    μ1 = mean(mat1, dims=1)[:]
    μ2 = mean(mat2, dims=1)[:]
    cov1 = cov(mat1)
    add_variance!(cov1)
    cov2 = cov(mat2)
    add_variance!(cov2)
    if !isposdef(cov1)
        println("cov1")
        println(cov1)
    end
    if !isposdef(cov2)
        println("cov2")
        println(cov2)
    end
    return intersects(μ1, μ2, cov1, cov2)
end

function add_variance!(x)
    if any(x -> isapprox(x, 0; atol = 1e-10), diag(x))
        x[diagind(x)] .= eps()
    end
    return nothing 
end

to_matrix(x) = mapreduce(permutedims, vcat, x.all_parms)

"""
    get_group_indices(chains, chain_indices)

Sorts chains of the same pattern into non-overlapping groups. The vector 
[[1,2],[3,4]] indices chains 1 and 2 are located in region R₁ and chains 
3 and 4 are located in region R₂.

# Arguments

- `chains`: a vector of chains 
- `group`: a vector of indices corresponding to chains with the same pattern
"""
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
                merge_chains!(chains[indices[g][1]], chains[c])
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

"""
    remove_redundant_chains!(chains, indices)

Removes chains that have the same pattern and location in the parameter space. 
For example, in the vector `indices` [[1,3],[34,50]], chains indexed in positions 1 and 3 have the same 
pattern and location, as do chains indexed at positions 34 and 50. Only the first element of each sub-vector
is retained (i.e., [1,50])

# Arguments

- `chains`: a vector containing all chain objects 
- `indices`: a nested vector that maps chains of the same pattern and location to `chains`
"""
function remove_redundant_chains!(chains, indices)
    k_indices = Int[]
    for i in indices
        push!(k_indices, i[1])
    end
    r_indices = setdiff(vcat(indices...), k_indices)
    sort!(r_indices)
    deleteat!(chains, r_indices)
    return nothing
end

"""
    group_by_pattern(chains)

Groups chains according to pattern and returns a nested vector of chain indices. The vector 
[[1,2],[3,4]] indices chains 1 and 2 are located in region R₁ and chains 
3 and 4 are located in region R₂. 

# Arguments

- `chains`: a vector of all chains 

"""
function group_by_pattern(chains)
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

function remove_nonposdef!(chains)
    n_chains = length(chains)
    test_vals = fill(false, n_chains)
    for c in 1:n_chains
        mat = to_matrix(chains[c])
        covar = cov(mat)
        test_vals[c] = !isposdef(covar)
    end
    deleteat!(chains, test_vals)
    return nothing
end

"""
    make_unique!(chains, options)

This function sorts chains by pattern and merges chains that are in the same region 

# Arguments

- `chains`: a vector of all chains 
- `options`: an Options configuration object for the PSP algorithm
"""
function make_unique!(chains, options)
    remove_nonposdef!(chains)
    chain_indices = group_by_pattern(chains)
    all_indices = Vector{Vector{Int}}()
    for c in chain_indices
        temp = get_group_indices(chains, c)
        push!(all_indices, temp...)
    end
    # merge_chains!(chains, all_indices, options)
    remove_redundant_chains!(chains, all_indices)
    return nothing
end

function merge_chains!(chains, indices, options) 
    return merge_chains!(chains, indices, options.max_merge)
end

"""
    merge_chains!(chains, indices, max_merge::Int)  

Merges a set of chains if they have the same pattern and their regions in the parameter space intersect.

# Arguments

- `chains`: a vector of chains with the same pattern 
- `indices`: a vector of indices that map the position of each chain to a vector of all chains 
"""
function merge_chains!(chains, indices, max_merge::Int)
    max_merge == 0 ? (return nothing) : nothing 
    for idx in indices 
        length(idx) == 1 ? (continue) : nothing
        n = min(length(idx), max_merge + 1)
        for c in 2:n 
            merge_chains!(chains[idx[1]], chains[idx[c]])
        end
    end
    return nothing
end

"""
    merge_chains!(chain1, chain2)

Merges chain2 into chain1 on fields `all_parms`, `acceptance`, and `radii`.

# Arguments

- `chain1`: a chain object 
- `chain2`: a chain object 
"""
function merge_chains!(chain1, chain2)
    push!(chain1.all_parms, chain2.all_parms...)
    push!(chain1.acceptance, chain2.acceptance...)
    push!(chain1.radii, chain2.radii...)
    return nothing 
end