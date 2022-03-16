# https://www.youtube.com/watch?v=OPSCKXXvWiM&ab_channel=TheOrganicChemistryTutor
# https://math.stackexchange.com/questions/1447730/drawing-ellipse-from-eigenvalue-eigenvector
# https://stats.stackexchange.com/questions/164741/how-to-find-the-maximum-axis-of-ellipsoid-given-the-covariance-matrix
"""
    intersects(μ1, μ2, cov1, cov2)

Tests whether two hyperellipsoids intersect. The test returns true if 
the hyperellipsoids intersection and false otherwise. 

# Arguments

- `μ1`: centroid of ellipsoid 1 
- `μ2`: centroid of ellipsoid 2
- `cov1`: covariance matrix of ellipsoid 1
- `cov2`: covariance matrix of ellipsoid 2
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