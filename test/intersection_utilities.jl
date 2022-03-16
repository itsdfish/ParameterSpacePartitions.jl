function get_axes(Σ, c=2)
    eig_val,_ = eigen(Σ)
    return 2 * c * sqrt.(eig_val)
end

function rand_cov_mat(n_dims)
    x = randn(n_dims, n_dims)
    return x' * x
end