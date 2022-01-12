module TestModels

    using LinearAlgebra

    export p_fun, 
        model,
        set_locations 

    export HyperCube,
        Polytope,
        HyperSphere

    """
        model(parms, args...; kwargs...) 
    
    Model function for test models. It is an identity function. 

    # Arguments

    - `parms`: a proposed location 

    """
    model(parms, args...; kwargs...) = parms

    """
        HyperCube

    A hypercube oject divided into partitions. 

    # Field 

    - `p_bounds`: a vector of boundaries of partitions
    """
    struct HyperCube
        p_bounds::Vector{Float64}
    end

    """
        p_fun(location, hypercube::HyperCube, args...; kwargs...)

    Returns the qualitative pattern associated with a hypercube. 

    # Arguments

    - 'location`: location of chain in the hypercube
    - `p_bounds`: boundaries of partitions for each dimension 

    """
    function p_fun(location, hypercube::HyperCube, args...; kwargs...)
        p_bounds = hypercube.p_bounds
        nb = length(p_bounds)
        nd = length(location)
        vals = fill(-100, nd)
        for j in 1:nd
            for i in 1:(nb-1) 
                if (location[j] ≥ p_bounds[i]) && (location[j] ≤ p_bounds[i+1])
                    vals[j] = i 
                    continue
                end
            end
        end
        return vals
    end

    """
        Polytope
    
    A polytope object containing centroid locations for each partition 

    # Fields

    - `location`: location of centroid for each partition
    """
    struct Polytope
        location::Vector{Float64}
    end


    """
        p_fun(location, points::Vector{Polytope}, args...; kwargs...)

    Returns the qualitative pattern associated with a polytope. 

    # Arguments

    - 'location`: location of chain in the hypercube
    - `points`: a vector of polytopes containing locations

    """
    function p_fun(location, points::Vector{Polytope}, args...; kwargs...)
        distances = map(p -> norm(location .- p.location), points)
        _,idx = findmin(distances)
        return idx
    end

    """
        HyperSphere

    A hypersphere object with arbitrary dimensionality. 

    # Fields 

    - `location`: a vector representing the centroid of the hypersphere
    - `radius`: the radius of the hypersphere
    """
    mutable struct HyperSphere
        location::Vector{Float64}
        radius::Float64
    end
    
    function HyperSphere(;location, radius)
        return HyperSphere(location, radius)
    end
    
    function in_contact(v1::HyperSphere, v2::HyperSphere)
        dist = norm(v1.location - v2.location)
        return dist ≤ max(v1.radius, v2.radius)
    end
    
    function in_contact(v1::HyperSphere, proposal::Vector{<:Number})
        dist = norm(v1.location - proposal)
        return dist ≤ v1.radius
    end
    
    not_valid(v, vs)  = any(x -> in_contact(v, x), vs) || out_of_bounds(v)
    
    out_of_bounds(v) = any((v.location .- v.radius) .< 0) || any((v.location .+ v.radius) .> 1)

    """
        set_locations(r_dist, n_dims, n_obj)
    
    Sets the location of hypspheres within a unit hypercube
    
    # Arguments

    - `r_dist`: a function that samples a value for the radius 
    - `n_dims`: the number of dimensions of the hypersphere 
    - `n_obj`: the number of hyperspheres

    """
    function set_locations(r_dist, n_dims, n_obj)
        hs = [HyperSphere(;location = rand(n_dims), radius = r_dist())]
        for i in 1:(n_obj-1)
            h = HyperSphere(;location = rand(n_dims), radius = r_dist())
            while not_valid(h, hs)
                h.location = rand(n_dims)
            end
            push!(hs, h)
        end
        return hs
    end

    """
        p_fun(location, hyperspheres::Vector{HyperSphere}, args...; kwargs...)

    Returns the qualitative data pattern for hypspheres embedded in a hypercube. -100 is
    returned if the proposed location is not in a hypsphere. Otherwise, an index of the hypersphere 
    in which the proposed location is contained. 

    # Arguments

    - `location`: a proposed location in a unit hypercube
    - `hyperspheres`: a vector of `HyperSphere` objects
    """
    function p_fun(location, hyperspheres::Vector{HyperSphere}, args...; kwargs...)
        for (i,h) in enumerate(hyperspheres)
            in_contact(h, location) ? (return i) : nothing
        end
        # invalid area (not in a hypersphere)
        return -100 
    end
end