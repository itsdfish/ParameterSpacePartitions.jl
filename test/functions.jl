function p_fun(data, p_bounds)
    nb = length(p_bounds)
    nd = length(data)
    vals = fill(-100, nd)
    for j in 1:nd
        for i in 1:(nb-1) 
            if (data[j] ≥ p_bounds[i]) && (data[j] ≤ p_bounds[i+1])
                vals[j] = i 
                continue
            end
        end
    end
    return vals
end