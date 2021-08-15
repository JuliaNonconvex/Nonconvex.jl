# Funcs
function is_triangular(n)
    if (n < 0)
        return false
    end
    c = (-2 * n)
    a, b = 1, 1
    d = (b * b) - (4 * a * c)
    if (d < 0)
        return false
    end
    root1 = ( -b + sqrt(d)) / (2 * a)
    root2 = ( -b - sqrt(d)) / (2 * a)
 
    if (root1 > 0 && floor(root1) == root1)
        return true
    elseif (root2 > 0 && floor(root2) == root2)
        return true
    else
        return false
    end
end

function length_to_dim(l)
    @assert isinteger(l) "Dimension numbers shoule be a integer. "
    trunc(Int, sqrt(l))
end

function lowertriangind(mat::Matrix; lower::Bool=true, separate::Bool=false)
    indices = [i for i in CartesianIndices(mat) if i[1]>i[2]]
    return LinearIndices(mat)[indices]
end
