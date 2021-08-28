# Funcs
function lowertriangind(mat::Matrix)
    indices = [i for i in CartesianIndices(mat) if i[1]>i[2]]
    return LinearIndices(mat)[indices]
end
