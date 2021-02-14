"""
    getdim(f::AbstractFunction)

Returns the dimension of the function `f`, i.e. the number of outputs.
"""
function getdim(::AbstractFunction) end

"""
    getfunction(f::AbstractFunction)

Returns any wrapped function inside the function `f`.
"""
function getfunction end
