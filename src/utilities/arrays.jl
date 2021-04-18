"""
Judge if an array forms a matrix of m*n numberf
I.e. Judge if the array consists of a series of m points from n-dimensional Euclid space
"""
is_matrix(arr::AbstractArray) = length(arr) == 0 || 
                                            ((arr isa Union{Vector{T}, Vector{Vector{T}}} where {T<:Number}) 
                                                &&
                                            all_size_equals(arr, size(arr[1])))


"""
Alias of is_matrix_of_number
"""
const is_array_of_points = is_matrix

"""
Given an array "arr" and a size "s"(number or tuple), judge if all elements of the array equals to this size
"""
function all_size_equals(arr::AbstractArray, s::T) where {T<:Union{Signed, Tuple}}
    all(e -> size(e) == s, arr)
end
