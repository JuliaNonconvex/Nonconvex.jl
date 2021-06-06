# Adapted from ParameterHandling.jl with the following license.
#=
Copyright (c) 2020 Invenia Technical Computing Corporation

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

"""
    flatten(x)

Returns a "flattened" representation of `x` as a vector of real numbers, and a function
`unflatten` that takes a vector of reals of the same length and returns an object of the
same type as `x`.

`unflatten` is the inverse of `flatten`, so
```julia
julia> x = (randn(5), 5.0, (a=5.0, b=randn(2, 3)));

julia> v, unflatten = flatten(x);

julia> x == unflatten(v)
true
```
"""
function flatten end

maybeflatten(x::Real) = x
maybeflatten(x) = flatten(x)

function flatten(x::Real)
    v = [x]
    unflatten_to_Real(v) = only(v)
    return v, unflatten_to_Real
end

flatten(x::Vector{<:Real}) = (x, identity)

function flatten(x::AbstractVector)
    x_vecs_and_backs = map(val -> flatten(val), x)
    x_vecs, backs = first.(x_vecs_and_backs), last.(x_vecs_and_backs)
    function Vector_from_vec(x_vec)
        sz = _cumsum(map(_length, x_vecs))
        x_Vec = [backs[n](x_vec[sz[n] - _length(x_vecs[n]) + 1:sz[n]]) for n in eachindex(x)]
        return oftype(x, x_Vec)
    end
    return reduce(vcat, x_vecs), Vector_from_vec
end

function flatten(x::AbstractArray)
    x_vec, from_vec = flatten(vec(x))
    Array_from_vec(x_vec) = oftype(x, reshape(from_vec(x_vec), size(x)))
    return x_vec, Array_from_vec
end

function flatten(x::SparseMatrixCSC)
    x_vec, from_vec = flatten(x.nzval)
    Array_from_vec(x_vec) = SparseMatrixCSC(x.m, x.n, x.colptr, x.rowval, from_vec(x_vec))
    return x_vec, Array_from_vec
end

_length(x) = length(x)
_length(::Nothing) = 1

function flatten(x::Tuple)
    x_vecs_and_backs = map(val -> flatten(val), x)
    x_vecs, x_backs = first.(x_vecs_and_backs), last.(x_vecs_and_backs)
    lengths = map(_length, x_vecs)
    sz = _cumsum(lengths)
    function unflatten_to_Tuple(v)
        map(x_backs, lengths, sz) do x_back, l, s
            return x_back(v[s - l + 1:s])
        end
    end
    return reduce(vcat, x_vecs), unflatten_to_Tuple
end

function flatten(x::NamedTuple)
    x_vec, unflatten = flatten(values(x))
    function unflatten_to_NamedTuple(v)
        v_vec_vec = unflatten(v)
        return typeof(x)(v_vec_vec)
    end
    return x_vec, unflatten_to_NamedTuple
end

function flatten(d::AbstractDict, ks = collect(keys(d)))
    _d = OrderedDict(k => d[k] for k in ks)
    d_vec, unflatten = flatten(collect(values(_d)))
    function unflatten_to_Dict(v)
        v_vec_vec = unflatten(v)
        return OrderedDict(key => v_vec_vec[n] for (n, key) in enumerate(keys(_d)))
    end
    return d_vec, unflatten_to_Dict
end
function ChainRulesCore.rrule(::typeof(flatten), d::AbstractDict, ks)
    _d = OrderedDict(k => d[k] for k in ks)
    d_vec, un = flatten(_d, ks)
    return (d_vec, un), Δ -> begin
        (NO_FIELDS, un(Δ[1]), nothing)
    end
end

struct Unflatten{F} <: Function
    unflatten::F
end
(f::Unflatten)(x) = f.unflatten(x)

function _merge(d1, d2::AbstractDict)
    @assert eltype(values(d1)) <: Real
    _d = OrderedDict(k => zero(v) for (k, v) in d1)
    return sort!(merge(_d, OrderedDict(d2)))
end
_merge(::Any, d2) = d2

function ChainRulesCore.rrule(un::Unflatten, v)
    x = un(v)
    return x, Δ -> begin
        _Δ = _merge(x, Δ)
        return (NO_FIELDS, flatten(_Δ)[1])
    end
end

function flatten(::Nothing)
    return [0.0], _ -> nothing
end
function flatten(::ZeroTangent)
    return [0.0], _ -> ZeroTangent()
end
function flatten(::Tuple{})
    return Float64[], _ -> ()
end
function flatten(x)
    v, un = flatten(ntfromstruct(x))
    return identity.(v), Unflatten(y -> structfromnt(typeof(x), un(y)))
end

macro constructor(T)
    return flatten_expr(T, T)
end
macro constructor(T, C)
    return flatten_expr(T, C)
end
flatten_expr(T, C) = quote
    function flatten(x::$(esc(T)))
        v, un = flatten(ntfromstruct(x))
        return identity.(v), Unflatten(y -> structfromnt($(esc(C)), un(y)))
    end
end

_cumsum(x) = cumsum(x)
if VERSION < v"1.5"
    _cumsum(x::Tuple) = (_cumsum(collect(x))..., )
end
