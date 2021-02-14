"""
    abstract type AbstractFunction <: Function end

An abstract function type.
"""
abstract type AbstractFunction <: Function end

"""
    length(f::AbstractFunction)

Returns the number of outputs of `f`.
"""
Base.length(f::AbstractFunction) = getdim(f)

"""
```
struct FunctionWrapper <: AbstractFunction
    f::Function
    dim::Int
end
```

A function wrapper type that wraps the function `f` where the function `f` is declared to return an output of dimension `dim`.
"""
@params struct FunctionWrapper <: AbstractFunction
    f::Function
    dim::Int
end

"""
    (f::FunctionWrapper)(args...; kwargs...)

Calls the wrapped function in `f` with arguments `args` and keyword arguments `kwargs` and returns the output.
"""
function (f::FunctionWrapper)(args...; kwargs...)
    out = f.f(args...; kwargs...)
    @assert length(out) == getdim(f)
    return out
end

"""
    getfunction(f::FunctionWrapper)

Returns the function wrapped in `f`.
"""
getfunction(f::FunctionWrapper) = f.f

getdim(f::FunctionWrapper) = f.dim

"""
```
struct VectorOfFunctions <: AbstractFunction
    fs
end
```

A struct for a collection of instances of `AbstractFunction`. The dimension of this function is the sum of the dimensions of the individual functions and the outputs are concatenated in a vector.
"""
@params struct VectorOfFunctions <: AbstractFunction
    fs::Union{AbstractVector{<:AbstractFunction}, Tuple{Vararg{AbstractFunction}}}
end

"""
    VectorOfFunctions(fs)

Constructs an instance of `VectorOfFunctions` made of functions `fs`. `fs` can be a vector or tuple of instances of `AbstractFunction`. If a function in `fs` is a `VectorOfFunctions`, it will be unwrapped.
"""
function VectorOfFunctions(fs) end

function VectorOfFunctions(
    fs::Union{
        AbstractVector{>:VectorOfFunctions},
        Tuple{Vararg{>:VectorOfFunctions}},
    },
)
    @assert length(fs) > 0
    new_fs = mapreduce(vcat, fs) do f
        if f isa VectorOfFunctions
            f.fs
        else
            f
        end
    end
    return VectorOfFunctions(new_fs)
end

"""
    (f::VectorOfFunctions)(args...; kwargs...)

Calls the wrapped functions in `f` with arguments `args` and keyword arguments `kwargs` and returns the concatenated output.
"""
function (f::VectorOfFunctions)(args...; kwargs...)
    ys = map(f.fs) do f
        f(args...; kwargs...)
    end
    return vcat(ys...)
end

function getdim(c::VectorOfFunctions)
    if length(c.fs) == 0
        return 0
    else
        return sum(getdim, c.fs)
    end
end

"""
    getfunction(f::VectorOfFunctions, i::Integer)

Returns the `i`^th wrapped function in `f`.
"""
getfunction(f::VectorOfFunctions, i::Integer) = getfunctions(f)[i]

"""
    getfunctions(f::VectorOfFunctions)

Returns the vector or tuple of functions wrapped in `f`.
"""
getfunctions(f::VectorOfFunctions) = f.fs

"""
    abstract type AbstractConstraint <: AbstractFunction end

An abstract constraint type.
"""
abstract type AbstractConstraint <: AbstractFunction end

"""
```
struct Objective <: AbstractFunction
    f
    multiple
end
```

The objective function to be optimized. The objective is always assumed to be of dimension 1. `multiple` is of type `Ref{<:Number}`. When an instance of `Objective` is called, its output will be scaled by a multiplier `multiple[]`. This is 1 by default.
"""
@params struct Objective <: AbstractFunction
    f
    multiple::Base.RefValue
end

"""
    Objective(f, multiple::Number = 1.0)

Constructs an instance of `Objective` wrapping `f`. When an instance of `Objective`, it calls `f` returning the output multiplied by `multiple`. `f` must return a number.
"""
Objective(f, multiple::Number = 1.0) = Objective(f, Ref(multiple))

"""
    (obj::Objective)(args...; kwargs...)

Calls the wrapped function in `obj` with arguments `args` and keyword arguments `kwargs` returning the output multiplied by `obj.multiple[]`. The output of the wrapped function must be a number.
"""
function (o::Objective)(args...; kwargs...)
    out = o.f(args...; kwargs...) * o.multiple[]
    @assert out isa Number
    return out
end

getdim(c::Objective) = 1

"""
    getfunction(obj::Objective)

Returns the function wrapped in `obj`.
"""
getfunction(obj::Objective) = obj.f

"""
```
struct IneqConstraint <: AbstractConstraint
    f::Function
    rhs
    dim::Int
end
```

A struct for an inequality constraint of the form `f(x) .- rhs .<= 0`. The dimension of the constraint is `dim`. Calling the struct will return `f(x) .- rhs`.
"""
@params struct IneqConstraint <: AbstractConstraint
    f::Function
    rhs
    dim::Int
end

"""
    IneqConstraint(f, rhs)

Constructs an instance of `IneqConstraint` with a left-hand-side function `f` and right-hand-side bound `rhs`. If `f` is an instance of `AbstractFunction`, the dimension of the constraint will be `getdim(f)`, otherwise it is assumed to be `length(rhs)` by default.
"""
IneqConstraint(f, rhs) = IneqConstraint(f, rhs, length(rhs))
IneqConstraint(f::AbstractFunction, rhs) = IneqConstraint(f, rhs, getdim(f))

"""
    (c::IneqConstraint)(args...; kwargs...)

Calls the wrapped function in the constraint `c`, `c.f`, with arguments `args` and keyword arguments `kwargs` returning the constraint violation `c.f(args...; kwargs...) .- c.rhs`.
"""
function (c::IneqConstraint)(args...; kwargs...)
    out = c.f(args...; kwargs...) .- c.rhs
    @assert length(out) == getdim(c)
    return out
end

getdim(c::IneqConstraint) = c.dim

"""
    getfunction(f::IneqConstraint)

Returns the function wrapped in `f`.
"""
getfunction(f::IneqConstraint) = f.f

"""
```
struct EqConstraint <: AbstractConstraint
    f::Function
    rhs
    dim::Int
end
```

A struct for an equality constraint of the form `f(x) .- rhs .== 0`. The dimension of the constraint is `dim`. Calling the struct will return `f(x) .- rhs`.
"""
@params struct EqConstraint <: AbstractConstraint
    f::Function
    rhs
    dim::Int
end

"""
    EqConstraint(f, rhs)

Constructs an instance of `EqConstraint` with a left-hand-side function `f` and right-hand-side bound `rhs`. If `f` is an instance of `AbstractFunction`, the dimension of the constraint will be `getdim(f)`, otherwise it is assumed to be `length(rhs)` by default.
"""
EqConstraint(f, rhs) = EqConstraint(f, rhs, length(rhs))
EqConstraint(f::AbstractFunction, rhs) = EqConstraint(f, rhs, getdim(f))

"""
    (c::EqConstraint)(args...; kwargs...)

Calls the wrapped function in the constraint `c`, `c.f`, with arguments `args` and keyword arguments `kwargs` returning the constraint violation `c.f(args...; kwargs...) .- c.rhs`.
"""
function (c::EqConstraint)(args...; kwargs...)
    out = c.f(args...; kwargs...) .- c.rhs
    @assert length(out) == getdim(c)
    return out
end

getdim(c::EqConstraint) = c.dim

"""
    getfunction(f::EqConstraint)

Returns the function wrapped in `f`.
"""
getfunction(f::EqConstraint) = f.f
