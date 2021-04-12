abstract type AbstractAggregation <: AbstractFunction end
getdim(::AbstractAggregation) = 1

@params struct NonNegSumOfSquares <: AbstractAggregation
    c::Real
end
(f::NonNegSumOfSquares)(x::AbstractVector) = f.c * sum(max.(0, x).^2)

@params struct NonNegSumOfRoots <: AbstractAggregation
    c::Real
end
(f::NonNegSumOfRoots)(x::AbstractVector) = f.c * sum(sqrt.(max.(0, x)))

@params struct NonNegSum <: AbstractAggregation
    c::Real
end
(f::NonNegSum)(x::AbstractVector) = f.c * sum(max.(0, x))

@params struct WeightedSum <: AbstractAggregation
    weights::AbstractVector
end
(f::WeightedSum)(x::AbstractVector) = dot(f.weights, x)
