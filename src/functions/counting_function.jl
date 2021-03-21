struct CountingFunction{F} <: AbstractFunction
    counter::Base.RefValue{Int}
    f::F
end
getdim(f::CountingFunction) = getdim(f.f)
CountingFunction(f::Function) = CountingFunction(Ref(0), f)
function (f::CountingFunction)(x)
    f.counter[] += 1
    return f.f(x)
end
