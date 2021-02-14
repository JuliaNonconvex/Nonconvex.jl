"""
    @params struct_def

A macro that changes all fields' types to type parameters while respecting the type bounds specificied by the user. For example:
```
@params struct MyType{T}
    f1::T
    f2::AbstractVector{T}
    f3::AbstractVector{<:Real}
    f4
end
```
will define:
```
struct MyType{T, T1 <: T, T2 <: AbstractVector{T}, T3 <: AbstractVector{<:Real}, T4}
    f1::T1
    f2::T2
    f3::T3
    f4::T4
end
```
The default non-parameteric constructor, e.g. `MyType(f1, f2, f3, f4)`, will always work if all the type parameters used in the type header are used in the field types. Using a type parameter in the type header that is not used in the field types such as:
```
struct MyType{T}
    f1
    f2
end
```
is not recommended.
"""
macro params(struct_expr)
    header = struct_expr.args[2]
    fields = @view struct_expr.args[3].args[2:2:end]
    params = []
    for i in 1:length(fields)
        x = fields[i]
        T = gensym()
        if x isa Symbol
            push!(params, T)
            fields[i] = :($x::$T)
        elseif x.head == :(::)
            abstr = x.args[2]
            var = x.args[1]
            push!(params, :($T <: $abstr))
            fields[i] = :($var::$T)
        end
    end
    if header isa Symbol && length(params) > 0
        struct_expr.args[2] = :($header{$(params...)})
    elseif header.head == :curly
        append!(struct_expr.args[2].args, params)
    elseif header.head == :<:
        if struct_expr.args[2].args[1] isa Symbol
            name = struct_expr.args[2].args[1]
            struct_expr.args[2].args[1] = :($name{$(params...)})
        elseif header.head == :<: && struct_expr.args[2].args[1] isa Expr
            append!(struct_expr.args[2].args[1].args, params)
        else
            error("Unidentified type definition.")
        end
    else
        error("Unidentified type definition.")
    end
    esc(struct_expr)
end
