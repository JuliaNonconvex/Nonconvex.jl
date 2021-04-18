"""
******************************************************************
*  Hyperparameter optimization using Hyperopt.jl
******************************************************************
"""


"""
    HyperoptResult

When using multiple x0 in [`optimize`](@ref), return this result, including following fields: 
- `optima`: The best result among all searches, which is a [`AbstractResult`](@ref)
- `search_result`: An instance of ['Hyperoptimizer'](@ref) in Hyperopt.jl, which stores all results and initial values of optimization process
- `results`: All search results, including optima and non-optima. Will not be saved if `keepall` is set to false.
"""
struct HyperoptResult <: AbstractResult
    optima::AbstractResult
    optimizer::Hyperoptimizer
    results::Array{AbstractResult}
end
function Base.getproperty(result::HyperoptResult, v::Symbol)
    optima = getfield(result, :optima)
    if hasfield(typeof(optima), v)
        getfield(optima, v)
    else
        getfield(result, v)
    end
end


"""
    HyperSearcher

The abstract struct that stores options to perform hyperparameter search.
"""
abstract type AbstractHyperoptOptions end


"""
    X0OptOptions: Options performing starting point optimization

- `x0_lb`: Lower bound of starting point, if don't specify it, the default value will be `nothing`, 
            then will end up be replaced by the lower bound of optimization problem.
- `x0_rb`: Hier bound of starting point, same as above. 
- `searchspace_size::Integer`: How many potential starting points we generate.
- `iters::Integer`: Among all generated potential starting points, how many of them will be evaluated. 
- `sampler::Hyperopt.Sampler`: An instance of ['Hyperopt.Sampler'](@ref), which decides search algorithm. 
- `verbose::Bool`: Printing information or not.
- `keepall::Bool`: Keep all the search result of selected potential starting points or not.
"""
struct X0OptOptions <: AbstractHyperoptOptions
    x0_lb
    x0_rb
    searchspace_size::Integer
    iters::Integer
    sampler::Hyperopt.Sampler
    verbose::Bool
    keepall::Bool
    function X0OptOptions(x0_lb, x0_rb, searchspace_size, iters, sampler, verbose, keepall)
        @assert x0_lb === nothing || x0_rb === nothing || size(x0_lb) == size(x0_rb) "x0_lb and x0_rb should be same size. "
        return new(x0_lb, x0_rb, searchspace_size, iters, sampler, verbose, keepall)
    end
end
function X0OptOptions(;x0_lb=nothing, x0_rb=nothing, searchspace_size=1000, iters=20, 
                        sampler=RandomSampler(), 
                        verbose=true,
                        keepall=true)
    X0OptOptions(x0_lb, x0_rb, searchspace_size, iters, sampler, verbose, keepall)
end
function X0OptOptions(x0_lb=nothing, x0_rb=nothing; searchspace_size=1000, iters=20, 
                                sampler=RandomSampler(), 
                                verbose=true,
                                keepall=true)
    X0OptOptions(x0_lb, x0_rb, searchspace_size, iters, sampler, verbose, keepall)
end


"""
    hypersearch

General macro to perform hyperparameter search.
- usage: 
    @hypesearch options::AbstractHyperoptOptions optimize(args; kwargs)

Currently supports searching x0(starting point) through putting a [`X0OptOptions`](@ref) at first position. 
If you want to have more choices, please implement more [`AbstractHyperoptOptions`](@ref) with corresponding [`hypersearch`](@ref) function. 
"""
macro hypersearch(search_expr)
    @assert valid_search_expr(search_expr) "Invalid hypersearch expression"
    options, optimize_expr = search_expr.args[1], search_expr.args[2]
    expr = get_hypersearch_expr(optimize_expr, options)
    esc(expr)
end


"""
    search_x0

A syntactic sugar that easily perform x0 search. 
- usage: 
    1. @search_x0 optimize(args; kwargs)    ------>    Searching x0 using default X0OptOptions
    2. @search_x0 options::AbstractHyperoptOptions optimize(args; kwargs)    ------>    Same as call @hypersearch
"""
macro search_x0(search_expr)
    @assert valid_search_expr(search_expr) ||  search_expr.head == :call "Invalid hypersearch expression"
    if search_expr.head == :call
        options, optimize_expr = :(X0OptOptions()), search_expr
    elseif search_expr.head == :tuple
        options, optimize_expr = search_expr.args[1], search_expr.args[2]
    end
    expr = get_hypersearch_expr(optimize_expr, options)
    esc(expr)
end


"""
    valid_search_expr

Check hyperparameter search expression passed to macro is valid or not. 
A valid expression should be of the form :(options::AbstractHyperoptOptions, optimize(args...; kwargs...))
"""
function valid_search_expr(search_expr)
    search_expr.head == :tuple && length(search_expr.args) == 2
end


"""
    get_hypersearch_expr

Concatenate hypersearch expression that will be returned by macro
"""
function get_hypersearch_expr(optimize_expr, options)
    @assert optimize_expr.args[1] == :(Nonconvex.optimize) or optimize_expr.args[1] == :(optimize) "Invalid hypersearch method. "
    expr = Expr(:call, :(Nonconvex.hypersearch), options, optimize_expr.args[2:end]...)
end

"""
    hypersearch

Main function that perform hyperparameter search
- `options::X0OptOptions`: An instance of ['X0OptOptions'](@ref) 
- `optim_args...`:  Arguments that will pass to `optimize`
- `optim_kwargs...`: keyword arguments that will pass to `optimize`
"""
function hypersearch(options::X0OptOptions, optim_args...; optim_kwargs...)
    # @show optim_args
    model = optim_args[1]
    @unpack x0_lb, x0_rb = options
    @unpack searchspace_size, iters, sampler = options
    @unpack verbose, keepall = options
    @info "Searching starting point... "
    # Warn if there is an x0 
    if optim_args[3] !== nothing
        @warn "Due to you are performing hyperparameter optimization on starting point, the x0 passed to function `optimize` will be ignored. "
    end
    # Process lower bound and upper bound, if the bounds are not specified, then use bounds of model instead
    @assert !(x0_lb === nothing && isempty(model.box_min)) || (x0_rb === nothing && isempty(model.box_max)) 
            "Model must have corresponding box constraint if you don't speficy search bound. "
    if x0_lb === nothing
        x0_lb = model.box_min
    end
    if x0_rb === nothing
        x0_rb = model.box_max
    end
    # Generate a search space using Sobol Sequence
    @assert searchspace_size >= 1 "searchspace_size must be a positive integer. "
    sobol_seq_generator = SobolSeq(x0_lb, x0_rb)
    searchspace = reduce(append!, [next!(sobol_seq_generator)] for _ = 1:searchspace_size)
    # Perform search on x0
    result_arr = AbstractResult[]
    print_fmt = string("Iteration: %4i, _x0: [", repeat("%-10.5f, ", length(searchspace[1]))[1:end-2], "], minimum: %-20.10f\n")
    ######################################
    # Printf has a bug that can not use predefined string format, (temporarily) define a closure to address this. 
    @eval printf_result(args...) = @printf($print_fmt, args...)
    ######################################
    ho = @hyperopt for i = iters, sampler = sampler, _x0 = searchspace
        _optim_args = (optim_args[begin:2]..., _x0)
        result = optimize(_optim_args...;optim_kwargs...)
        push!(result_arr, result)
        result.minimum
        if verbose === true
            ######################################
            # Invoke function in above closure
            Base.invokelatest(printf_result, i, _x0..., result.minimum)
            ######################################
        end
    end
    optima = reduce((result1, result2) -> result2.minimum < result1.minimum ? result2 : result1, result_arr)
    HyperoptResult(optima, ho, keepall ? result_arr : nothing)
end
