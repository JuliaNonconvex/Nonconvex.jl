"""
    PlottingCallback

A callback function that plots a number of convergence metrics using [ConvergencePlots.jl](https://github.com/mohamed82008/ConvergencePlots.jl). The following are the fields of the struct:
 - `plot`: an instance of `ConvergencePlots.ConvergencePlot`.
 - `n`: the number of history points to keep track of.
 - `iter`: the number of times the callback was called.
"""
mutable struct PlottingCallback <: Function
    plot::Union{Nothing, ConvergencePlot}
    n::Int
    iter::Int
end
PlottingCallback(n = 1000) = PlottingCallback(nothing, n, 0)

"""
    (callback::PlottingCallback)(solution::Solution)

Updates the convergence plots using the current solution.
"""
function (callback::PlottingCallback)(solution::Solution)
    @unpack convstate = solution
    if callback.iter == 0
        callback.plot = ConvergencePlot(
            callback.n,
            names = ["KKT residual", "|Δf|", "|Δx|_∞", "Infeasibility"],
            #names = ["KKT residual", "IPOPT residual", "|Δf|", "|Δx|_∞", "Infeasibility"],
            options = Dict(
                "KKT residual" => (color = "black",),
                #"IPOPT residual" => (color = "blue",),
                "|Δf|" => (color = "green",),
                "|Δx|_∞" => (color = "violet",),
                "Infeasibility" => (color = "red",),
            ),
        )
    else
        @unpack plot = callback
        addpoint!(
            plot,
            Dict(
                "KKT residual" => convstate.kkt_residual,
                #"IPOPT residual" => convstate.ipopt_residual,
                "|Δf|" => convstate.Δf,
                "|Δx|_∞" => convstate.Δx,
                "Infeasibility" => convstate.infeas,
            ),
        )
    end
    callback.iter += 1
end


"""
    LazyPlottingCallback

A callback function that plots a number of convergence metrics using [ConvergencePlots.jl](https://github.com/mohamed82008/ConvergencePlots.jl) only when `show` is called. The following are the fields of the struct:
 - `plot`: an instance of `ConvergencePlots.ConvergencePlot`.
 - `n`: the number of history points to keep track of.
 - `iter`: the number of times the callback was called.
"""
mutable struct LazyPlottingCallback <: Function
    plot::Union{Nothing, ConvergencePlot}
    n::Int
    iter::Int
    show_plot::Bool
    save_plot::Union{Nothing, String}
end
LazyPlottingCallback(n = 10^10; show_plot = true, save_plot = nothing) = LazyPlottingCallback(nothing, n, 0, show_plot, save_plot)

"""
    (callback::LazyPlottingCallback)(solution::Solution)

Updates the convergence plots using the current solution.
"""
function (callback::LazyPlottingCallback)(solution::Solution; update = false)
    @unpack convstate = solution
    if callback.iter == 0
        callback.plot = ConvergencePlot(
            callback.n,
            names = ["KKT residual", "|Δf|", "|Δx|_∞", "Infeasibility"],
            #names = ["KKT residual", "IPOPT residual", "|Δf|", "|Δx|_∞", "Infeasibility"],
            options = Dict(
                "KKT residual" => (color = "black",),
                #"IPOPT residual" => (color = "blue",),
                "|Δf|" => (color = "green",),
                "|Δx|_∞" => (color = "violet",),
                "Infeasibility" => (color = "red",),
            ),
            show = update && callback.show_plot,
        )
    else
        @unpack plot = callback
        addpoint!(
            plot,
            Dict(
                "KKT residual" => convstate.kkt_residual,
                #"IPOPT residual" => convstate.ipopt_residual,
                "|Δf|" => convstate.Δf,
                "|Δx|_∞" => convstate.Δx,
                "Infeasibility" => convstate.infeas,
            ),
            show = callback.show_plot,
            update = update,
            filename = callback.save_plot,
        )
    end
    callback.iter += 1
end

struct NoCallback <: Function end
(::NoCallback)(args...; kwargs...) = nothing