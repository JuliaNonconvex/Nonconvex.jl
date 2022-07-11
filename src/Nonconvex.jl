module Nonconvex

using Reexport
@reexport using NonconvexCore
@reexport using NonconvexUtils

macro load(algo)
    esc(_load(string(algo)))
end
macro load(algos...)
    exprs = map(algos) do algo
        :(Nonconvex.@load $algo)
    end
    return Expr(:block, exprs...)
end

function _load(algo)
    if algo in ("MMA", "GCMMA", "MMA87", "MMA02")
        return install_and_load_module(:NonconvexMMA)
    elseif algo in ("Ipopt", "IpoptAlg")
        return install_and_load_module(:NonconvexIpopt)
    elseif algo in ("NLopt", "NLoptAlg")
        return install_and_load_module(:NonconvexNLopt)
    elseif algo in ("Juniper", "JuniperIpopt", "JuniperIpoptAlg")
        return install_and_load_module(:NonconvexJuniper)
    elseif algo in ("Pavito", "PavitoIpoptCbc", "PavitoIpoptCbcAlg")
        return install_and_load_module(:NonconvexPavito)
    elseif algo in ("Percival", "AugLag", "PercivalAlg")
        return install_and_load_module(:NonconvexPercival)
    elseif algo == "AugLag2"
        return install_and_load_module(:NonconvexAugLagLab)
    elseif algo in ("Bayesian", "BayesOpt", "BayesOptAlg")
        return install_and_load_module(:NonconvexBayesian)
    elseif algo in ("SDP", "Semidefinite", "SDPBarrier", "SDPBarrierAlg")
        return install_and_load_module(:NonconvexSemidefinite)
    elseif algo in ("Search", "MTS", "LS1", "MTSAlg", "LS1Alg")
        return install_and_load_module(:NonconvexSearch)
    elseif algo in ("Hyperopt", "Deflated", "Multistart", "HyperoptAlg", "DeflatedAlg")
        return install_and_load_module(:NonconvexMultistart)
    elseif algo == "TOBS"
        return install_and_load_module(:NonconvexTOBS)
    elseif algo == "Metaheuristics"
        return install_and_load_module(:NonconvexMetaheuristics)
    elseif algo == "NOMAD"
        return install_and_load_module(:NonconvexNOMAD)
    else
        throw("Unsupported algorithm. Please check the documentation of Nonconvex.jl.")
    end
end
function install_and_load_module(mod)
    quote
        using Pkg
        modname = $(Meta.quot(mod))
        try
            @info "Attempting to load the package $modname."
            using $mod
            @info "Loading succesful."
        catch
            @info "Couldn't find the package $modname. Attempting to install it."
            try
                Pkg.add(string(modname))
            catch
                @info "Package installation failed! Please report an issue."
                return
            end
            @info "$modname installed."
            @info "Attempting to load the package $modname."
            using $mod
            @info "Loading succesful."
        end
    end
end

end
