using Test, Nonconvex

@test_throws ArgumentError using NonconvexIpopt
Nonconvex.@load Ipopt
IpoptAlg()

@test_throws ArgumentError using NonconvexNLopt
Nonconvex.@load NLopt
NLoptAlg(:LD_MMA)

@test_throws ArgumentError using NonconvexJuniper
Nonconvex.@load Juniper
JuniperIpoptAlg()
IpoptAlg()

@test_throws ArgumentError using NonconvexMMA
Nonconvex.@load MMA
MMA87()
MMA02()

@test_throws ArgumentError using NonconvexPavito
Nonconvex.@load Pavito
PavitoIpoptCbcAlg()
IpoptAlg()

@test_throws ArgumentError using NonconvexPercival
Nonconvex.@load AugLag
AugLag()

@test_throws ArgumentError using NonconvexBayesian
Nonconvex.@load Bayesian
BayesOptAlg(IpoptAlg())

#@test_throws ArgumentError using NonconvexAugLagLab
#Nonconvex.@load AugLag2
#AugLag2()

@test_throws ArgumentError using NonconvexSemidefinite
Nonconvex.@load Semidefinite
SDPBarrierAlg(IpoptAlg())

@test_throws ArgumentError using NonconvexSearch
Nonconvex.@load Search
MTSAlg()
LS1Alg()

@test_throws ArgumentError using NonconvexMultistart
Nonconvex.@load Multistart
HyperoptAlg(IpoptAlg())
DeflatedAlg(IpoptAlg())
