# NEWS

## 4.0.0

- Support NonconvexCore 2 (`NonconvexCore = "2"` compat); the
  NonconvexUtils dependency is dropped since its utilities were absorbed
  into NonconvexCore 2.
- The installation/group tests load Ipopt before constructing IpoptAlg
  sub-solvers, and the group testsets run independently (GROUP=All runs
  every group).
- Aqua and JET quality gates added to the test targets.

## 3.2.0

- `Nonconvex.@load` now accepts `RCD` and `RCDAlg` as aliases for the
  random coordinate descent meta-algorithm in NonconvexSearch.jl. The
  `mts.md` docs page documents `RCDAlg` and `RCDOptions`.

## 3.1.1

- The `@load` macro only triggers `Pkg.add` for genuinely missing packages;
  other load failures (version conflicts, precompilation errors, package bugs)
  are rethrown instead of mutating the active environment.
- Unsupported algorithms throw a proper `ErrorException` instead of a bare
  `String`.
- The `modname` binding no longer leaks into the caller's scope, and the
  vararg `@load` form works without a `Nonconvex` binding in the caller.

## 3.1.0

- Requires Julia 1.10+ (matching NonconvexCore 1.6 / NonconvexUtils 1.0).
- Widened NonconvexCore and NonconvexUtils compatibility bounds.
