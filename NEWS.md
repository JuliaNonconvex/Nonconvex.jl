# NEWS

## 4.0.0

- Support NonconvexCore 2 (`NonconvexCore = "2"` compat); the
  NonconvexUtils dependency is dropped since its utilities were absorbed
  into NonconvexCore 2.
- The documentation is now assembled from the documentation fragments of
  NonconvexCore and of all the algorithm packages. The Core section became a
  General section hosting the Core docs together with a new `@load` guide
  (accepted names, the packages they install and load), and the
  Gradients, Jacobians and Hessians section moved into NonconvexCore.
- The landing page is an actual overview (features, installation and a quick
  start) with a small centered logo; the same small-logo treatment applies
  to every package page in the assembled docs.
- The downstream integration tests run in a temporary environment, and the
  Installation matrix precompiles the test environment with the test flags
  instead of running through `Pkg.test`'s sandbox, which also fixes the
  Julia 1.13 pkgimage rejections.
- The generated `@load` code splices the captured `Pkg` module instead of
  loading `Pkg` in the caller, so it works inside `Pkg.test`-style sandboxes.
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
