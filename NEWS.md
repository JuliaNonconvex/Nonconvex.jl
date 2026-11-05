# NEWS

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
