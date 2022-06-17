# Problem definition

There are 3 ways to define a model in Nonconvex.jl:
1. `Model` which assumes all the variables are indexed by an integer index starting from 1. The decision variables are therefore a vector.
2. `DictModel` which assumes each variable has a name. The decision variables are stored in an `OrderedDict`, an ordered dictionary data structure.
3. Start from `JuMP.Model` and convert it to `DictModel`. This is convenient to make use of `JuMP`'s user-friendly macros for variable and linear expression, objective or constraint definitions.

```@contents
Pages = ["model.md", "dict_model.md", "queries.md"]
Depth = 3
```
