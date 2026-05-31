# Polynomial ring helpers for SemistandardTableaux.jl
# Replaces the old DoublePolyRing abstraction with a plain MPolyRing
# plus name-based variable extraction.

export ssyt_ring, extract_vars

"""
    extract_vars(R::MPolyRing; varname::Symbol) -> Vector

Return the generators of `R` whose printed name begins with `varname`, in the order
`gens(R)` returns them. NOTE: relies on the ring having been built with the intended
generators in index order (x1,x2,...,xn then y1,...,ym). Prefix match means `varname=:x`
would also catch any `x`-prefixed family; fine for this package's rings.
"""
function extract_vars(R::MPolyRing; varname::Symbol)
    pre = string(varname)
    z = elem_type(R)[]
    for v in gens(R)
        startswith(string(v), pre) && push!(z, v)
    end
    return z
end

"""
    ssyt_ring(n, m=0; coeff=ZZ, xname=:x, yname=:y) -> MPolyRing

Polynomial ring with variables x1..xn then y1..ym over `coeff` (default Nemo ZZ).
Use `extract_vars` to recover the x- and y-families.
"""
function ssyt_ring(n::Int, m::Int=0; coeff=ZZ, xname::Symbol=:x, yname::Symbol=:y)
    names = vcat(["$(xname)$(i)" for i in 1:n], ["$(yname)$(j)" for j in 1:m])
    R, _ = polynomial_ring(coeff, names)
    return R
end

# --- deprecation shim: keep old call sites limping for one minor version ---
"""
    xy_ring(args...; kwargs...)  (DEPRECATED)

Returns `(R, x_vars, y_vars)`. The `DoublePolyRing` type is removed; migrate
`R.x_vars` -> `extract_vars(R; varname=:x)` and `R.y_vars` -> `extract_vars(R; varname=:y)`.
"""
function xy_ring(n::Int, m::Int=0; kwargs...)
    Base.depwarn("`xy_ring` and `DoublePolyRing` are deprecated; use `ssyt_ring` " *
                 "(returns an MPolyRing) with `extract_vars`.", :xy_ring)
    R = ssyt_ring(n, m; kwargs...)
    return R, extract_vars(R; varname=:x), extract_vars(R; varname=:y)
end
export xy_ring  # keep exported during the deprecation window
