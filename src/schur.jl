# Schur polynomials from semistandard tableaux
# David Anderson, June 2024.

	# * * Schur Polynomials * * #

export schur_poly

# return the product of binomials for a single tableau
function tab2bin( tab::Tableau; double, ring, xoffset = 0, yoffset = 0 )
  len = length(tab.t)

  x=extract_vars(ring; varname=:x)
  y=extract_vars(ring; varname=:y)

  n = length(x)
  double ? m = length(y) : m = 0 # if not double, ignore y-variables even if present

  bin = one(ring)

  for i=1:len
    for j=1:length( tab.t[i] )
      tt = tab.t[i][j]
      if tt>0
        p=zero(ring)
        if tt+xoffset<=n
          p=p+x[tt+xoffset]
        end
        if tt+j-i+yoffset<=m && tt+j-i+yoffset>0
          p=p+y[tt+j-i+yoffset]
        end
        bin = bin*p
      end
    end
  end

  return bin

end


# sum of binomials for a set of tableaux
function ssyt2pol( tabs; double=false, ring, xoffset=0, yoffset=0 )

  pol=zero(ring)

  for tab in tabs
    pol = pol + tab2bin( tab; double=double, ring=ring, xoffset=xoffset, yoffset=yoffset )
  end

  return pol

end



"""
    schur_poly(la, ff=length(la); double=false, ring=nothing, coeff=ZZ, mu=Int[], xoffset=0, yoffset=0, rowmin=false)

Compute the Schur polynomial corresponding to a given (skew) partition `la/mu` and a flag `ff`. The polynomial is constructed as an enumerator of semistandard Young tableaux of (skew) shape `la/mu` and bounded by the flagging condition `ff`.

By default (`ring=nothing`, `double=false`) the **ordinary** (single) Schur polynomial is returned, in `x`-variables `x1..xN` over `coeff` (default Nemo `ZZ`), where `N` is the largest flag entry. To obtain the **double/factorial** Schur polynomial, pass a `ring` that also contains `y`-variables (e.g. from `ssyt_ring(n, m)` with `m > 0`); the `x`- and `y`-families are recovered by name via `extract_vars`.

## Arguments
- `la::Vector{Int}`: A partition represented as a vector of integers, specifying the shape of the Young diagram.
- `ff::Union{Int,Vector{Int},Vector{Vector{Int}}}`: A flag specifying bounds on the tableaux. If `ff` is given as a single integer, it bounds the entries of the tableaux.  If `ff` is a vector of integers, it must be of length at least that of `la`; then `ff[i]` bounds the entries in the `i`th row of the tableaux.  If `ff` is a vector of vectors, it is interpreted as a tableaux whose shape is assumed to contain `la`; then the entries of `ff` bound the tableaux entrywise. Defaults to `length(la)` (the number of rows).

## Keywords
- `double::Bool`: When `true`, compute the double/factorial Schur polynomial. If `ring` is not supplied, defaults to `ssyt_ring(length(la), length(la)+la[1])`. Default `false`.
- `ring::Union{Nothing,MPolyRing}`: The polynomial ring to build the answer in. `nothing` (default) builds an `x`-only ring for the ordinary Schur polynomial; supply a ring with `y`-variables for the double/factorial version.
- `coeff`: The coefficient ring used when `ring` is not supplied. Defaults to Nemo `ZZ`.
- `mu::Vector{Int}`: A subpartition of `la`, for skew Schur polynomials. Defaults to the empty vector (straight shape `la`).
- `xoffset::Int`, `yoffset::Int`: Offsets for the `x`- and `y`-variable indices. Default `0`.
- `rowmin::Bool`: When `true`, require entries in row `i` to be at least `i` (nontrivial only for skew shapes). Default `false`.

## Returns
- An element of the polynomial ring (e.g. `ZZMPolyRingElem`): the Schur polynomial.

# Examples
```julia-repl
# Ordinary Schur polynomial of shape [2,1] with entries up to 3
julia> schur_poly([2, 1], 3)
x1^2*x2 + x1^2*x3 + x1*x2^2 + 2*x1*x2*x3 + x1*x3^2 + x2^2*x3 + x2*x3^2

# Over the rationals
julia> schur_poly([2, 1], 3; coeff=QQ)
x1^2*x2 + x1^2*x3 + x1*x2^2 + 2*x1*x2*x3 + x1*x3^2 + x2^2*x3 + x2*x3^2

# Double/factorial version: supply a ring with y-variables
julia> R = ssyt_ring(3, 5);

julia> schur_poly([2, 1], 3; ring=R);
```
"""
function schur_poly( la, ff::Vector{Vector{Int}};
                     double::Bool=false,
                     coeff=ZZ,
                     ring::MPolyRing=(double ? ssyt_ring(length(la), length(la)+la[1]; coeff=coeff) : ssyt_ring(length(la), 0; coeff=coeff)),
                     mu::Vector{Int}=Int[],
                     xoffset::Int=0, yoffset::Int=0, rowmin::Bool=false )
  if length(la)==0
    return one(ring)
  end

  if !double
    # default: ordinary single Schur polynomial in x-variables over `coeff`
    return schur_poly_single( la, ff; ring=ring, coeff=coeff, mu=mu, xoffset=xoffset, rowmin=rowmin )
  end

  x = extract_vars(ring; varname=:x)
  y = extract_vars(ring; varname=:y)

  if length(y)==0
     return schur_poly_single( la, ff; ring=ring, x=x, mu=mu, xoffset=xoffset, rowmin=rowmin )
  end

  tbs = ssyt( la, ff; mu=mu, rowmin=rowmin )

  return ssyt2pol( tbs; double=double, ring=ring, xoffset=xoffset, yoffset=yoffset )
end

###
schur_poly( la, ff::Vector{Int}; kwargs... ) =
  schur_poly( la, Vector{Vector{Int}}([fill(ff[i], la[i]) for i=1:length(la)]); kwargs... )

###
schur_poly( la, ff::Int=length(la); kwargs... ) =
  schur_poly( la, Vector{Int}(fill(ff, length(la))); kwargs... )


# faster polynomial constructor for single (one-set-of-variables) Schur
# polynomials: stream weights straight off the shared enumeration core into a
# build context, with no intermediate Tableau allocation.
function schur_poly_single(lambda::Vector{Int}, ff::Vector{Vector{Int}};
                           ring::MPolyRing=ssyt_ring(length(lambda), 0; coeff=coeff), coeff=ZZ,
                           mu::Vector{Int}=Int[], x=nothing, xoffset::Int=0, rowmin::Bool=false)
  if isempty(lambda)
    return one(ring)
  end

  R = ring
  xx = isnothing(x) ? extract_vars(R; varname=:x) : x

  S = base_ring(R)
  sf = MPolyBuildCtx(R)

  s = ssyt_state(lambda, ff; mu=mu, rowmin=rowmin)
  # exponent vectors must span every variable of `R` (the ring may also carry
  # y-variables we are ignoring for the single polynomial); x are the first gens
  count = zeros(Int, length(gens(R)))

  while next_ssyt!(s)
    count .= 0
    for i = 1:s.len
      for j = s.mu[i]+1:s.lambda[i]
        v = s.tab[i][j] + xoffset
        if v <= length(xx)
          count[v] += 1
        end
      end
    end
    push_term!(sf, S(1), count)
  end

  return finish(sf)
end
