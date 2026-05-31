# Schur polynomials from semistandard tableaux
# David Anderson, June 2024.

	# * * Schur Polynomials * * #

export schur_poly

# return the product of binomials for a single tableau
function tab2bin( tab::Tableau, x, y; ring, xoffset = 0, yoffset = 0 )
  len = length(tab.t)

  n = length(x)
  m = length(y)

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
function ssyt2pol( tabs, x, y; ring, xoffset=0, yoffset=0 )

  pol=zero(ring)

  for tab in tabs
    pol = pol + tab2bin( tab, x, y; ring=ring, xoffset=xoffset, yoffset=yoffset )
  end

  return pol

end



"""
    schur_poly(la, ff, R=ssyt_ring(length(la), length(la)+la[1]); mu=[], xoffset=0, yoffset=0, rowmin=false)

Compute the Schur polynomial corresponding to a given (skew) partition `la/mu` and a flag `ff`, in an optionally specified polynomial ring `R`. The polynomial is constructed as an enumerator of semistandard Young tableaux of (skew) shape `la/mu` and bounded by the flagging condition `ff`.

## Arguments
- `la::Vector{Int}`: A partition represented as a vector of integers, specifying the shape of the Young diagram.
- `ff::Union{Int,Vector{Int},Vector{Vector{Int}}}`: A flag specifying bounds on the tableaux. If `ff` is given as a single integer, it bounds the entries of the tableaux.  If `ff` is a vector of integers, it must be of length at least that of `la`; then `ff[i]` bounds the entries in the `i`th row of the tableaux.  If `ff` is a vector of vectors, it is interpreted as a tableaux whose shape is assumed to contain `la`; then the entries of `ff` bound the tableaux entrywise.
- `R::MPolyRing`: An optional argument specifying the polynomial ring to use for constructing the Schur polynomial. The `x`- and `y`-variable families are recovered by name via `extract_vars`; if there are no `y`-variables the ordinary (single) Schur polynomial is returned, otherwise the double/factorial version. Defaults to a ring constructed by `ssyt_ring` based on the size of `la`.
- `mu::Vector{Int}`: An optional argument specifying a subpartition of `la`, for skew Schur polynomials. Defaults to an empty vector, for the straight shape `la`.
- `xoffset::Int`: An optional argument specifying an offset value for the x-variable indices in the polynomial. Defaults to 0.
- `yoffset::Int`: An optional argument specifying an offset value for the y-variable indices in the polynomial. Defaults to 0.
- `rowmin::Bool`: An optional argument specifying whether to use row-minimal tableau, i.e., to require that entries in row `i` be at least `i`.  (This is a nontrivial condition only for skew shapes.) Defaults to `false`.

## Returns
- An element of the polynomial ring `R` (e.g. `ZZMPolyRingElem`): the Schur polynomial.

# Examples
```julia-repl
# Specify a partition
julia> la = [2, 1]

# Specify a bound for the x-variables
julia> ff = 3

# Compute the Schur polynomial
julia> poly = schur_poly(la, ff)


### To get the single Schur polynomial, use a ring with no y-variables
julia> R = ssyt_ring(3, 0);

julia> poly1 = schur_poly(la, ff, R)
x1^2*x2 + x1^2*x3 + x1*x2^2 + 2*x1*x2*x3 + x1*x3^2 + x2^2*x3 + x2*x3^2

```
"""
function schur_poly( la, ff::Vector{Vector{Int}}, R::MPolyRing=ssyt_ring( length(la) , length(la)+la[1] ); mu = Int[], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return one(R)
  end

  x = extract_vars(R; varname=:x)
  y = extract_vars(R; varname=:y)

  if length(y)==0
     return schur_poly_single( la, ff, R; mu=mu, x=x, xoffset=xoffset, rowmin=rowmin )
  end

  tbs = ssyt( la, ff, mu=mu, rowmin=rowmin )

  pol = ssyt2pol( tbs, x, y; ring=R, xoffset=xoffset, yoffset=yoffset )

  return pol

end

###
function schur_poly( la, ff::Vector{Int}, R::MPolyRing=ssyt_ring( length(la) , length(la)+la[1] ); mu = Int[], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return one(R)
  end

  return schur_poly( la, Vector{Vector{Int}}([fill(ff[i],la[i]) for i=1:length(la)]), R; mu = mu, xoffset=xoffset, yoffset=yoffset, rowmin=rowmin )

end

###
function schur_poly( la, ff::Int=length(la), R::MPolyRing=ssyt_ring( length(la) , length(la)+la[1] ); mu = Int[], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return one(R)
  end

  return schur_poly( la, Vector{Int}(fill(ff,length(la))), R; mu = mu, xoffset=xoffset, yoffset=yoffset, rowmin=rowmin )
end


# faster polynomial constructor for single (one-set-of-variables) Schur
# polynomials: stream weights straight off the shared enumeration core into a
# build context, with no intermediate Tableau allocation.
function schur_poly_single(lambda::Vector{Int}, ff::Vector{Vector{Int}}, R::MPolyRing=ssyt_ring(max(max(ff...)...),0); mu::Vector{Int}=Int[], x=extract_vars(R; varname=:x), xoffset::Int=0, rowmin::Bool=false)
  if isempty(lambda)
    return one(R)
  end

  S = base_ring(R)
  sf = MPolyBuildCtx(R)
  xx = x

  s = ssyt_state(lambda, ff; mu=mu, rowmin=rowmin)
  count = zeros(Int, length(xx))

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
