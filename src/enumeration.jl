# Enumeration of semistandard Young tableaux
# David Anderson, June 2024.
#
# A single in-place "advance" core (next_ssyt!) drives every consumer: the
# materializing `ssyt`, the lazy `ssyt_iterator`, and `schur_poly_single`
# (see schur.jl). Method based on the OSCAR version by Ulrich Thiel and
# collaborators.

export ssyt, ssyt_iterator


# make the superstandard tableau from partition shape la(/mu)
function sstab( la::Vector{Int}, mu=[]; rowmin=false )
  tab=[]

  if rowmin
    for i in 1:length(la)
      push!(tab, fill(i, la[i]) )
    end

    for i=1:length(mu)
      for j=1:mu[i]
        tab[i][j]=0
      end
    end
    return Tableau(tab)
  end


  for i in 1:length(la)
    tr=[]
    if i<=length(mu)
      push!(tr,fill(0,mu[i])...)
      j=mu[i]+1
    else
      j=1
    end
    while j<=la[i]
      k=1
      while k<i && tab[i-k][j]>0
        k+=1
      end
      push!(tr,k)
    j+=1
    end
    push!(tab,tr)
  end

  return Tableau(tab)
end


	#- - - - - - - - - - - - - - - - - - - -#
	#	Shared in-place enumeration core	#
	#- - - - - - - - - - - - - - - - - - - -#

"""
    SSYTState

Mutable cursor over the semistandard Young tableaux of a fixed (skew) shape and
flag.  Holds the current tableau (`tab`, mutated in place) and the bookkeeping
needed by `next_ssyt!`.  `mu` is stored extended to length `len`.
"""
mutable struct SSYTState
    lambda::Vector{Int}
    mu::Vector{Int}
    ff::Vector{Vector{Int}}
    rowmin::Bool
    tab::Vector{Vector{Int}}
    m::Int
    n::Int
    len::Int
    started::Bool
    live::Bool
end

# initialize from sstab; cursor sits on the first (super)standard filling.
function ssyt_state(lambda::Vector{Int}, ff::Vector{Vector{Int}};
                    mu::Vector{Int}=fill(0, length(lambda)), rowmin::Bool=false)
    len = length(lambda)
    tab = sstab(lambda, mu, rowmin=rowmin).t
    muext = vcat(mu, [0 for s = length(mu)+1:len])  # extend mu by 0 if necessary
    n = len == 0 ? 0 : lambda[len]
    return SSYTState(lambda, muext, ff, rowmin, tab, len, n, len, false, true)
end

# is the current tableau within the flag everywhere?
function _ssyt_valid(s::SSYTState)::Bool
    for i = 1:s.len
        for j = s.mu[i]+1:s.lambda[i]
            if s.tab[i][j] > s.ff[i][j]
                return false
            end
        end
    end
    return true
end

# advance the tableau in place to the next filling in enumeration order,
# regardless of flag validity; return false when the shape is exhausted.
# This is the original "raise one element / minimize trailing" logic, ONCE.
function _ssyt_advance!(s::SSYTState)::Bool
    s.len == 0 && return false

    lambda = s.lambda; mu = s.mu; ff = s.ff; tab = s.tab
    len = s.len; rowmin = s.rowmin
    m = s.m; n = s.n

    #raise one element by 1
    while !(tab[m][n] < ff[m][n] &&
            (n == lambda[m] || tab[m][n] < tab[m][n + 1]) &&
            (m == len || lambda[m + 1] < n || tab[m][n] + 1 < tab[m + 1][n]))
      if n > mu[m]+1
        n -= 1
      elseif m > 1
        m -= 1
        n = lambda[m]
      else
        s.m = m; s.n = n
        return false
      end
    end

    if tab[m][n]>0
       tab[m][n] += 1
    end

    #minimize trailing elements
    if n < lambda[m]
      i = m
      j = n + 1
    else
      i = m + 1
      i<=len ? j = mu[i]+1 : j=1
    end
    while (i <= len && j <= lambda[i])
      if i == 1
        tab[1][j] = tab[1][j - 1]  #if i==1 then j!=1 by initialization
      elseif j == 1
        if rowmin  # ensure tab entry is at least row index
           tab[i][1] = max(i,tab[i - 1][1] + 1) #likewise if j==1 then i!=1
        else
           tab[i][1] = tab[i - 1][1] + 1 #likewise if j==1 then i!=1
        end
      else
        if rowmin  # ensure tab entry is at least row index
           tab[i][j] = max(tab[i][j - 1], tab[i - 1][j] + 1, i) #likewise if j==1 then i!=1
        else
           tab[i][j] = max(tab[i][j - 1], tab[i - 1][j] + 1) #likewise if j==1 then i!=1
        end
      end
      if j < lambda[i]
        j += 1
      else
        i += 1
        i<=len ? j = mu[i]+1 : j=1
      end
    end

    s.m = len
    s.n = lambda[len]
    return true
end

"""
    next_ssyt!(s::SSYTState) -> Bool

Advance `s` in place to the next *valid* (flag-respecting) tableau, leaving it in
`s.tab`.  Returns `false` when the enumeration is exhausted.

Faithful to the original loop: validity is "sticky".  The first time the cursor
lands on a tableau that exceeds the flag, the enumeration stops — matching the
old code, where `valid` was set once and never reset.
"""
function next_ssyt!(s::SSYTState)::Bool
    s.live || return false

    if !s.started
        s.started = true            # cursor already sits on the initial filling
    elseif !_ssyt_advance!(s)
        s.live = false
        return false
    end

    if _ssyt_valid(s)
        return true
    else
        s.live = false              # sticky: stop at the first invalid tableau
        return false
    end
end


	#- - - - - - - - - - - - - - - - - - - -#
	#	Lazy iterator						#
	#- - - - - - - - - - - - - - - - - - - -#

"""
    ssyt_iterator(la, ff; mu=fill(0,length(la)), rowmin=false)

Lazy iterator over the semistandard Young tableaux of (skew) shape `la/mu`
bounded by the flag `ff`, mirroring OSCAR's `semistandard_tableaux`.  Each
element is a freshly copied `Tableau`; `collect` materializes the full vector
(this is exactly what `ssyt` returns).  `ff` accepts the same `Int`,
`Vector{Int}`, and `Vector{Vector{Int}}` forms as `ssyt`.
"""
struct SSYTIterator
    lambda::Vector{Int}
    ff::Vector{Vector{Int}}
    mu::Vector{Int}
    rowmin::Bool
end

function ssyt_iterator(lambda::Vector{Int}, ff::Vector{Vector{Int}};
                       mu::Vector{Int}=fill(0, length(lambda)), rowmin::Bool=false)
    return SSYTIterator(lambda, ff, mu, rowmin)
end

# flag normalizers
function ssyt_iterator(lambda, ff::Vector{Int}; mu=fill(0, length(lambda)), rowmin=false)
    bd = Vector{Vector{Int}}([fill(ff[i], lambda[i]) for i = 1:length(lambda)])
    return ssyt_iterator(lambda, bd; mu=mu, rowmin=rowmin)
end

function ssyt_iterator(lambda, ff::Int=length(lambda); mu=fill(0, length(lambda)), rowmin=false)
    bd = Vector{Vector{Int}}([fill(ff, lambda[i]) for i = 1:length(lambda)])
    return ssyt_iterator(lambda, bd; mu=mu, rowmin=rowmin)
end

Base.IteratorSize(::Type{SSYTIterator}) = Base.SizeUnknown()
Base.eltype(::Type{SSYTIterator}) = Tableau

function Base.iterate(it::SSYTIterator)
    s = ssyt_state(it.lambda, it.ff; mu=it.mu, rowmin=it.rowmin)
    return iterate(it, s)
end

function Base.iterate(it::SSYTIterator, s::SSYTState)
    next_ssyt!(s) || return nothing
    return (Tableau(deepcopy(s.tab)), s)
end


	#- - - - - - - - - - - - - - - - - - - -#
	#	Materializing constructor			#
	#- - - - - - - - - - - - - - - - - - - -#

"""
    ssyt(la::Vector{Int}, ff::Union{Int, Vector{Int}, Vector{Vector{Int}}}; mu::Vector{Int}=fill(0,length(la)), rowmin::Bool=false) -> Vector{Tableau}

Constructs semistandard Young tableaux on a shape `la`, with given flagging conditions `ff`, optionally skewed by a subshape `mu`.  Method based on the OSCAR version by Ulrich Thiel and collaborators.

This is `collect(ssyt_iterator(la, ff; mu, rowmin))`; use `ssyt_iterator` directly for lazy enumeration.

## Arguments
- `la::Vector{Int}`: A partition represented as a nonincreasing vector of integers.

- `ff::Union{Int, Vector{Int}, Vector{Vector{Int}}}`: Flagging conditions for the tableaux. If provided as an integer, it specifies the largest value allowed for an entry. If given as a single vector of integers, it specifies a row flagging condition. If provided as a vector of vectors, it specifies a filling which bounds the tableaux entrywise.

- `mu::Vector{Int}`: An optional keyword argument, giving a subshape of the partition `la`. Defaults to an empty vector for a straight shape.

- `rowmin::Bool`: An optional boolean keyword. When set to `true`, indicates the entries in row `i` must be at least `i`. Condition is redundant for straight-shape tableaux. Defaults to `false`.

## Returns
- `Vector{Tableau}`: A vector of `Tableau` objects, each representing a semistandard Young tableau that satisfies the given shape and flagging conditions.

# Examples
```julia-repl
# Generate SSYTs for the partition [3, 2] with largest entry 3
julia> tabs = ssyt([3, 2], 3)

# Generate SSYTs for the shape [3, 2, 1]/[1] and row-minimal condition
julia> skewtabs = ssyt([3, 2, 1], [3, 3, 3], mu=[1], rowmin=true)
```
"""
function ssyt(lambda::Vector{Int}, ff::Vector{Vector{Int}};
              mu::Vector{Int}=fill(0, length(lambda)), rowmin::Bool=false)
    return collect(ssyt_iterator(lambda, ff; mu=mu, rowmin=rowmin))
end

###
function ssyt(lambda, ff::Vector{Int}; mu=fill(0, length(lambda)), rowmin=false)
# case of flagged tableaux
    return collect(ssyt_iterator(lambda, ff; mu=mu, rowmin=rowmin))
end

###
function ssyt(lambda, ff::Int=length(lambda); mu=fill(0, length(lambda)), rowmin=false)
# uniformly bounded, default to number of rows of lambda
    return collect(ssyt_iterator(lambda, ff; mu=mu, rowmin=rowmin))
end
