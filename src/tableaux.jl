# Tableau type and basic operations
# David Anderson, June 2024.

	#- - - - - - - - - - - - - - - - - - - -#
	#	Type and Base methods				#
	#- - - - - - - - - - - - - - - - - - - -#

export Tableau

struct Tableau
  t::Vector{Vector{Int}}
end


# extract the (skew) shape of a tableau
function shape(tab::Tableau)::Tuple{Vector{Int},Vector{Int}}
  la=Int[]
  mu=Int[]

  len=length(tab.t)

  for i in 1:len
    push!(la,length(tab.t[i]))
    j=1
    while j<=length(tab.t[i]) && tab.t[i][j]==0
      j+=1
    end
    if j>1
      push!(mu,j-1)
    end
  end

  return la,mu
end


		# * * Base Methods * * #

# overload show to display Tableau
function Base.show(io::IO, tab::Tableau)
    println(io)
    for i in 1:length(tab.t)
        for j in 1:length(tab.t[i])
            tt=tab.t[i][j]
            if tt==0
              print(io, 'x', " ")
            else
              print(io, tt, " ")
            end
        end
        println(io)
    end
end

# overload identity for Tableau type
Base.:(==)(tab1::Tableau,tab2::Tableau) = tab1.t==tab2.t


	# * * Tableaux operations * * #

# returns the transpose of the Tableau T
function Base.transpose(T::Tableau)
	t = []
	if length(T.t) > 0
		for i in 1:length(T.t[1]) push!(t,[]) end
		for row in T.t for j in 1:length(row)
			push!(t[j],row[j])
		end end
	end
	return Tableau(t)
end
