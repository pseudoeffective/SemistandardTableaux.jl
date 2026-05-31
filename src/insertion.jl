# Insertion algorithms: RSK and Edelman-Greene
# David Anderson, June 2024.
#
# These are ring-independent operations on the Tableau type.

	# * * Insertion Algorithms * * #

export rsk_insert, rsk

# RSK insert the biletter (q,p) into the bitableau (P,Q).
function rsk_insert((P,Q)::Tuple{Tableau,Tableau},(q,p)::Tuple{Integer,Integer})
	local PP = deepcopy(P)
	local QQ = deepcopy(Q)
	local pp = p

	local i = 1
	while true
		if i > length(PP.t)
			push!(PP.t,[pp])
			push!(QQ.t,[q])
			break
		elseif pp >= last(PP.t[i])
			push!(PP.t[i],pp)
			push!(QQ.t[i],q)
			break
		end

		for j in 1:length(PP.t[i])
			if pp < PP.t[i][j]
				local temp = PP.t[i][j]
				PP.t[i][j] = pp
				pp = temp
				break
			end
		end

		i += 1
	end

	return (PP,QQ)
end

# RSK insert the letter p into the tableau P
function rsk_insert(P::Tableau,p::Integer)
	return rsk_insert((P,P),(1,p))[1]
end

# RSK insert a biword into the bitableau (P,Q).
function rsk_insert((P,Q)::Tuple{Tableau,Tableau},biword::Vector{<:Tuple{Vararg{Integer}}})
	local PP = deepcopy(P)
	local QQ = deepcopy(Q)
	for (q,p) in sort(biword)
		(PP,QQ) = rsk_insert((PP,QQ),(q,p))
	end
	return (PP,QQ)
end

# RSK insert a word into the tableau P.
function rsk_insert(P::Tableau,word::Vector{<:Integer})
	local PP = deepcopy(P)
	for p in word
		PP = rsk_insert(PP,p)
	end
	return PP
end

# sends a biword to its RSK tableau (P,Q).
function rsk(biword::Vector{<:Tuple{Vararg{Integer}}})
	local P = Tableau([])
	local Q = Tableau([])
	for (q,p) in sort(biword)
		(P,Q) = rsk_insert((P,Q),(q,p))
	end
	return (P,Q)
end

# sends a word to its RSK recording insertion P
function rsk(word::Vector{<:Integer})
	return rsk([(i,word[i]) for i in 1:length(word)])
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--

	julia> T0 = Tableau([[1,1,3],[2,4],[5,5]])
	1 1 3
	2 4
	5 5

	julia> T1 = rsk_insert(T0,2)
	1 1 2
	2 3
	4 5
	5

	julia> T2 = rsk_insert(T1,[4,3,4])
	1 1 2 3 4
	2 3 4
	4 5
	5

	julia> biword = [(1,2),(1,5),(2,3),(2,5),(3,1),(3,1),(3,4),(3,5),(4,3),(5,4)];

	julia> (P,Q) = rsk(biword)
	(
	1 1 3 4
	2 3 4 5
	5 5
	,
	1 1 2 3
	2 3 3 5
	3 4
	)

	julia> rsk_insert((P,Q),(6,2))
	(
	1 1 2 4
	2 3 3 5
	4 5
	5
	,
	1 1 2 3
	2 3 3 5
	3 4
	6
	)

+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#



export edelman_greene_insert, edelman_greene

# Edelman-Greene insert the biletter (q,p) into the bitableau (P,Q)
function edelman_greene_insert((P,Q)::Tuple{Tableau,Tableau},(q,p)::Tuple{Integer,Integer})
	local Pt = transpose(P)
	local Qt = transpose(Q)
	local pp = p

	local j = 1
	while true
		if j > length(Pt.t)
			push!(Pt.t,[pp])
			push!(Qt.t,[q])
			break
		elseif pp > last(Pt.t[j])
			push!(Pt.t[j],pp)
			push!(Qt.t[j],q)
			break
		end

		for i in 1:length(Pt.t[j])
			if pp < Pt.t[j][i]
				if Pt.t[j][i]==pp+1 && i>1 if Pt.t[j][i-1]==pp
					pp += 1
					break
				end end

				local temp = Pt.t[j][i]
				Pt.t[j][i] = pp
				pp = temp
				break
			end
		end

		j += 1
	end

	return (transpose(Pt),transpose(Qt))
end

# Edelman-Greene insert a biword into the bitableau (P,Q)
function edelman_greene_insert((P,Q)::Tuple{Tableau,Tableau},biword::Vector{<:Tuple{Vararg{Integer}}})
	local PP = deepcopy(P)
	local QQ = deepcopy(Q)
	for (q,p) in biword
		(PP,QQ) = edelman_greene_insert((PP,QQ),(q,p))
	end
	return (PP,QQ)
end

function edelman_greene(biword::Vector{<:Tuple{Vararg{Integer}}})
	local P = Tableau([])
	local Q = Tableau([])
	for (q,p) in biword
		(P,Q) = edelman_greene_insert((P,Q),(q,p))
	end
	return (P,Q)
end

#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =+
	--Example--

	julia> biword = [(1,4),(1,2),(2,5),(2,3),(2,2),(3,4),(3,3)];

	julia> (P,Q) = edelman_greene_insert(biword)
	(
	2 3 4
	3 4 5
	4
	,
	1 1 2
	2 2 3
	3
	)

+= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =#
