################################################################################
# SemistandardTableaux.jl
#
# A Semistandard Tableaux package for Julia.
#
# Copyright (C) 2026 Dave Anderson, pseudoeffective.github.io
################################################################################


module SemistandardTableaux

################################################################################
# Import
################################################################################

# Base
import Base:
	*, transpose

# Combinatorics
import Combinatorics:
	nthperm

# AbstractAlgebra
import AbstractAlgebra:
	base_ring, gen, gens, parent_type, nvars, polynomial_ring, MPolyBuildCtx, push_term!, finish

# Nemo
import Nemo:
	ZZ, QQ, libflint, ZZMPolyRing, ZZMPolyRingElem, evaluate, vars, coefficients

# LinearAlgebra for determinant
import LinearAlgebra: 
	det

# Memoization
import Memoization: 
	@memoize


# bpds
using BumplessPipeDreams

################################################################################
# Export (more exports are in the source files)
################################################################################

export
	ZZ, QQ, PolyRing, ZZMPolyRing, ZZMPolyRingElem, QQMPolyRing, QQMPolyRingElem
	



################################################################################
# source files
################################################################################

include("ssyt.jl")


end
