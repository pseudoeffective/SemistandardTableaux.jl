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

# AbstractAlgebra
import AbstractAlgebra:
	base_ring, gens, polynomial_ring, MPolyBuildCtx, push_term!, finish,
	elem_type, MPolyRing, one, zero

# Nemo (coefficient rings)
import Nemo:
	ZZ, QQ

################################################################################
# Export (more exports are in the source files)
################################################################################

export
	ZZ, QQ



################################################################################
# source files
################################################################################

include("rings.jl")
include("ssyt.jl")


end
