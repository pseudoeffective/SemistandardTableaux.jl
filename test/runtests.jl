using SemistandardTableaux
using Test

@testset "SemistandardTableaux.jl" begin

    # Phase 1 smoke tests: the package loads and the existing happy paths run.

    @testset "SSYT enumeration" begin
        # |SSYT(shape (2,1), entries <= 3)| = 8
        @test length(ssyt([2, 1], 3)) == 8
    end

    @testset "schur_poly happy paths" begin
        # default ring (with y-variables): double/factorial Schur polynomial
        p = schur_poly([2, 1], 3)
        @test !iszero(p)

        # x-only ring: ordinary single Schur polynomial
        R = ssyt_ring(3, 0)
        p1 = schur_poly([2, 1], 3, R)
        @test !iszero(p1)
        @test parent(p1) === R
    end

    @testset "deprecated xy_ring shim" begin
        # full deprecation-warning assertion lives in the Phase 4 suite;
        # here we just confirm the compat shim still returns (R, x, y).
        R, x, y = xy_ring(3, 2)
        @test length(x) == 3
        @test length(y) == 2
    end

end
