using SemistandardTableaux
using AbstractAlgebra: gens, evaluate, zero
using Test

@testset "SemistandardTableaux.jl" begin

    @testset "SSYT cardinality" begin
        @test length(ssyt([2, 1], 3)) == 8
        @test length(ssyt([3, 2, 1], 7)) == 2352
        # returned tableaux are independent copies, all distinct
        tabs = ssyt([2, 1], 3)
        @test length(unique(tabs)) == 8
    end

    @testset "lazy iterator agrees with ssyt" begin
        it = ssyt_iterator([3, 2], 3)
        @test eltype(it) == Tableau
        @test collect(it) == ssyt([3, 2], 3)
        # lazy access works without materializing
        @test first(ssyt_iterator([2, 1], 3)) isa Tableau
        # flag-form normalizers reach the iterator too
        @test collect(ssyt_iterator([2, 1], [3, 3])) == ssyt([2, 1], [3, 3])
    end

    @testset "skew and flagged enumeration" begin
        # skew shape (2,2)/(1): cells (1,2),(2,1),(2,2), entries <= 3
        @test length(ssyt([2, 2], 3; mu = [1])) == 8
        # row flag row1<=1, row2<=2 on shape (2,1): only [[1,1],[2]]
        @test length(ssyt([2, 1], [1, 2])) == 1
        @test ssyt([2, 1], [1, 2])[1] == Tableau([[1, 1], [2]])
    end

    @testset "single Schur expansion (ordinary by default)" begin
        R = ssyt_ring(3, 0)
        x = gens(R)
        expected = x[1]^2*x[2] + x[1]^2*x[3] + x[1]*x[2]^2 + 2*x[1]*x[2]*x[3] +
                   x[1]*x[3]^2 + x[2]^2*x[3] + x[2]*x[3]^2
        # bare call defaults to the ordinary single Schur polynomial
        @test schur_poly([2, 1], 3) == expected
        # explicit x-only ring agrees
        @test schur_poly([2, 1], 3; ring=R) == expected
    end

    @testset "keyword ring/coeff API" begin
        # coeff=QQ path constructs and evaluates
        RQ = ssyt_ring(3, 0; coeff=QQ)
        xq = gens(RQ)
        expectedQ = xq[1]^2*xq[2] + xq[1]^2*xq[3] + xq[1]*xq[2]^2 + 2*xq[1]*xq[2]*xq[3] +
                    xq[1]*xq[3]^2 + xq[2]^2*xq[3] + xq[2]*xq[3]^2
        @test schur_poly([2, 1], 3; coeff=QQ) == expectedQ
        # empty shape is handled before any la[1] access
        @test schur_poly(Int[], 0) == one(ssyt_ring(0, 0))
    end

    @testset "factorial reduces to ordinary (y => 0)" begin
        R1 = ssyt_ring(3, 0)
        g1 = gens(R1)
        R2 = ssyt_ring(3, 5)                  # x1..x3 then y1..y5
        p2 = schur_poly([2, 1], 3; ring=R2)   # double/factorial version

        # substitute x_i -> x_i of R1, all y_j -> 0
        sub = [startswith(string(v), "x") ? g1[parse(Int, string(v)[2:end])] : zero(R1)
               for v in gens(R2)]
        @test evaluate(p2, sub) == schur_poly([2, 1], 3; ring=R1)
    end

    @testset "deprecated xy_ring shim" begin
        R, x, y = xy_ring(3, 2)
        @test length(x) == 3
        @test length(y) == 2
    end

    @testset "insertion algorithms" begin
        T0 = Tableau([[1, 1, 3], [2, 4], [5, 5]])
        @test rsk_insert(T0, 2) == Tableau([[1, 1, 2], [2, 3], [4, 5], [5]])
    end

end
