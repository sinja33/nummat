using Test
using DN2
using LinearAlgebra
using QuadGK

# Kontrolni poligon iz naloge
const TOCKE = [(0,0),(1,1),(2,3),(1,4),(0,4),(-1,3),(0,1),(1,0)]
const B = Bezier(TOCKE)

@testset "DN2 — Ploščina zanke Bézierjeve krivulje" begin

    # -------------------------------------------------------
    # 1. Konstrukcija tipa Bezier
    # -------------------------------------------------------
    @testset "Konstrukcija Bezier" begin
        @test size(B.P) == (2, 8)
        @test B.P[:, 1] ≈ [0.0, 0.0]
        @test B.P[:, end] ≈ [1.0, 0.0]
    end

    # -------------------------------------------------------
    # 2. Evalvacija krivulje — robni pogoji
    # -------------------------------------------------------
    @testset "tocka — robni pogoji" begin
        # Bézierjeva krivulja vedno začne in konča v krajnih kontrolnih točkah
        @test tocka(B, 0.0) ≈ [0.0, 0.0] atol=1e-12
        @test tocka(B, 1.0) ≈ [1.0, 0.0] atol=1e-12
    end

    # -------------------------------------------------------
    # 3. Evalvacija krivulje — analitična rešitev
    # -------------------------------------------------------
    @testset "tocka — kvadratna krivulja (analitična rešitev)" begin
        # B(t) = (1-t)^2*(0,0) + 2t(1-t)*(1,2) + t^2*(2,0)
        # B(0.5) = 0.25*(0,0) + 0.5*(1,2) + 0.25*(2,0) = (1, 1)
        B2 = Bezier([(0.0,0.0),(1.0,2.0),(2.0,0.0)])
        @test tocka(B2, 0.5) ≈ [1.0, 1.0] atol=1e-12
        @test tocka(B2, 0.0) ≈ [0.0, 0.0] atol=1e-12
        @test tocka(B2, 1.0) ≈ [2.0, 0.0] atol=1e-12
    end

    # -------------------------------------------------------
    # 4. Odvod krivulje — linearna krivulja
    # -------------------------------------------------------
    @testset "odvod — linearna krivulja" begin
        # B(t) = (2t, 4t)  →  B'(t) = (2, 4) povsod
        Blin = Bezier([(0.0,0.0),(2.0,4.0)])
        for t in [0.0, 0.25, 0.5, 0.75, 1.0]
            @test odvod(Blin, t) ≈ [2.0, 4.0] atol=1e-12
        end
    end

    # -------------------------------------------------------
    # 5. Odvod krivulje — analitična rešitev
    # -------------------------------------------------------
    @testset "odvod — kvadratna krivulja (analitična rešitev)" begin
        # B'(0.5) = 2*(1-0.5)*[(1,2)-(0,0)] + 2*0.5*[(2,0)-(1,2)] = (2, 0)
        B2 = Bezier([(0.0,0.0),(1.0,2.0),(2.0,0.0)])
        @test odvod(B2, 0.5) ≈ [2.0, 0.0] atol=1e-12
    end

    # -------------------------------------------------------
    # 6. Gauss-Legendrova kvadratura — znani integrali
    # -------------------------------------------------------
    @testset "integriraj_gl — znani integrali" begin
        @test integriraj_gl(sin, 0.0, Float64(π), 50) ≈ 2.0   atol=1e-10
        @test integriraj_gl(x -> x^2, 0.0, 1.0, 10)  ≈ 1/3    atol=1e-12
        @test integriraj_gl(exp, 0.0, 1.0, 20)        ≈ ℯ - 1  atol=1e-12
        # Polinom stopnje 9 — GL s 5 točkami na enem intervalu mora biti točen
        @test integriraj_gl(x -> x^9, 0.0, 1.0, 1)   ≈ 1/10   atol=1e-12
    end

    # -------------------------------------------------------
    # 7. Iskanje zanke — krivulja z zanko
    # -------------------------------------------------------
    @testset "poisci_zanko — krivulja z zanko" begin
        rezultat = poisci_zanko(B)
        @test rezultat !== nothing
        t1, t2 = rezultat
        @test 0.0 ≤ t1 < t2 ≤ 1.0
        # Presečišče točno na 1e-10
        @test norm(tocka(B, t1) - tocka(B, t2)) < 1e-10
        # Simetrija: t1 + t2 ≈ 1 (kontrolni poligon je simetričen)
        @test t1 + t2 ≈ 1.0 atol=1e-8
    end

    # -------------------------------------------------------
    # 8. Iskanje zanke — krivulja brez zanke
    # -------------------------------------------------------
    @testset "poisci_zanko — krivulja brez zanke" begin
        # Enostavna kvadratna krivulja nima zanke
        Benostavna = Bezier([(0.0,0.0),(1.0,2.0),(2.0,0.0)])
        @test poisci_zanko(Benostavna) === nothing
    end

    # -------------------------------------------------------
    # 9. Ploščina zanke — konvergenca
    # -------------------------------------------------------
    @testset "ploscina_zanke — konvergenca" begin
        P1 = ploscina_zanke(B, n_intervalov=5)
        P2 = ploscina_zanke(B, n_intervalov=10)
        P3 = ploscina_zanke(B, n_intervalov=20)
        P4 = ploscina_zanke(B, n_intervalov=40)
        @test abs(P2 - P3) < abs(P1 - P2)
        @test abs(P3 - P4) < abs(P2 - P3)
        # Na koncu se stabilizira
        @test abs(P3 - P4) < 1e-8
    end

    # -------------------------------------------------------
    # 10. Ploščina zanke — primerjava z neodvisno metodo
    # -------------------------------------------------------
    @testset "ploscina_zanke — primerjava z QuadGK" begin
        t1, t2 = poisci_zanko(B)
        integrand = t -> begin
            xy  = tocka(B, t)
            dxy = odvod(B, t)
            xy[1] * dxy[2] - dxy[1] * xy[2]
        end
        P_ref, _ = quadgk(integrand, t1, t2, rtol=1e-12)
        P_ref = abs(P_ref) / 2
        @test ploscina_zanke(B) ≈ P_ref atol=1e-8
    end

    # -------------------------------------------------------
    # 11. Ploščina zanke — robni primer brez zanke
    # -------------------------------------------------------
    @testset "ploscina_zanke — brez zanke vrne 0" begin
        Benostavna = Bezier([(0.0,0.0),(1.0,2.0),(2.0,0.0)])
        @test ploscina_zanke(Benostavna) == 0.0
    end

end