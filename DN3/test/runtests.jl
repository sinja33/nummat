using Test
using DN3
using LinearAlgebra
using SpecialFunctions

# Referenčne vrednosti iz SpecialFunctions.jl
const Ai0_ref  = airyai(0.0)
const Aip0_ref = airyaiprime(0.0)

# Znane ničle (DLMF / Abramowitz & Stegun)
const ZNANE_NICLE = [
    -2.3381074104597670,
    -4.0879494441309706,
    -5.5205598280955260,
    -6.7867080900717660,
    -7.9441335871654310,
    -9.0226508533409880,
    -10.0401743415671360,
    -11.0085243037332600,
    -11.9360155632362590,
    -12.8287767528657330,
]

@testset "DN3 — Airyjeva funkcija" begin

    @testset "Začetni pogoji" begin
        Ai_val, Aip_val = airy_vrednost(0.0)
        @test abs(Ai_val  - Ai0_ref)  < 1e-10
        @test abs(Aip_val - Aip0_ref) < 1e-10
    end

    @testset "Magnusov korak — dimenzije" begin
        y0 = [Ai0_ref, Aip0_ref]
        y1 = magnus_korak(0.0, 0.1, y0)
        @test length(y1) == 2
        @test eltype(y1) == Float64
    end

    @testset "Vrednost v točki (primerjava z referenco)" begin
        # Primerjava z airyai iz SpecialFunctions.jl
        Ai1, _ = airy_vrednost(1.0; h=0.001)
        @test abs(Ai1 - airyai(1.0)) < 1e-6

        Aim1, _ = airy_vrednost(-1.0; h=0.001)
        @test abs(Aim1 - airyai(-1.0)) < 1e-5
    end

    @testset "Bisekcija — konvergenca" begin
        n, it, hist = bisekcija(sin, 3.0, 3.5)
        @test abs(n - π) < 1e-10
        @test it < 60

        napake = abs.(hist .- π)
        @test napake[end] < napake[1] / 2^(length(hist) - 1) * 2
    end

    @testset "Bisekcija — zahteva različna predznaka" begin
        @test_throws AssertionError bisekcija(x -> x^2 + 1, 0.0, 1.0)
    end

    @testset "Newtonova metoda — konvergenca" begin
        n, it, hist = newton_raphson(sin, cos, 3.0)
        @test abs(n - π) < 1e-12
        @test it < 10
    end

    @testset "Ničle Airyjeve funkcije" begin
        nicle, rez_b, rez_n = poisci_nicle(5; h=0.01, x_min=-15.0)

        for i in 1:min(3, length(nicle))
            @test abs(nicle[i] - ZNANE_NICLE[i]) < 1e-5
        end

        for nb in nicle[1:min(3, end)]
            val, _ = airy_vrednost(nb; h=0.001)
            @test abs(val) < 1e-4
        end
    end

    @testset "Newton vs Bisekcija — Newton potrebuje manj iteracij" begin
        nicle, rez_b, rez_n = poisci_nicle(3; h=0.01, x_min=-10.0)
        it_b = [r.iteracije for r in rez_b]
        it_n = [r.iteracije for r in rez_n]
        @test sum(it_n) < sum(it_b)
    end

end