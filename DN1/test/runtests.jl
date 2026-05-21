using Test
using Plots
include("../src/DN1.jl")
import .DN1: Zlepek, interpoliraj, vrednost, plot_zlepek

# -------------------------------------------------------
# Pomožna funkcija: preveri da zlepek gre skozi vse točke
# -------------------------------------------------------
function preveri_interpolacijo(Z, x, f; atol=1e-10)
    for i in 1:length(x)
        @test vrednost(Z, x[i]) ≈ f[i] atol=atol
    end
end

# -------------------------------------------------------
# 1. Osnovni test
# -------------------------------------------------------
@testset "Osnovni test" begin
    x = [0.0, 1.0, 2.0, 3.0]
    f = [0.0, 1.0, 0.0, 1.0]
    Z = interpoliraj(x, f)

    preveri_interpolacijo(Z, x, f)

    # Naravni robni pogoj: c[1] = c[n] = 0
    @test Z.c[1] ≈ 0.0
    @test Z.c[end] ≈ 0.0
end

# -------------------------------------------------------
# 2. Linearna funkcija — zlepek mora biti točen
#    (kubični zlepek je vsaj tako dober kot linearni)
# -------------------------------------------------------
@testset "Linearna funkcija" begin
    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    f = 2.0 .* x .+ 1.0  # f(x) = 2x + 1
    Z = interpoliraj(x, f)

    preveri_interpolacijo(Z, x, f)

    # Vrednosti med točkami morajo biti tudi točne
    @test vrednost(Z, 0.5) ≈ 2.0
    @test vrednost(Z, 1.5) ≈ 4.0
    @test vrednost(Z, 3.7) ≈ 8.4
end

# -------------------------------------------------------
# 3. Kvadratna funkcija — zlepek mora biti točen v interpolacijskih točkah
# -------------------------------------------------------
@testset "Kvadratna funkcija" begin
    x = collect(range(0.0, 4.0, length=5))
    f = x .^ 2
    Z = interpoliraj(x, f)

    # Točen je samo v interpolacijskih točkah
    preveri_interpolacijo(Z, x, f)

    # Naravni robni pogoj
    @test Z.c[1] ≈ 0.0
    @test Z.c[end] ≈ 0.0
end

# -------------------------------------------------------
# 4. Sinus z malo točkami
# -------------------------------------------------------
@testset "Sinus (6 točk)" begin
    x = collect(range(0.0, 2π, length=6))
    f = sin.(x)
    Z = interpoliraj(x, f)

    preveri_interpolacijo(Z, x, f)
end

# -------------------------------------------------------
# 5. Sinus z veliko točkami — napaka mora biti zelo majhna
# -------------------------------------------------------
@testset "Sinus (20 točk) — natančnost" begin
    x = collect(range(0.0, 2π, length=20))
    f = sin.(x)
    Z = interpoliraj(x, f)

    preveri_interpolacijo(Z, x, f)

    # Preverimo natančnost med točkami
    x_test = collect(range(0.0, 2π, length=100))
    for xt in x_test
        @test vrednost(Z, xt) ≈ sin(xt) atol=1e-4
    end
end

# -------------------------------------------------------
# 6. Eksponentna funkcija
# -------------------------------------------------------
@testset "Eksponentna funkcija" begin
    x = collect(range(0.0, 3.0, length=8))
    f = exp.(x)
    Z = interpoliraj(x, f)

    preveri_interpolacijo(Z, x, f)

    # Aproksimacija med točkami — z 8 točkami pričakujemo napako ~1e-3
    @test vrednost(Z, 1.5) ≈ exp(1.5) atol=1e-2
end

# -------------------------------------------------------
# 7. Neenakomerno razporejene točke
# -------------------------------------------------------
@testset "Neenakomerne točke" begin
    x = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    f = log.(x .+ 1.0)  # f(x) = ln(x+1)
    Z = interpoliraj(x, f)

    preveri_interpolacijo(Z, x, f)

    @test Z.c[1] ≈ 0.0
    @test Z.c[end] ≈ 0.0
end

# -------------------------------------------------------
# 8. Samo 3 točke (minimalen primer)
# -------------------------------------------------------
@testset "Minimalni primer (3 točke)" begin
    x = [0.0, 1.0, 2.0]
    f = [1.0, 3.0, 2.0]
    Z = interpoliraj(x, f)

    preveri_interpolacijo(Z, x, f)

    @test Z.c[1] ≈ 0.0
    @test Z.c[end] ≈ 0.0
end

# -------------------------------------------------------
# 9. Konstantna funkcija
# -------------------------------------------------------
@testset "Konstantna funkcija" begin
    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    f = fill(5.0, 5)
    Z = interpoliraj(x, f)

    preveri_interpolacijo(Z, x, f)

    # Vrednosti med točkami morajo biti tudi 5.0
    @test vrednost(Z, 0.5) ≈ 5.0 atol=1e-10
    @test vrednost(Z, 2.3) ≈ 5.0 atol=1e-10
end

# -------------------------------------------------------
# 10. Veliko točk — robustnost
# -------------------------------------------------------
@testset "Veliko točk (50)" begin
    x = collect(range(0.0, 2π, length=50))
    f = sin.(x) .+ 0.5 .* cos.(2 .* x)
    Z = interpoliraj(x, f)

    preveri_interpolacijo(Z, x, f)

    @test Z.c[1] ≈ 0.0
    @test Z.c[end] ≈ 0.0
end

# -------------------------------------------------------
# Vizualizacija
# -------------------------------------------------------
plots = []

# Osnovni primer
x1 = [0.0, 1.0, 2.0, 3.0]
f1 = [0.0, 1.0, 0.0, 1.0]
push!(plots, plot_zlepek(interpoliraj(x1, f1)))

# Sinus
x2 = collect(range(0.0, 2π, length=10))
f2 = sin.(x2)
push!(plots, plot_zlepek(interpoliraj(x2, f2)))

# Eksponentna
x3 = collect(range(0.0, 2.0, length=8))
f3 = exp.(x3)
push!(plots, plot_zlepek(interpoliraj(x3, f3)))

# Logaritem z neenakomernimi točkami
x4 = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
f4 = log.(x4 .+ 1.0)
push!(plots, plot_zlepek(interpoliraj(x4, f4)))

display(Plots.plot(plots..., layout=4, size=(900, 600)))