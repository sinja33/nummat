"""
    DN3 — Ničle Airyjeve funkcije

Paket implementira:
- Magnusovo metodo 4. reda za reševanje ODE y'(x) = A(x)y(x)
- Bisekcijo in Newtonovo metodo za iskanje ničel
- Primerjavo konvergence obeh metod
"""
module DN3

using LinearAlgebra
using Statistics: mean
using Printf
using SpecialFunctions

export magnus_korak, airy_resitev, airy_vrednost
export bisekcija, newton_raphson
export poisci_nicle, primerjaj_metodi

# ---------------------------------------------------------------------------
# Magnusova metoda 4. reda
# ---------------------------------------------------------------------------

"""
    A_matrix(x)

Matrika sistema za Airyjevo ODE. Enačba Ai''(x) = x·Ai(x) se prevede
v sistem prvega reda y' = A(x)y, kjer je y = [Ai(x), Ai'(x)]ᵀ.
"""
function A_matrix(x::Float64)
    return [0.0  1.0;
            x    0.0]
end

"""
    magnus_korak(x, h, y)

En korak Magnusove metode 4. reda za ODE y'(x) = A(x)y(x).

Gaussovi točki znotraj koraka so:
    c₁ = x + (1/2 - √3/6)h
    c₂ = x + (1/2 + √3/6)h

Magnus eksponent:
    σ = h/2·(A₁ + A₂) - √3/12·h²·[A₁, A₂]

kjer je [A, B] = AB - BA komutator matrik.

# Argumenti
- `x`: trenutna točka
- `h`: dolžina koraka
- `y`: trenutni vektor stanja [Ai(x), Ai'(x)]ᵀ

# Vrne
Vektor stanja y v točki x + h.
"""
function magnus_korak(x::Float64, h::Float64, y::Vector{Float64})
    sq3 = sqrt(3.0)
    c1 = x + (0.5 - sq3/6) * h
    c2 = x + (0.5 + sq3/6) * h

    A1 = A_matrix(c1)
    A2 = A_matrix(c2)

    # Komutator [A1, A2] = A1·A2 - A2·A1
    kom = A1 * A2 - A2 * A1

    # Magnus eksponent σ
    sigma = (h/2) * (A1 + A2) - (sq3/12) * h^2 * kom

    # Naslednji korak: y_{k+1} = exp(σ) · y_k
    return exp(sigma) * y
end

"""
    airy_resitev(x_konec; h=0.01, x_zac=-30.0)

Izračuna vrednosti Airyjeve funkcije Ai in njenega odvoda Ai' na
enakomerni mreži od `x_zac` do `x_konec` z Magnusovo metodo 4. reda.

Začetni pogoji pri x=0 (vrednosti so znane eksaktno):
    Ai(0)  = 1 / (2^(1/3) · Γ(2/3))
    Ai'(0) = -1 / (2^(1/3) · Γ(1/3))

Ker integriramo v negativno smer (x_zac < 0), začnemo pri x=0
in hodimo nazaj z negativnim korakom.

# Argumenti
- `x_konec`: končna točka integracije
- `h`: dolžina koraka (privzeto 0.01)
- `x_zac`: začetna točka (privzeto -30.0)

# Vrne
`(xs, Ai_vrednosti, Aip_vrednosti)` — vektori točk in vrednosti.
"""
function airy_resitev(x_konec::Float64=0.0; h::Float64=0.01, x_zac::Float64=-30.0)
    # Začetni pogoji pri x=0 (Abramowitz & Stegun)
    Ai0  =  1.0 / (3^(2/3) * gamma(2/3))
    Aip0 = -1.0 / (3^(1/3) * gamma(1/3))  # Ai'(0) < 0

    # Integriramo nazaj: od 0 do x_zac z negativnim korakom
    n_neg = round(Int, abs(x_zac) / h)
    xs_neg   = zeros(n_neg + 1)
    Ai_neg   = zeros(n_neg + 1)
    Aip_neg  = zeros(n_neg + 1)

    xs_neg[1]  = 0.0
    y = [Ai0, Aip0]
    Ai_neg[1]  = y[1]
    Aip_neg[1] = y[2]

    for i in 2:(n_neg + 1)
        y = magnus_korak(xs_neg[i-1], -h, y)  # negativen korak
        xs_neg[i]  = xs_neg[i-1] - h
        Ai_neg[i]  = y[1]
        Aip_neg[i] = y[2]
    end

    # Obrni, da dobimo naraščajoče x
    reverse!(xs_neg); reverse!(Ai_neg); reverse!(Aip_neg)

    # Integriramo naprej: od 0 do x_konec (za pozitivne x, kjer Ai→0)
    n_pos = round(Int, max(x_konec, 0.0) / h)
    if n_pos > 0
        xs_pos   = zeros(n_pos)
        Ai_pos   = zeros(n_pos)
        Aip_pos  = zeros(n_pos)
        y = [Ai0, Aip0]
        for i in 1:n_pos
            y = magnus_korak((i-1)*h, h, y)
            xs_pos[i]  = i * h
            Ai_pos[i]  = y[1]
            Aip_pos[i] = y[2]
        end
        xs   = vcat(xs_neg, xs_pos)
        Ai   = vcat(Ai_neg, Ai_pos)
        Aip  = vcat(Aip_neg, Aip_pos)
    else
        xs  = xs_neg
        Ai  = Ai_neg
        Aip = Aip_neg
    end

    return xs, Ai, Aip
end

"""
    airy_vrednost(x; h=0.001)

Izračuna vrednost Ai(x) in Ai'(x) v točki x z Magnusovo metodo.
Za negativne x integrira od 0 nazaj, za pozitivne x naprej.

# Vrne
`(Ai_x, Aip_x)` — vrednost funkcije in odvoda v točki x.
"""
function airy_vrednost(x::Float64; h::Float64=0.001)
    Ai0  =  1.0 / (3^(2/3) * gamma(2/3))
    Aip0 = -1.0 / (3^(1/3) * gamma(1/3))

    y = [Ai0, Aip0]
    if x == 0.0
        return y[1], y[2]
    end

    n = round(Int, abs(x) / h)
    hk = x < 0 ? -h : h
    x_curr = 0.0
    for _ in 1:n
        y = magnus_korak(x_curr, hk, y)
        x_curr += hk
    end
    # Morebitni preostanek
    preostanek = x - x_curr
    if abs(preostanek) > 1e-14
        y = magnus_korak(x_curr, preostanek, y)
    end
    return y[1], y[2]
end

# ---------------------------------------------------------------------------
# Bisekcija
# ---------------------------------------------------------------------------

"""
    bisekcija(f, a, b; tol=1e-12, max_iter=100)

Poišče ničlo funkcije `f` na intervalu `[a, b]` z bisekcijo.
Zahteva, da velja f(a)·f(b) < 0.

# Vrne
`(nicla, iteracije, history)` kjer `history` vsebuje zaporedje
polovišč za analizo konvergence.
"""
function bisekcija(f, a::Float64, b::Float64; tol::Float64=1e-12, max_iter::Int=100)
    fa = f(a)
    fb = f(b)
    @assert fa * fb < 0 "Funkcija mora imeti različna predznaka na robovih intervala"

    history = Float64[]
    c = a
    for i in 1:max_iter
        c = (a + b) / 2
        push!(history, c)
        fc = f(c)
        if abs(b - a) / 2 < tol || abs(fc) < tol
            return c, i, history
        end
        if fa * fc < 0
            b = c; fb = fc
        else
            a = c; fa = fc
        end
    end
    return c, max_iter, history
end

# ---------------------------------------------------------------------------
# Newtonova metoda
# ---------------------------------------------------------------------------

"""
    newton_raphson(f, df, x0; tol=1e-12, max_iter=50)

Poišče ničlo funkcije `f` z Newtonovo metodo, začenši pri `x0`.
Potrebuje odvod `df`.

Za Airyjevo funkcijo: df(x) = Ai'(x), ki jo že računamo z Magnusovo metodo.

# Vrne
`(nicla, iteracije, history)` kjer `history` vsebuje zaporedje približkov.
"""
function newton_raphson(f, df, x0::Float64; tol::Float64=1e-12, max_iter::Int=50)
    history = Float64[x0]
    x = x0
    for i in 1:max_iter
        fx  = f(x)
        dfx = df(x)
        if abs(dfx) < 1e-15
            error("Odvod je (skoraj) nič pri x = $x")
        end
        x_new = x - fx / dfx
        push!(history, x_new)
        if abs(x_new - x) < tol
            return x_new, i, history
        end
        x = x_new
    end
    return x, max_iter, history
end

# ---------------------------------------------------------------------------
# Iskanje ničel Airyjeve funkcije
# ---------------------------------------------------------------------------

"""
    poisci_nicle(n_nicel; h=0.005, x_min=-30.0)

Poišče prvih `n_nicel` ničel Airyjeve funkcije na negativni osi.

Strategija:
1. Izračunaj Ai(x) na gosti mreži z Magnusovo metodo.
2. Zaznaj menjave predznaka → okrog vsake je ena ničla.
3. Za vsako ničlo poženi obe metodi (bisekcijo in Newton) in primerjaj.

# Vrne
`(nicle, rezultati_bisekcije, rezultati_newtona)` — vektori ničel
in podrobnosti konvergence vsake metode.
"""
function poisci_nicle(n_nicel::Int; h::Float64=0.005, x_min::Float64=-30.0, h_nicla::Float64=0.0005)
    xs, Ai_vals, Aip_vals = airy_resitev(0.0; h=h, x_zac=x_min)

    intervali = Tuple{Float64,Float64}[]
    for i in 2:length(xs)
        if Ai_vals[i-1] * Ai_vals[i] < 0
            push!(intervali, (xs[i-1], xs[i]))
        end
    end

    intervali = reverse(intervali)
    intervali = intervali[1:min(n_nicel, length(intervali))]

    nicle = Float64[]
    rez_bisekcija = []
    rez_newton    = []

    for (a, b) in intervali
        f_ai(x)  = airy_vrednost(x; h=h_nicla)[1]
        df_ai(x) = airy_vrednost(x; h=h_nicla)[2]

        nicla_b, it_b, hist_b = bisekcija(f_ai, a, b; tol=1e-11)
        x0 = (a + b) / 2
        nicla_n, it_n, hist_n = newton_raphson(f_ai, df_ai, x0; tol=1e-11)

        push!(nicle, nicla_b)
        push!(rez_bisekcija, (nicla=nicla_b, iteracije=it_b, history=hist_b))
        push!(rez_newton,    (nicla=nicla_n, iteracije=it_n, history=hist_n))
    end

    return nicle, rez_bisekcija, rez_newton
end

"""
    primerjaj_metodi(n_nicel=10; h=0.005, x_min=-30.0)

Primerja bisekcijo in Newtonovo metodo za iskanje ničel Airyjeve funkcije.
Izpiše tabelo z rezultati in vrne podatke za vizualizacijo.

# Vrne
`(nicle, bisekcija_iters, newton_iters, bisekcija_napake, newton_napake)`
"""
function primerjaj_metodi(n_nicel::Int=10; h::Float64=0.005, x_min::Float64=-30.0)
    nicle, rez_b, rez_n = poisci_nicle(n_nicel; h=h, x_min=x_min)

    println("\n" * "="^65)
    println("Primerjava bisekcije in Newtonove metode za ničle Ai(x)")
    println("="^65)
    println(@sprintf("%-4s  %-16s  %-8s  %-8s  %-12s  %-12s",
        "#", "Ničla", "It. bis.", "It. New.", "Err bis.", "Err New."))
    println("-"^65)

    bisekcija_iters = Int[]
    newton_iters    = Int[]
    bisekcija_nap   = Float64[]
    newton_nap      = Float64[]

    for i in eachindex(nicle)
        nb  = rez_b[i].nicla
        nn  = rez_n[i].nicla
        itb = rez_b[i].iteracije
        itn = rez_n[i].iteracije
        # Referenca: privzamemo, da je bisekcija točnejša za napako Newtona
        ref = nb
        en  = abs(nn - ref)
        eb  = 0.0  # bisekcija je referenca

        push!(bisekcija_iters, itb)
        push!(newton_iters,    itn)
        push!(bisekcija_nap,   abs(airy_vrednost(nb; h=0.0001)[1]))
        push!(newton_nap,      abs(airy_vrednost(nn; h=0.0001)[1]))

        println(@sprintf("%-4d  %-16.10f  %-8d  %-8d  %-12.2e  %-12.2e",
            i, nb, itb, itn, bisekcija_nap[end], newton_nap[end]))
    end
    println("="^65)
    println(@sprintf("Povprečje iteracij:  bisekcija = %.1f,  Newton = %.1f",
        mean(bisekcija_iters), mean(newton_iters)))
    println()

    return nicle, bisekcija_iters, newton_iters, bisekcija_nap, newton_nap
end

end # module DN3
