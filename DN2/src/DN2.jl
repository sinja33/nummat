module DN2

using LinearAlgebra

export Bezier, tocka, odvod, integriraj_gl, poisci_zanko, ploscina_zanke

"""
    Bezier

Tip ki predstavlja Bézierjevo krivuljo stopnje n, določeno z n+1 kontrolnimi točkami.

Polja:
- `P`: matrika kontrolnih točk velikosti (2, n+1), vsak stolpec je ena točka [x; y]
"""
struct Bezier
    P::Matrix{Float64}
end

"""
    Bezier(points)

Ustvari Bézierjevo krivuljo iz seznama točk podanih kot vektor parov.

# Primer
```julia
pts = [(0,0),(1,1),(2,3),(1,4),(0,4),(-1,3),(0,1),(1,0)]
B = Bezier(pts)
```
"""
function Bezier(points::Vector{<:Tuple})
    P = hcat([[Float64(p[1]), Float64(p[2])] for p in points]...)
    Bezier(P)
end

"""
    tocka(B, t)

Izračuna točko na Bézierjevi krivulji `B` pri parametru `t ∈ [0, 1]`
z de Casteljaujevim algoritmom.

De Casteljaujev algoritem rekurzivno interpolira med kontrolnimi točkami:
na vsakem koraku r zamenjamo n+1-r točk z n-r linearnimi interpolanti.
Po n korakih ostane ena sama točka — vrednost krivulje pri t.

# Primer
```julia
B = Bezier([(0,0),(1,1),(2,0)])
tocka(B, 0.5)
```
"""
function tocka(B::Bezier, t::Float64)
    Q = copy(B.P)
    n = size(Q, 2) - 1
    for r in 1:n
        for k in 1:(n - r + 1)
            Q[:, k] = (1 - t) * Q[:, k] + t * Q[:, k + 1]
        end
    end
    return Q[:, 1]
end

"""
    odvod(B, t)

Izračuna odvod Bézierjeve krivulje `B` pri parametru `t`.

Odvod Bézierjeve krivulje stopnje n je spet Bézierjeva krivulja stopnje n-1,
z novimi kontrolnimi točkami: P'_k = n * (P_{k+1} - P_k).

# Primer
```julia
B = Bezier([(0,0),(1,1),(2,0)])
odvod(B, 0.5)
```
"""
function odvod(B::Bezier, t::Float64)
    n = size(B.P, 2) - 1
    dP = n * diff(B.P, dims=2)
    return tocka(Bezier(dP), t)
end

"""
    integriraj_gl(f, a, b, n_intervalov)

Numerično integrira funkcijo `f` na intervalu `[a, b]` z sestavljeno
Gauss-Legendrovo kvadraturo s 5 točkami na vsakem od `n_intervalov` podintervalov.

Z m točkami je GL kvadratura točna za polinome stopnje do 2m-1.
S 5 točkami je torej točna za polinome stopnje do 9.

# Primer
```julia
integriraj_gl(sin, 0.0, Float64(π), 50)  # ≈ 2.0
```
"""
function integriraj_gl(f, a::Float64, b::Float64, n_intervalov::Int=100)
    # Vozlišča in uteži za 5-točkovno GL kvadraturo na [-1, 1]
    vozlisca = [-0.9061798459386640,
                -0.5384693101056831,
                 0.0,
                 0.5384693101056831,
                 0.9061798459386640]
    utezi    = [ 0.2369268850561891,
                 0.4786286704993665,
                 0.5688888888888889,
                 0.4786286704993665,
                 0.2369268850561891]

    h = (b - a) / n_intervalov
    vsota = 0.0

    for i in 0:(n_intervalov - 1)
        a_i = a + i * h
        b_i = a_i + h
        sredina = (a_i + b_i) / 2
        pol_dolzina = (b_i - a_i) / 2

        vsota_i = 0.0
        for (xi, wi) in zip(vozlisca, utezi)
            vsota_i += wi * f(sredina + pol_dolzina * xi)
        end
        vsota += pol_dolzina * vsota_i
    end

    return vsota
end

"""
    poisci_zanko(B; st_tock=200, tol=1e-10)

Poišče parametra t1 < t2 kjer se Bézierjeva krivulja `B` seče sama s seboj.
Vrne `nothing` če zanka ne obstaja (najbližji par točk je preveč oddaljen).

Postopek:
1. Vzorčimo krivuljo na `st_tock` točkah in poiščemo par ki je si najbližje
2. Če je razdalja večja od praga, zanka ne obstaja → vrnemo `nothing`
3. Z Newtonovo metodo izpopolnimo rešitev sistema B(t1) = B(t2)

Newtonova metoda reši F(t1,t2) = B(t1) - B(t2) = 0,
jakobijan je J = [B'(t1), -B'(t2)].

# Primer
```julia
B = Bezier([(0,0),(1,1),(2,3),(1,4),(0,4),(-1,3),(0,1),(1,0)])
rezultat = poisci_zanko(B)  # vrne (t1, t2)
```
"""
function poisci_zanko(B::Bezier; st_tock::Int=200, tol::Float64=1e-10, prag::Float64=0.5)
    t_vzorci = range(0.0, 1.0, length=st_tock)
    vzorcene_tocke = [tocka(B, t) for t in t_vzorci]

    # Poišči par točk ki sta si najbližje (preskočimo sosednje)
    najboljsi_d = Inf
    najboljsi_i = 0
    najboljsi_j = 0
    for i in 1:st_tock
        for j in (i + 10):st_tock
            d = norm(vzorcene_tocke[i] - vzorcene_tocke[j])
            if d < najboljsi_d
                najboljsi_d = d
                najboljsi_i = i
                najboljsi_j = j
            end
        end
    end

    # Če je najbližji par preveč oddaljen, zanka ne obstaja
    if najboljsi_d > prag
        return nothing
    end

    t1 = Float64(t_vzorci[najboljsi_i])
    t2 = Float64(t_vzorci[najboljsi_j])

    # Newtonova metoda za B(t1) = B(t2)
    for _ in 1:50
        F = tocka(B, t1) - tocka(B, t2)
        if norm(F) < tol
            break
        end
        J = hcat(odvod(B, t1), -odvod(B, t2))
        delta = J \ F
        t1 = clamp(t1 - delta[1], 0.0, 1.0)
        t2 = clamp(t2 - delta[2], 0.0, 1.0)
    end
    # Končna preveritev — če Newton ni konvergiral, vrnemo nothing
    if norm(tocka(B, t1) - tocka(B, t2)) > sqrt(tol)
        return nothing
    end
    if abs(t1 - t2) < 0.05   # t1 in t2 morata biti različna parametra
        return nothing
    end

    return min(t1, t2), max(t1, t2)
end


"""
    ploscina_zanke(B; n_intervalov=100)

Izračuna ploščino zanke ki jo omejuje Bézierjeva krivulja `B`.
Vrne 0.0 če krivulja nima zanke.

Uporablja Greenovo formulo:
    P = (1/2) |∫_{t1}^{t2} (x(t)ẏ(t) - ẋ(t)y(t)) dt|

kjer sta t1, t2 parametra samopreseče krivulje.

# Primer
```julia
pts = [(0,0),(1,1),(2,3),(1,4),(0,4),(-1,3),(0,1),(1,0)]
B = Bezier(pts)
ploscina_zanke(B)  # ≈ 2.2537
```
"""
function ploscina_zanke(B::Bezier; n_intervalov::Int=100)
    rezultat = poisci_zanko(B)
    if rezultat === nothing
        return 0.0
    end
    t1, t2 = rezultat

    integrand = t -> begin
        xy  = tocka(B, t)
        dxy = odvod(B, t)
        xy[1] * dxy[2] - dxy[1] * xy[2]
    end

    return abs(integriraj_gl(integrand, t1, t2, n_intervalov)) / 2
end

end # module