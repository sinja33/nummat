module DN1

export Zlepek, thomas, interpoliraj, vrednost, plot_zlepek

using Plots

"""
    Zlepek
 
Tip ki hrani koeficiente naravnega kubičnega zlepka.
 
Na vsakem intervalu [x[i], x[i+1]] je zlepek definiran kot kubični polinom:
 
    S_i(x) = a[i] + b[i]*(x - x[i]) + c[i]*(x - x[i])^2 + d[i]*(x - x[i])^3
 
# Polja
- `x::Vector{Float64}`: interpolacijske točke, dolžine n
- `a::Vector{Float64}`: konstantni členi, dolžine n-1
- `b::Vector{Float64}`: linearni členi, dolžine n-1
- `c::Vector{Float64}`: kvadratni členi, dolžine n (vključno z robovima c[1]=c[n]=0)
- `d::Vector{Float64}`: kubični členi, dolžine n-1
"""
struct Zlepek
    x::Vector{Float64}  # n interpolacijskih točk
    a::Vector{Float64}  # n-1 koeficientov (konstantni členi)
    b::Vector{Float64}  # n-1 koeficientov (linearni členi)
    c::Vector{Float64}  # n koeficientov (kvadratni členi, vključno z obema roboma c[1]=c[n]=0)
    d::Vector{Float64}  # n-1 koeficientov (kubični členi)
end


"""
    thomas(l, d, u, b)
 
Reši tridiagonalni sistem linearnih enačb z Thomasovim algoritmom.
 
Sistem je oblike A*x = b kjer je A tridiagonalna matrika z:
- glavno diagonalo `d`
- spodnjo diagonalo `l`
- zgornjo diagonalo `u`
 
Algoritem deluje v dveh korakih:
1. Naprej (forward sweep): eliminacija spodnje diagonale
2. Nazaj (back substitution): reševanje za neznanke
 
Časovna zahtevnost je O(n), prostorska O(n).
 
# Argumenti
- `l::Vector{Float64}`: spodnja diagonala, dolžine n-1
- `d::Vector{Float64}`: glavna diagonala, dolžine n
- `u::Vector{Float64}`: zgornja diagonala, dolžine n-1
- `b::Vector{Float64}`: desna stran sistema, dolžine n
 
# Rezultat
Vrne vektor rešitev `x` dolžine n.
 
# Primer
```julia
# Reši sistem s tridiagonalno matriko
l = [1.0, 1.0]
d = [4.0, 4.0, 4.0]
u = [1.0, 1.0]
b = [1.0, 2.0, 3.0]
x = thomas(l, d, u, b)
```
"""
function thomas(l::Vector{Float64}, d::Vector{Float64}, 
                u::Vector{Float64}, b::Vector{Float64})
    m = length(d)
    
    # Kopije da ne spreminjamo originalnih vektorjev
    d = copy(d)
    b = copy(b)
    
    # Forward sweep
    for i in 2:m
        faktor = l[i-1] / d[i-1]
        d[i] -= faktor * u[i-1]
        b[i] -= faktor * b[i-1]
    end
    
    # Back substitution
    x = zeros(m)
    x[m] = b[m] / d[m]
    for i in m-1:-1:1
        x[i] = (b[i] - u[i] * x[i+1]) / d[i]
    end
    
    return x
end


"""
    interpoliraj(x, f)
 
Izračuna koeficiente naravnega kubičnega zlepka za podane interpolacijske točke.
 
Zlepek izpolnjuje naslednje pogoje:
1. Interpolacija: S(x[i]) = f[i] za vse i
2. Zveznost vrednosti, prvega in drugega odvoda na stikih
3. Naravni robni pogoj: S''(x[1]) = S''(x[n]) = 0
 
Koeficienti se izračunajo z reševanjem tridiagonalnega sistema za kvadratne
člene c[i], iz katerih sledijo linearni b[i] in kubični d[i] členi.
Za reševanje tridiagonalnega sistema se uporabi Thomasov algoritem z
časovno zahtevnostjo O(n).
 
# Argumenti
- `x::Vector{Float64}`: interpolacijske točke, morajo biti urejene naraščajoče
- `f::Vector{Float64}`: vrednosti funkcije v interpolacijskih točkah
 
# Rezultat
Vrne element tipa `Zlepek` s koeficienti kubičnih polinomov na vsakem intervalu.
 
# Primer
```julia
x = [0.0, 1.0, 2.0, 3.0]
f = [0.0, 1.0, 0.0, 1.0]
Z = interpoliraj(x, f)
```
"""
function interpoliraj(x::Vector{Float64}, f::Vector{Float64})
    n = length(x)
    h = diff(x)  # h[i] = x[i+1] - x[i]

    m = n - 2

    # Stara implementacija (A \ b):
    # A = zeros(m, m)
    # b = zeros(m)
    # for i in 1:m
    #     b[i] = 3 * ((f[i+2] - f[i+1]) / h[i+1] - (f[i+1] - f[i]) / h[i])
    #     A[i, i] = 2 * (h[i] + h[i+1])
    #     if i > 1; A[i, i-1] = h[i] end
    #     if i < m; A[i, i+1] = h[i+1] end
    # end
    # c_inner = A \ b

    # Thomasov algoritem ker imamo tridiagonalno strukturo, O(n) namesto O(n^3)
    dl = [h[i] for i in 2:m]
    dd = [2 * (h[i] + h[i+1]) for i in 1:m]
    du = [h[i+1] for i in 1:m-1]
    bb = [3 * ((f[i+2] - f[i+1]) / h[i+1] - (f[i+1] - f[i]) / h[i]) for i in 1:m]

    c_inner = thomas(dl, dd, du, bb)

    # dodana robna pogoja
    c = [0.0; c_inner; 0.0]

    # koeficienti a, b, d
    a = f[1:n-1]
    bv = [(f[i+1] - f[i]) / h[i] - h[i] / 3 * (2 * c[i] + c[i+1]) for i in 1:n-1]
    d = [(c[i+1] - c[i]) / (3 * h[i]) for i in 1:n-1]

    return Zlepek(x, a, bv, c, d)
end


"""
    vrednost(Z, x)
 
Izračuna vrednost zlepka v točki x.
 
Najprej poišče interval [x[i], x[i+1]] v katerem leži x, nato
izračuna vrednost kubičnega polinoma na tem intervalu:
 
    S_i(x) = a[i] + b[i]*dx + c[i]*dx^2 + d[i]*dx^3
 
kjer je dx = x - x[i].
 
# Argumenti
- `Z::Zlepek`: zlepek izračunan z `interpoliraj`
- `x::Float64`: točka v kateri računamo vrednost
 
# Rezultat
Vrne vrednost zlepka v točki x.
 
# Primer
```julia
x = [0.0, 1.0, 2.0, 3.0]
f = [0.0, 1.0, 0.0, 1.0]
Z = interpoliraj(x, f)
y = vrednost(Z, 1.5)  # vrednost zlepka pri x=1.5
```
"""
function vrednost(Z::Zlepek, x::Float64)
    # iskanje pravega intervala
    i = searchsortedlast(Z.x, x)
    
    # robni primeri
    if i == 0
        i = 1
    elseif i == length(Z.x)
        i = length(Z.x) - 1
    end
    
    # racunanje vrednosti zlepka
    dx = x - Z.x[i]
    return Z.a[i] + Z.b[i]*dx + Z.c[i]*dx^2 + Z.d[i]*dx^3
end


"""
    plot_zlepek(Z)
 
Nariše graf zlepka z izmenično rdečo in modro barvo za posamezne odseke.
 
Interpolacijske točke so označene s črnimi pikami.
 
# Argumenti
- `Z::Zlepek`: zlepek izračunan z `interpoliraj`
 
# Rezultat
Vrne objekt tipa `Plots.Plot`.
 
# Primer
```julia
x = [0.0, 1.0, 2.0, 3.0]
f = [0.0, 1.0, 0.0, 1.0]
Z = interpoliraj(x, f)
p = plot_zlepek(Z)
display(p)
```
"""
function plot_zlepek(Z::Zlepek)
    colors = [:red, :blue]
    p = Plots.plot(legend=false)
    
    for i in 1:(length(Z.x)-1)
        xs = range(Z.x[i], Z.x[i+1], length=100)
        ys = [vrednost(Z, x) for x in xs]
        Plots.plot!(p, xs, ys, color=colors[mod1(i, 2)])
    end
    
    Plots.scatter!(p, Z.x, [vrednost(Z, x) for x in Z.x], color=:black)
    
    return p
end

end # module DN1
