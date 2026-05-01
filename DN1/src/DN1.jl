module DN1

export Zlepek, thomas, interpoliraj, vrednost, plot_zlepek

using Plots

# Tip ki hrani koeficiente kubičnega zlepka
struct Zlepek
    x::Vector{Float64}  # n interpolacijskih točk
    a::Vector{Float64}  # n-1 koeficientov (konstantni členi)
    b::Vector{Float64}  # n-1 koeficientov (linearni členi)
    c::Vector{Float64}  # n koeficientov (kvadratni členi, vključno z obema roboma c[1]=c[n]=0)
    d::Vector{Float64}  # n-1 koeficientov (kubični členi)
end

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
