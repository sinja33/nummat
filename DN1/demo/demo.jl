using Plots
include("../src/DN1.jl")
import .DN1: Zlepek, interpoliraj, vrednost, plot_zlepek

println("=" ^ 60)
println("Demonstracija naravnega kubičnega zlepka")
println("=" ^ 60)

# -------------------------------------------------------
# 1. Osnovni primer
# -------------------------------------------------------
println("\n1. Osnovni primer")

x = [0.0, 1.0, 2.0, 3.0]
f = [0.0, 1.0, 0.0, 1.0]
Z = interpoliraj(x, f)

println("   Interpolacijske točke: ", x)
println("   Vrednosti:             ", f)
println("   Vrednost zlepka v x=0.5: ", round(vrednost(Z, 0.5), digits=6))
println("   Vrednost zlepka v x=1.5: ", round(vrednost(Z, 1.5), digits=6))
println("   Vrednost zlepka v x=2.7: ", round(vrednost(Z, 2.7), digits=6))

p1 = plot_zlepek(Z)
title!(p1, "Osnovni primer")
savefig(p1, "demo_osnovni.png")
println("   Graf shranjen: demo_osnovni.png")

# -------------------------------------------------------
# 2. Primerjava z sinusom — kako dobra je aproksimacija?
# -------------------------------------------------------
println("\n2. Primerjava zlepka s sinusom")

for n in [5, 10, 20]
    x = collect(range(0.0, 2π, length=n))
    f = sin.(x)
    Z = interpoliraj(x, f)

    x_test = collect(range(0.0, 2π, length=1000))
    napaka = maximum(abs(vrednost(Z, xt) - sin(xt)) for xt in x_test)
    println("   n=$n točk → max napaka = ", round(napaka, sigdigits=3))
end

# Graf primerjave z n=10
x = collect(range(0.0, 2π, length=10))
f = sin.(x)
Z = interpoliraj(x, f)
 
tocka = π / 3
println("   Vrednost zlepka v x=π/3:  ", round(vrednost(Z, tocka), digits=6))
println("   Točna vrednost sin(π/3):  ", round(sin(tocka), digits=6))
 
x_fine = collect(range(0.0, 2π, length=500))
p2 = plot_zlepek(Z)
plot!(p2, x_fine, sin.(x_fine), color=:black, linestyle=:dash, linewidth=2)
title!(p2, "Zlepek vs sin(x), n=10")
savefig(p2, "demo_sinus.png")
println("   Graf shranjen: demo_sinus.png")

# -------------------------------------------------------
# 3. Neenakomerne točke — logaritem
# -------------------------------------------------------
println("\n3. Neenakomerne interpolacijske točke")

x = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
f = log.(x .+ 1.0)
Z = interpoliraj(x, f)

println("   Točke: ", x)
println("   Vrednost zlepka v x=3.0: ", round(vrednost(Z, 3.0), digits=6))
println("   Točna vrednost ln(4.0):  ", round(log(4.0), digits=6))
println("   Vrednost zlepka v x=7.0: ", round(vrednost(Z, 7.0), digits=6))
println("   Točna vrednost ln(8.0):  ", round(log(8.0), digits=6))

p3 = plot_zlepek(Z)
x_fine = collect(range(0.0, 10.0, length=500))
plot!(p3, x_fine, log.(x_fine .+ 1.0), color=:black, linestyle=:dash, label="ln(x+1)")
title!(p3, "Zlepek vs ln(x+1), neenakomerne točke")
savefig(p3, "demo_log.png")
println("   Graf shranjen: demo_log.png")

# -------------------------------------------------------
# 4. Konvergenca — napaka pada z večanjem n
# -------------------------------------------------------
println("\n4. Konvergenca napake")

ns = [4, 6, 8, 10, 15, 20, 30]
napake = Float64[]

for n in ns
    x = collect(range(0.0, 2π, length=n))
    f = sin.(x)
    Z = interpoliraj(x, f)

    x_test = collect(range(0.0, 2π, length=1000))
    napaka = maximum(abs(vrednost(Z, xt) - sin(xt)) for xt in x_test)
    push!(napake, napaka)
    println("   n=$n → max napaka = ", round(napaka, sigdigits=3))
end

p4 = plot(ns, napake, 
    yscale=:log10,
    xlabel="število točk n", 
    ylabel="max napaka (log skala)",
    title="Konvergenca napake zlepka za sin(x)",
    marker=:circle,
    color=:blue,
    legend=false)
savefig(p4, "demo_konvergenca.png")
println("   Graf shranjen: demo_konvergenca.png")

# -------------------------------------------------------
# Združen prikaz vseh grafov
# -------------------------------------------------------
p_vse = plot(p1, p2, p3, p4, layout=4, size=(1000, 800))
savefig(p_vse, "demo_vse.png")
println("\nVsi grafi shranjeni v: demo_vse.png")
println("=" ^ 60)
