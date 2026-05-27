ENV["GKSwstype"] = "100"
using DN2
using Plots
using LinearAlgebra

# ── Pomožna funkcija za risanje ───────────────────────────────────────────────

function narisi_krivuljo(B, naslov, ime_datoteke)
    rezultat = poisci_zanko(B)

    t_vse = range(0.0, 1.0, length=500)
    xy = [tocka(B, t) for t in t_vse]
    xs = [p[1] for p in xy]
    ys = [p[2] for p in xy]

    cx = [B.P[1, i] for i in 1:size(B.P, 2)]
    cy = [B.P[2, i] for i in 1:size(B.P, 2)]

    if rezultat === nothing
        plt = plot(xs, ys, color=:royalblue, linewidth=2, label="krivulja")
        plot!(plt, cx, cy, color=:lightgray, linewidth=1, linestyle=:dash,
              marker=:circle, markersize=4, label="kontrolni poligon")
        title!(plt, "$naslov — ni zanke")
    else
        t1, t2 = rezultat
        P = ploscina_zanke(B)
        i1 = findfirst(t -> t >= t1, collect(t_vse))
        i2 = findfirst(t -> t >= t2, collect(t_vse))

        plt = plot(xs[1:i1], ys[1:i1], color=:gray, linewidth=2, label="zunaj zanke")
        plot!(plt, xs[i1:i2], ys[i1:i2], color=:royalblue, linewidth=2.5, label="zanka")
        plot!(plt, xs[i2:end], ys[i2:end], color=:gray, linewidth=2, label=false)
        plot!(plt, cx, cy, color=:lightgray, linewidth=1, linestyle=:dash,
              marker=:circle, markersize=4, label="kontrolni poligon")
        px = tocka(B, t1)
        scatter!(plt, [px[1]], [px[2]], color=:red, markersize=8, label="presečišče")
        title!(plt, "$naslov — ploščina = $(round(P, digits=6))")
    end

    xlabel!(plt, "x"); ylabel!(plt, "y")
    plot!(plt, aspect_ratio=:equal, legend=:topright)
    savefig(plt, ime_datoteke)
    println("Graf shranjen: $ime_datoteke")
    return plt
end

# ── 1. EVALVACIJA — različne stopnje krivulj ──────────────────────────────────

println("=" ^ 50)
println("EVALVACIJA KRIVULJ")

# Linearna (stopnja 1)
B_lin = Bezier([(0.0,0.0),(2.0,1.0)])
println("Linearna: tocka(B, 0.5) = $(tocka(B_lin, 0.5))")  # pričakovano: (1, 0.5)

# Kvadratna (stopnja 2)
B_kv = Bezier([(0.0,0.0),(1.0,2.0),(2.0,0.0)])
println("Kvadratna: tocka(B, 0.5) = $(tocka(B_kv, 0.5))")  # pričakovano: (1, 1)
println("Kvadratna: tocka(B, 0.0) = $(tocka(B_kv, 0.0))")  # pričakovano: (0, 0)
println("Kvadratna: tocka(B, 1.0) = $(tocka(B_kv, 1.0))")  # pričakovano: (2, 0)

# Kubična (stopnja 3)
B_kub = Bezier([(0.0,0.0),(0.5,2.0),(1.5,2.0),(2.0,0.0)])
println("Kubična: tocka(B, 0.5) = $(tocka(B_kub, 0.5))")

# Nariši vse tri
t_vse = range(0.0, 1.0, length=300)
plt_eval = plot(title="Evalvacija — linearna, kvadratna, kubična krivulja",
                xlabel="x", ylabel="y", aspect_ratio=:equal)

for (B, ime, barva) in [(B_lin, "linearna", :green),
                         (B_kv,  "kvadratna", :royalblue),
                         (B_kub, "kubična",   :orange)]
    xy = [tocka(B, t) for t in t_vse]
    plot!(plt_eval, [p[1] for p in xy], [p[2] for p in xy],
          color=barva, linewidth=2, label=ime)
    cx = [B.P[1,i] for i in 1:size(B.P,2)]
    cy = [B.P[2,i] for i in 1:size(B.P,2)]
    plot!(plt_eval, cx, cy, color=barva, linewidth=1, linestyle=:dash,
          marker=:circle, markersize=4, label=false)
end
savefig(plt_eval, "bezier_evalvacija.png")
println("Graf shranjen: bezier_evalvacija.png")

# ── 2. ODVOD ──────────────────────────────────────────────────────────────────

println()
println("=" ^ 50)
println("ODVOD")

# Linearna krivulja — odvod mora biti konstanten (2, 1)
for t in [0.0, 0.25, 0.5, 0.75, 1.0]
    println("  odvod(B_lin, $t) = $(odvod(B_lin, t))")  # vedno (2, 1)
end

# Kvadratna — odvod v t=0.5
println("  odvod(B_kv, 0.5) = $(odvod(B_kv, 0.5))")  # pričakovano: (2, 0)

# ── 3. ISKANJE ZANKE ──────────────────────────────────────────────────────────

println()
println("=" ^ 50)
println("ISKANJE ZANKE")

# Krivulja brez zanke
println("Kvadratna krivulja (brez zanke):")
println("  poisci_zanko = $(poisci_zanko(B_kv))")  # nothing

# Krivulja z zanko — figura petica
B_petica = Bezier([(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0),
                   (0.0,0.5),(1.0,0.5),(1.0,0.0),(0.5,-0.5)])
println("Figura petica (z zanko):")
rez = poisci_zanko(B_petica)
println("  poisci_zanko = $rez")
println("  B(t1) = $(tocka(B_petica, rez[1]))")
println("  B(t2) = $(tocka(B_petica, rez[2]))")
println("  razlika = $(norm(tocka(B_petica, rez[1]) - tocka(B_petica, rez[2])))")

narisi_krivuljo(B_kv, "Kvadratna krivulja", "bezier_kvadrat.png")
narisi_krivuljo(B_petica, "Figura petica", "bezier_petica.png")

# ── 4. PLOSCINA — glavna naloga ───────────────────────────────────────────────

println()
println("=" ^ 50)
println("PLOŠČINA ZANKE — krivulja iz naloge")

tocke_naloga = [(0,0),(1,1),(2,3),(1,4),(0,4),(-1,3),(0,1),(1,0)]
B_naloga = Bezier(tocke_naloga)
t1, t2 = poisci_zanko(B_naloga)
println("  Presečišče: t1=$t1, t2=$t2")
println("  Točka presečišča: $(tocka(B_naloga, t1))")
println("  Ploščina zanke: $(ploscina_zanke(B_naloga))")
narisi_krivuljo(B_naloga, "Krivulja iz naloge", "bezier_naloga.png")

# ── 5. KONVERGENČNI TEST ──────────────────────────────────────────────────────

println()
println("=" ^ 50)
println("KONVERGENČNI TEST:")
for n in [10, 50, 100, 500, 1000]
    P_test = ploscina_zanke(B_naloga, n_intervalov=n)
    println("  n=$n  →  P=$P_test")
end
