using DN3
using Plots
using Printf
using SpecialFunctions

# Znane ničle iz tabel (DLMF / Abramowitz & Stegun)
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

println("="^60)
println("DN3 — Demo: Ničle Airyjeve funkcije")
println("="^60)

# -----------------------------------------------------------------------
# 1. Izračun in izris Airyjeve funkcije
# -----------------------------------------------------------------------
println("\n[1/6] Računam Airyjevo funkcijo z Magnusovo metodo...")

xs, Ai_vals, Aip_vals = airy_resitev(1.5; h=0.005, x_zac=-12.0)

p1 = plot(xs, Ai_vals,
    label="Ai(x)",
    color=:royalblue,
    lw=2,
    xlabel="x",
    ylabel="vrednost",
    title="Airyjeva funkcija Ai(x)",
    legend=:topright,
    ylims=(-0.6, 0.6)
)
plot!(p1, xs, Aip_vals,
    label="Ai'(x)",
    color=:darkorange,
    lw=2,
    ls=:dash
)
hline!(p1, [0.0], color=:gray, lw=0.8, label="")

savefig(p1, "airy_funkcija.png")
println("  → Shranjen: airy_funkcija.png")

# -----------------------------------------------------------------------
# 2. Iskanje ničel
# -----------------------------------------------------------------------
println("\n[2/6] Iščem ničle Airyjeve funkcije...")

nicle, rez_b, rez_n = poisci_nicle(10; h=0.005, x_min=-30.0, h_nicla=0.0005)

println("\n  Prvih 10 ničel Airyjeve funkcije:")
println("  " * "-"^42)
println(@sprintf("  %-4s  %-18s  %-12s", "#", "Ničla", "|Ai(ničla)|"))
println("  " * "-"^42)
for (i, nb) in enumerate(nicle)
    val, _ = airy_vrednost(nb; h=0.0005)
    println(@sprintf("  %-4d  %-18.12f  %-12.2e", i, nb, abs(val)))
end
println("  " * "-"^42)

# -----------------------------------------------------------------------
# 3. Primerjava z znanimi ničlami
# -----------------------------------------------------------------------
println("\n[3/6] Primerjava z znanimi ničlami (DLMF)...")

println("\n  Primerjava z znanimi vrednostmi:")
println("  " * "-"^55)
println(@sprintf("  %-4s  %-18s  %-18s  %-10s", "#", "Znana ničla", "Izračunana", "Napaka"))
println("  " * "-"^55)
for i in eachindex(nicle)
    napaka = abs(nicle[i] - ZNANE_NICLE[i])
    println(@sprintf("  %-4d  %-18.12f  %-18.12f  %-10.2e", i, ZNANE_NICLE[i], nicle[i], napaka))
end
println("  " * "-"^55)

# -----------------------------------------------------------------------
# 4. Primerjava natančnosti za različne korake h
# -----------------------------------------------------------------------
println("\n[4/6] Primerjava natančnosti za h=0.0005 vs h=0.0001...")

nicle_h1, _, _ = poisci_nicle(10; h=0.005, x_min=-30.0, h_nicla=0.0005)
nicle_h2, _, _ = poisci_nicle(10; h=0.005, x_min=-30.0, h_nicla=0.0001)

napake_h1 = [abs(nicle_h1[i] - ZNANE_NICLE[i]) for i in eachindex(ZNANE_NICLE)]
napake_h2 = [abs(nicle_h2[i] - ZNANE_NICLE[i]) for i in eachindex(ZNANE_NICLE)]

println("\n  h = 0.0005 vs h = 0.0001:")
println("  " * "-"^65)
println(@sprintf("  %-4s  %-18s  %-14s  %-14s", "#", "Znana ničla", "Napaka h=0.0005", "Napaka h=0.0001"))
println("  " * "-"^65)
for i in eachindex(ZNANE_NICLE)
    println(@sprintf("  %-4d  %-18.12f  %-14.2e  %-14.2e", i, ZNANE_NICLE[i], napake_h1[i], napake_h2[i]))
end
println("  " * "-"^65)

# Graf primerjave napak za obe vrednosti h
indeksi = 1:length(ZNANE_NICLE)
p4 = plot(indeksi, napake_h1,
    label="h = 0.0005",
    marker=:circle, ms=5,
    color=:royalblue,
    lw=2,
    yscale=:log10,
    xlabel="Ničla #",
    ylabel="Absolutna napaka",
    title="Vpliv koraka h na natančnost ničel",
    legend=:topright,
    minorgrid=true
)
plot!(p4, indeksi, napake_h2,
    label="h = 0.0001",
    marker=:diamond, ms=5,
    color=:crimson,
    lw=2
)
hline!(p4, [1e-10], color=:gray, lw=1, ls=:dash, label="zahtevana natančnost")

savefig(p4, "napake_h.png")
println("  → Shranjen: napake_h.png")

# -----------------------------------------------------------------------
# 5. Graf konvergence bisekcije in Newtonove metode
# -----------------------------------------------------------------------
println("\n[5/6] Rišem graf konvergence...")

a1, b1 = nicle[1] - 0.1, nicle[1] + 0.2
f_ai(x)  = airy_vrednost(x; h=0.0002)[1]
df_ai(x) = airy_vrednost(x; h=0.0002)[2]
ref = nicle[1]

_, _, hist_b = bisekcija(f_ai, a1, b1; tol=1e-12)
_, _, hist_n = newton_raphson(f_ai, df_ai, (a1 + b1)/2; tol=1e-12)

napake_b = max.(abs.(hist_b .- ref), 1e-16)
napake_n = max.(abs.(hist_n .- ref), 1e-16)

p2 = plot(1:length(napake_b), napake_b,
    label="Bisekcija",
    marker=:circle, ms=5,
    color=:royalblue,
    lw=2,
    yscale=:log10,
    xlabel="Iteracija",
    ylabel="Absolutna napaka",
    title="Konvergenca pri 1. ničli Ai(x) ≈ $(round(ref, digits=6))",
    legend=:topright,
    minorgrid=true
)
plot!(p2, 1:length(napake_n), napake_n,
    label="Newton",
    marker=:diamond, ms=5,
    color=:crimson,
    lw=2
)

savefig(p2, "konvergenca.png")
println("  → Shranjen: konvergenca.png")

# -----------------------------------------------------------------------
# 6. Primerjava števila iteracij za vse ničle
# -----------------------------------------------------------------------
println("\n[6/6] Rišem primerjavo iteracij...")

it_b = [r.iteracije for r in rez_b]
it_n = [r.iteracije for r in rez_n]

p3 = bar(indeksi .- 0.2, it_b,
    bar_width=0.35,
    label="Bisekcija",
    color=:royalblue,
    alpha=0.8,
    xlabel="Ničla #",
    ylabel="Število iteracij",
    title="Primerjava iteracij: bisekcija vs. Newton",
    legend=:topright
)
bar!(p3, indeksi .+ 0.2, it_n,
    bar_width=0.35,
    label="Newton",
    color=:crimson,
    alpha=0.8
)

savefig(p3, "primerjava_iteracij.png")
println("  → Shranjen: primerjava_iteracij.png")

# -----------------------------------------------------------------------
# Povzetek
# -----------------------------------------------------------------------
println("\n" * "="^60)
println("POVZETEK")
println("="^60)
println(@sprintf("  Bisekcija — povp. iteracij: %.1f", sum(it_b)/length(it_b)))
println(@sprintf("  Newton    — povp. iteracij: %.1f", sum(it_n)/length(it_n)))
println(@sprintf("  Newton je ~%.1fx hitrejši (po iteracijah)", sum(it_b)/sum(it_n)))
println()
println("  Ničle Airyjeve funkcije (do 10. decimalke):")
for (i, n) in enumerate(nicle)
    println(@sprintf("    a_%d = %.10f", i, n))
end
println("="^60)