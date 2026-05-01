# Numerična matematika — Domače naloge

Repozitorij vsebuje rešitve treh domačih nalog pri predmetu Numerična matematika, implementiranih v programskem jeziku Julia.

## Struktura

Vsaka domača naloga je implementirana kot samostojen Julia paket v svojem direktoriju.

## Naloge

### Domača naloga 1 (poglavje 18.1) — Naravni kubični zlepek

Implementacija naravnega kubičnega zlepka za interpolacijo podatkov. Zlepek je definiran po kosih — na vsakem intervalu med sosednjima točkama je kubični polinom. Pogoji zveznosti vrednosti, prvega in drugega odvoda na stikih ter naravni robni pogoj (S''(x₁) = S''(xₙ) = 0) enolično določijo koeficiente.

Implementirano:
- Tip `Zlepek` ki hrani koeficiente kubičnega zlepka
- Funkcija `interpoliraj(x, f)` ki izračuna koeficiente z Thomasovim algoritmom
- Funkcija `vrednost(Z, x)` ki vrne vrednost zlepka v dani točki
- Funkcija `plot_zlepek(Z)` ki nariše graf zlepka z izmenično rdečo in modro barvo

### Domača naloga 2 (poglavje 18.2)
Naloge s poudarkom na numeričnem računanju vrednosti funkcij na 10 decimalk natančno.

### Domača naloga 3 (poglavje 18.3)
Naloge s poudarkom na numeričnem reševanju diferencialnih enačb na 10 decimalk natančno.

## Zahteve

- Julia 1.x
- Paketi so navedeni v vsakem podprojektu posebej

## Uporaba

Vsako nalogo poženete z:

```julia
include("ime_naloge/src/ime_naloge.jl")
```

ali pa vstopite v direktorij naloge in poženete demonstracijsko skripto.