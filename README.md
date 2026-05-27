# Numerična matematika — Domače naloge

Repozitorij vsebuje rešitve treh domačih nalog pri predmetu Numerična matematika, implementiranih v programskem jeziku Julia.

## Struktura

Vsaka domača naloga je implementirana kot samostojen Julia paket v svojem direktoriju.

## Naloge

### Domača naloga 1 (poglavje 18.1) — Naravni kubični zlepek

Implementacija naravnega kubičnega zlepka za interpolacijo podatkov. Zlepek je definiran po kosih — na vsakem intervalu med sosednjima točkama je kubični polinom. Pogoji zveznosti vrednosti, prvega in drugega odvoda na stikih ter naravni robni pogoj ($S''(x_1) = S''(x_n) = 0$) enolično določijo koeficiente.

Implementirano:
- Tip `Zlepek` ki hrani koeficiente kubičnega zlepka
- Funkcija `interpoliraj(x, f)` ki izračuna koeficiente z Thomasovim algoritmom
- Funkcija `vrednost(Z, x)` ki vrne vrednost zlepka v dani točki
- Funkcija `plot_zlepek(Z)` ki nariše graf zlepka z izmenično rdečo in modro barvo

### Domača naloga 2 (poglavje 18.2) — Ploščina zanke Bézierjeve krivulje

Izračun ploščine zanke, ki jo omejuje Bézierjeva krivulja dana s kontrolnim poligonom osmih točk. Krivulja se seče sama s seboj in tvori zanko, katere ploščino izračunamo z Greenovo formulo.

Implementirano:
- Tip `Bezier` ki hrani kontrolne točke krivulje
- Funkcija `tocka(B, t)` ki izračuna točko na krivulji z de Casteljaujevim algoritmom
- Funkcija `odvod(B, t)` ki izračuna odvod krivulje analitično
- Funkcija `integriraj_gl(f, a, b, n)` ki numerično integrira z Gauss-Legendrovo kvadraturo
- Funkcija `poisci_zanko(B)` ki poišče presečišče krivulje z Newtonovo metodo
- Funkcija `ploscina_zanke(B)` ki izračuna ploščino zanke z Greenovo formulo

### Domača naloga 3 (poglavje 18.3)

TBD

## Uporaba

Vsako nalogo aktivirate in poženete v Julia REPL-u iz direktorija naloge:

```julia
import Pkg; Pkg.activate(".")
using ImeNaloge
```

ali pa poženete demonstracijsko skripto:

```julia
import Pkg; Pkg.activate(".")
cd("demo")
include("demo.jl")
```

## Poganjanje testov

Teste za posamezno nalogo poženete v Julia REPL-u iz direktorija naloge:

```julia
import Pkg; Pkg.activate(".")
include("test/runtests.jl")
```

ali z Julia testnim sistemom:

```julia
] activate .
] test
```