# DN2 — Ploščina zanke Bézierjeve krivulje

**Avtor:** Sinja Kočica 

## Opis naloge

Naloga zahteva izračun ploščine zanke, ki jo omejuje Bézierjeva krivulja dana s kontrolnim poligonom:

```
(0, 0), (1, 1), (2, 3), (1, 4), (0, 4), (−1, 3), (0, 1), (1, 0)
```

## Struktura paketa

```
DN2/
├── src/
│   └── DN2.jl          ← implementacija
├── test/
│   └── runtests.jl     ← avtomatski testi
├── demo/
│   └── demo.jl         ← demonstracijska skripta
├── report/
│   └── porocilo.pdf    ← poročilo
│   └── main.tex        ← izvorna koda poročila
├── Project.toml
├── Manifest.toml
└── README.md
```

## Uporaba

V terminalu odpremo Julia REPL:

```
julia
```

Nato aktiviramo paket in preizkusimo glavne funkcije:

```julia
import Pkg; Pkg.activate(".")
using DN2

# Definiraj krivuljo s kontrolnimi točkami
tocke = [(0,0),(1,1),(2,3),(1,4),(0,4),(-1,3),(0,1),(1,0)]
B = Bezier(tocke)

# Vrednost krivulje pri parametru t
tocka(B, 0.5)

# Odvod krivulje pri parametru t
odvod(B, 0.5)

# Poišči presečišče zanke (vrne nothing če zanka ne obstaja)
rezultat = poisci_zanko(B)
if rezultat !== nothing
    t1, t2 = rezultat
end

# Izračunaj ploščino zanke (vrne 0.0 če zanka ne obstaja)
P = ploscina_zanke(B)
```

## Poganjanje testov

V Julia REPL-u iz mape `DN2/`:

```julia
import Pkg; Pkg.activate(".")
include("test/runtests.jl")
```

ali z Julia testnim sistemom:

```julia
] activate .
] test
```

Pred prvim zagonom je potrebno dodati pakete:

```julia
] activate .
] add Test
] add QuadGK
] test
```

## Poganjanje demo skripte

Pred prvim zagonom je potrebno dodati paket Plots:

```julia
] activate .
] add Plots
```

Nato zaženemo demo skripto iz mape `DN2/`:

```julia
import Pkg; Pkg.activate(".")
cd("demo")
include("demo.jl")
```

Skripta izriše grafe za različne primere krivulj in jih shrani kot PNG datoteke v mapi `DN2/demo/`.