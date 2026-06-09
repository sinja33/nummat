# DN3 — Ničle Airyjeve funkcije

**Avtor:** Sinja Kočica

## Opis naloge

Naloga zahteva iskanje ničel Airyjeve funkcije `Ai(x)`, ki je rešitev diferencialne enačbe:

```
Ai''(x) = x · Ai(x)
```

Ničle so poiskane na 10 decimalk natančno z Magnusovo metodo reda 4 za reševanje ODE ter bisekcijo in Newtonovo metodo za iskanje ničel.

## Struktura paketa

```
DN3/
├── src/
│   └── DN3.jl          ← implementacija
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
using DN3

# En korak Magnusove metode
y0 = [0.3550280538878172, -0.2588194037928068]
y1 = magnus_korak(0.0, 0.1, y0)

# Vrednost Ai(x) in Ai'(x) v točki x
Ai_val, Aip_val = airy_vrednost(-2.0)

# Mreža vrednosti Ai(x) na intervalu [-12, 1.5]
xs, Ai_vals, Aip_vals = airy_resitev(1.5; h=0.005, x_zac=-12.0)

# Poišči prvih 10 ničel
nicle, rez_b, rez_n = poisci_nicle(10)

# Primerjaj bisekcijo in Newtonovo metodo
primerjaj_metodi(10)
```

## Poganjanje testov

V Julia REPL-u iz mape `DN3/`:

```julia
] activate .
] test
```

ali:

```julia
import Pkg; Pkg.activate(".")
include("test/runtests.jl")
```

## Poganjanje demo skripte

Iz mape `DN3/`:

```julia
import Pkg; Pkg.activate(".")
cd("demo")
include("demo.jl")
```

Skripta izračuna in izriše:
- Graf Airyjeve funkcije in njenega odvoda
- Tabelo prvih 10 ničel z napakami
- Primerjavo napak za različne korake h
- Graf konvergence bisekcije in Newtonove metode
- Primerjavo števila iteracij

Grafi so shranjeni kot PNG datoteke v mapi `DN3/demo/`.

