# DN1 — Naravni kubični zlepek

Avtor: Sinja Kočica

## Opis naloge

Naloga implementira naravni interpolacijski kubični zlepek za danih `n` interpolacijskih točk `(x_i, f_i)`.

Zlepek `S` je funkcija ki izpolnjuje naslednje pogoje:
1. `S(x_i) = f_i` za vse i — gre skozi vse interpolacijske točke
2. `S` je kubični polinom na vsakem podintervalu `[x_i, x_{i+1}]`
3. `S` je dvakrat zvezno odvedljiva na celotnem interpolacijskem intervalu
4. `S''(x_1) = S''(x_n) = 0` — naravni robni pogoj

## Struktura paketa

```
DN1/
├── demo/
│   └── demo.jl       ← demonstracijska skripta
├── report/
│   └── porocilo.pdf  ← poročilo
├── src/
│   └── DN1.jl        ← glavna koda
├── test/
│   └── runtests.jl   ← avtomatski testi
├── Manifest.toml
├── Project.toml
└── README.md
```

## Uporaba

V terminalu odpremo Julia REPL:
```
julia
```

Nato pa lahko preizkusimo glavne funkcije:
```julia
include("src/DN1.jl")
import .DN1: interpoliraj, vrednost, plot_zlepek

# Definiraj interpolacijske točke
x = [0.0, 1.0, 2.0, 3.0]
f = [0.0, 1.0, 0.0, 1.0]

# Izračunaj zlepek
Z = interpoliraj(x, f)

# Vrednost zlepka v točki
y = vrednost(Z, 1.5)

# Nariši graf
p = plot_zlepek(Z)
display(p)
```

## Poganjanje testov

V Julia REPL-u iz glavne mape repozitorija:

```julia
include("DN1/test/runtests.jl")
```

ali z Julia testnim sistemom:

```julia
] activate DN1
] test
```

Pred prvim zagonom je potrebno dodati še paket Test:
```julia
] activate DN1
] add Test
] test
```

## Poganjanje demo skripte

V Julia REPL-u iz glavne mape repozitorija:

```julia
cd("DN1/demo")
include("demo.jl")
```

Skripta ustvari grafe in jih shrani kot PNG datoteke v mapi `DN1/demo/`.