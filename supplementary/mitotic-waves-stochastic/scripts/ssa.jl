using Printf

include("../src/gillespie.jl")


_VOLUMES = 10.0 .^ Vector(-2:0.1:1)

_N_REPEAT_PER_CONDITION = 50

# Synthesis
_K = 1.5

# Cdc25
_A1 = 0.8
_B1 = 4.0
_E1 = 35.0
_N1 = 11.0

# Wee1
_A2 = 0.4
_B2 = 2.0
_E2 = 30.0
_N2 = 3.5

# Degradation
_A3 = 0.01
_B3 = 0.06
_E3 = 32.0
_N3 = 17.0

for v = _VOLUMES
    vstr = @sprintf "%0.3f" v
    for i = 1:_N_REPEAT_PER_CONDITION
        istr = @sprintf "%03d" i

        cell_cycle_2ode(
            [0.0, 1440.0],
            v,
            saveat=0.2,
            saveas="data/ssa_$(vstr)_$(istr).csv",
            k=_K,
            a1=_A1,
            b1=_B1,
            E1=_E1,
            n1=_N1,
            a2=_A2,
            b2=_B2,
            E2=_E2,
            n2=_N2,
            a3=_A3,
            b3=_B3,
            E3=_E3,
            n3=_N3,
        )
    end
end
