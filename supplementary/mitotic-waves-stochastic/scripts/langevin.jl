using Printf

include("../src/wave.jl")


_NOISES = 10.0 .^ Vector([-1, -0.7, 0, 0.3, 0.5])
_NGRID = 1024

_N_REPEAT_PER_CONDITION = 10

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

# +noise, -difusion, -source
for n = _NOISES
    nstr = @sprintf "%.3f" n
    cell_cycle_2ode(
        [0.0, 500.0],
        [0.0, 100.0],
        saveat=0.2,
        saveas="data/langevin_-d-s_$(nstr).csv",
        k=k(x) = _K,
        a1=a1(x) = _A1,
        b1=b1(x) = _B1,
        E1=E1(x) = _E1,
        n1=n1(x) = _N1,
        a2=a2(x) = _A2,
        b2=b2(x) = _B2,
        E2=E2(x) = _E2,
        n2=n2(x) = _N2,
        a3=a3(x) = _A3,
        b3=b3(x) = _B3,
        E3=E3(x) = _E3,
        n3=n3(x) = _N3,
        D_u=D_u(x) = 0.0,
        D_v=D_v(x) = 0.0,
        N_u=N_u(x) = n,
        N_v=N_v(x) = n,
        ngrid=_N_REPEAT_PER_CONDITION,
    )
end

# +noise, +difusion, -source
for n = _NOISES
    nstr = @sprintf "%.3f" n
    for i = 1:_N_REPEAT_PER_CONDITION
        istr = @sprintf "%03d" i
        cell_cycle_2ode(
            [0.0, 600.0],
            [0.0, 5000.0],
            saveat=0.2,
            saveas="data/langevin_+d-s_$(nstr)_$(istr).csv",
            k=k(x) = _K,
            a1=a1(x) = _A1,
            b1=b1(x) = _B1,
            E1=E1(x) = _E1,
            n1=n1(x) = _N1,
            a2=a2(x) = _A2,
            b2=b2(x) = _B2,
            E2=E2(x) = _E2,
            n2=n2(x) = _N2,
            a3=a3(x) = _A3,
            b3=b3(x) = _B3,
            E3=E3(x) = _E3,
            n3=n3(x) = _N3,
            D_u=D_u(x) = 240.0,
            D_v=D_v(x) = 240.0,
            N_u=N_u(x) = n,
            N_v=N_v(x) = n,
            ngrid=_NGRID,
        )
    end
end

# +noise, +difusion, +source
k = k(x) = _K * (1 + 0.3 * 0.5 * (sign(25.0 - abs(x - 2500.0)) + 1))

for n = _NOISES
    nstr = @sprintf "%.3f" n
    for i = 1:_N_REPEAT_PER_CONDITION
        istr = @sprintf "%03d" i

        cell_cycle_2ode(
            [0.0, 1440.0],
            [0.0, 5000.0],  # 5 mm
            saveat=0.2,
            saveas="data/langevin_+d+s_$(nstr)_$(istr).csv",
            k=k,
            a1=a1(x) = _A1,
            b1=b1(x) = _B1,
            E1=E1(x) = _E1,
            n1=n1(x) = _N1,
            a2=a2(x) = _A2,
            b2=b2(x) = _B2,
            E2=E2(x) = _E2,
            n2=n2(x) = _N2,
            a3=a3(x) = _A3,
            b3=b3(x) = _B3,
            E3=E3(x) = _E3,
            n3=n3(x) = _N3,
            D_u=D_u(x) = 240.0,
            D_v=D_v(x) = 240.0,
            N_u=N_u(x) = n,
            N_v=N_v(x) = n,
            ngrid=_NGRID,
        )
    end
end
