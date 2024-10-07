using Base

using CSV
using Tables

using SciMLBase
using DifferentialEquations
using DomainSets
using MethodOfLines
using ModelingToolkit
using SymbolicUtils


mutable struct Waves1DSolution
    _sol
end


function Base.getindex(s::Waves1DSolution, n::Num)
    return getindex(s._sol, n)
end


mutable struct Waves1DProblem
    _p::SciMLBase.AbstractDEProblem
    _n::Int
    _t
    _x
    _ivs
    _grid
end


function Waves1DProblem(
    eqs::AbstractVector{Equation},
    inits::AbstractVector{Function},
    tspan,
    xspan,
    p::Any = SciMLBase.NullParameters();
    ngrid::Int = 100,
    noise = nothing,
)
    N = ngrid + 1
    # boundaries are excluded

    eqs_ = deepcopy(eqs)
    # Must prevent modification (simplification) of the original equation

    t_ = Num(SymbolicUtils.arguments(SymbolicUtils.arguments(eqs_[1].lhs)[1])[1])
    x_ = Num(SymbolicUtils.arguments(SymbolicUtils.arguments(eqs_[1].lhs)[1])[2])
    
    Dx_ = Differential(x_)

    #get_variables
    vars = []
    bc = []

    for (eq, init) in zip(eqs_, inits)
        u = Num(SymbolicUtils.arguments(eq.lhs)[1])
        push!(vars, u)

        # Initial conditions
        push!(bc, SymbolicUtils.substitute(u, Dict([t_ => tspan[1]])) ~ init(x_))

        # No flux boundary conditions
        push!(bc, Dx_(SymbolicUtils.substitute(u, Dict([x_ => xspan[1]]))) ~ 0.0)
        push!(bc, Dx_(SymbolicUtils.substitute(u, Dict([x_ => xspan[2]]))) ~ 0.0)
    end

    domains = [t_ ∈ Interval(tspan...), x_ ∈ Interval(xspan...)]

    dx = (xspan[2] - xspan[1]) / N
    grid = [i for i in xspan[1]:dx:xspan[2]][2:end - 1]
    discretization = MOLFiniteDifference([x_ => dx], t_)
    # The second arguments (usually time) must be wrapped in a Symbolics.Num variable

    @named pde = PDESystem(eqs_, bc, domains, [t_, x_], vars, p)
    prob = discretize(pde, discretization)

    if isnothing(noise)
        # If noise is not defined, return an ODE problem
        op = ODEProblem(prob.f, prob.u0, prob.tspan, prob.p)

        return Waves1DProblem(op, length(eqs_), t_, x_, vars, grid)
    else
        # If noise is defined
        eq_noise = []
        bc_noise = []

        Dt_ = Differential(t_)
        for (u, n) in zip(vars, noise)
            push!(eq_noise, Dt_(u) ~ n(x_))

            push!(bc_noise, SymbolicUtils.substitute(u, Dict([t_ => tspan[1]])) ~ n(x_))
            push!(bc_noise, SymbolicUtils.substitute(u, Dict([x_ => xspan[1]])) ~ n(xspan[1]))
            push!(bc_noise, SymbolicUtils.substitute(u, Dict([x_ => xspan[2]])) ~ n(xspan[2]))
        end

        @named pdesys_noise = PDESystem(eq_noise, bc_noise, domains, [t_, x_], vars, p)
        prob_noise = discretize(pdesys_noise, discretization)

        sp = SDEProblem(prob.f, prob_noise.f, prob.u0, prob.tspan, prob.p)

        return Waves1DProblem(sp, length(eqs_), t_, x_, vars, grid)
    end
end


function solve(problem::Waves1DProblem, solver; kwargs...)
    if typeof(problem._p) <: SciMLBase.AbstractODEProblem
        # Deterministic
    elseif typeof(problem._p) <: SciMLBase.AbstractSDEProblem
        # Stochastic
    else
        throw("Cannot identify the equation type.")
    end

    sol = DifferentialEquations.solve(problem._p, solver; kwargs...)

    dat = Dict{Num, Any}([problem._t => sol.t, problem._x => problem._grid])

    res = permutedims(hcat(sol.u...))
    chunk_size = Int(size(res)[2] / problem._n)

    @assert chunk_size == size(res)[2] / problem._n

    for (i, u) in enumerate(problem._ivs)
        merge!(dat, Dict([u => res[:, (i - 1) * chunk_size + 1:i * chunk_size]]))
    end

    return Waves1DSolution(dat)
end


function cell_cycle_2ode(
    tspan,
    xspan;
    saveat = [],
    saveas = "",
    k::Function = k(x) = 1.0,
    a1::Function = a1(x) = 0.16,
    b1::Function = b1(x) = 0.80,
    E1::Function = E1(x) = 35.0,
    n1::Function = n1(x) = 11.0,
    a2::Function = a2(x) = 0.08,
    b2::Function = b2(x) = 0.40,
    E2::Function = E2(x) = 30.0,
    n2::Function = n2(x) = 3.5,
    a3::Function = a3(x) = 0.01,
    b3::Function = b3(x) = 0.04,
    E3::Function = E3(x) = 32.0,
    n3::Function = n3(x) = 17.0,
    D_u::Function = D_u(x) = 600.0,
    D_v::Function = D_v(x) = 600.0,
    init_u::Function = init_u(x) = 0.1,
    init_v::Function = init_v(x) = 0.1,
    N_u::Union{Function, Nothing} = nothing,
    N_v::Union{Function, Nothing} = nothing,
    ngrid::Int = 100,
)
    @variables t x u(..) v(..)

    Dt = Differential(t)
    Dxx = Differential(x)^2

    function hill(x, E, n)
        # The following form is safe from exponential divergence and possible numerical instability
        return 1 / (abs(E / x) ^ n + 1)
        # small x Eand instability
    end

    function hillr(x, E, n)
        return 1 / (abs(x / E) ^ n + 1)
    end

    eqs = [
        Dt(u(t, x)) ~ D_u(x) * Dxx(u(t, x)) + k(x) - (a3(x) + b3(x) * hill(v(t, x), E3(x), n3(x))) * u(t, x),
        Dt(v(t, x)) ~ D_v(x) * Dxx(v(t, x)) + k(x) + (a1(x) + b1(x) * hill(v(t, x), E1(x), n1(x))) * (u(t, x) - v(t, x)) - (a2(x) + b2(x) * hillr(v(t, x), E2(x), n2(x))) * v(t, x) - (a3(x) + b3(x) * hill(v(t, x), E3(x), n3(x))) * v(t, x),
    ]

    inits = [init_u, init_v]

    noise = [N_u, N_v]

    # Solve
    if isnothing(N_u) && isnothing(N_v)
        op = Waves1DProblem(eqs, inits, tspan, xspan; ngrid=ngrid)
        sol = solve(op, RK4(); saveat=saveat)
    elseif typeof(N_u) <: Function && typeof(N_v) <: Function
        sp = Waves1DProblem(eqs, inits, tspan, xspan; noise=noise, ngrid=ngrid)
        sol = solve(sp, SOSRA(); saveat=saveat)
    else
        throw("Inconsistent definition of noise.")
    end

    if isempty(saveas)
        return sol[v(t, x)]
    else
        CSV.write(saveas, Tables.table(sol[v(t, x)]), header=false)
    end
end
