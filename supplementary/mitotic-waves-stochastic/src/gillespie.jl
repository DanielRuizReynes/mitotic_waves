using CSV
using Tables

using DifferentialEquations
using Distributions


UNIT_VOLUME = 1e-12  # pL
UNIT_CONCENTRATION = 1e-9  # nM

N_A = 6.022e23  # Avogadro's number


function cell_cycle_2ode(
    tspan,
    vol;
    saveat = [],
    saveas = "",
    k::Number = 1.0,
    a1::Number = 0.16,
    b1::Number = 0.80,
    E1::Number = 35.0,
    n1::Number = 11.0,
    a2::Number = 0.08,
    b2::Number = 0.40,
    E2::Number = 30.0,
    n2::Number = 3.5,
    a3::Number = 0.01,
    b3::Number = 0.04,
    E3::Number = 32.0,
    n3::Number = 17.0,
    init_u::Number = 0.0,
    init_v::Number = 0.0,
)
    tspan = Tuple(tspan)

    function hill(x, E, n)
        # The following form is safe from exponential divergence and possible numerical instability
        return 1 / ((E / x) ^ n + 1)
        # small x Eand instability
    end

    function hillr(x, E, n)
        return 1 / ((x / E) ^ n + 1)
    end

    n_a = N_A * UNIT_VOLUME * UNIT_CONCENTRATION

    r_active_cdk1_synthesis(u, p, t) = k * n_a * vol

    function active_cdk1_synthesis!(int)
        int.u[1] += 1
    end

    r_cdk1_inactivation(u, p, t) = (a2 + b2 * hillr(u[1], E2 * n_a * vol, n2)) * u[1]

    function cdk1_inactivation!(int)
        int.u[1] += -1
        int.u[2] += 1
    end

    r_cdk1_activation(u, p, t) = (a1 + b1 * hill(u[1], E1 * n_a * vol, n1)) * u[2]

    function cdk1_activation!(int)
        int.u[1] += 1
        int.u[2] += -1
    end

    r_active_cdk1_degradation(u, p, t) = (a3 + b3 * hill(u[1], E3 * n_a * vol, n3)) * u[1]

    function active_cdk1_degradation!(int)
        int.u[1] += -1
    end

    r_inactive_cdk1_degradation(u, p, t) = (a3 + b3 * hill(u[1], E3 * n_a * vol, n3)) * u[2]

    function inactive_cdk1_degradation!(int)
        int.u[2] += -1
    end

    inits = [init_u, init_v]

    prob = DiscreteProblem(inits, tspan, nothing)

    jump_prob = JumpProblem(
        prob,
        Direct(),
        ConstantRateJump(r_active_cdk1_synthesis, active_cdk1_synthesis!),
        ConstantRateJump(r_cdk1_inactivation, cdk1_inactivation!),
        ConstantRateJump(r_cdk1_activation, cdk1_activation!),
        ConstantRateJump(r_active_cdk1_degradation, active_cdk1_degradation!),
        ConstantRateJump(r_inactive_cdk1_degradation, inactive_cdk1_degradation!),
    )

    if isempty(saveat)
        @time sol = solve(jump_prob, SSAStepper())

        sol = hcat(sol.t, permutedims(hcat(sol.u...)))
    else
        # Saving interval should be much longer than the inverse of the reaction rates
        if typeof(saveat) <: Number
            saveat = tspan[1]:saveat:tspan[2]
        else
            # saveat given as an array
        end

        sol = zeros(length(saveat), 1 + length(inits))
        sol[1, :] = hcat(0.0, inits...)

        cnt = 2
        for i in 2:size(sol)[1]
            prob = DiscreteProblem(inits, (saveat[i - 1], saveat[i]), nothing)

            jump_prob = JumpProblem(
                prob,
                Direct(),
                ConstantRateJump(r_active_cdk1_synthesis, active_cdk1_synthesis!),
                ConstantRateJump(r_cdk1_inactivation, cdk1_inactivation!),
                ConstantRateJump(r_cdk1_activation, cdk1_activation!),
                ConstantRateJump(r_active_cdk1_degradation, active_cdk1_degradation!),
                ConstantRateJump(r_inactive_cdk1_degradation, inactive_cdk1_degradation!),
            )
            
            sol_ = solve(jump_prob, SSAStepper())

            sol[cnt, :] = hcat(sol_.t[end], sol_.u[end]...)

            inits = sol_.u[end]

            cnt += 1
        end
    end

    if isempty(saveas)
        return sol
    else
        CSV.write(saveas, Tables.table(sol), header=false)
    end
end
