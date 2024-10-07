A Julia script for stochastic simulation of cell cycle oscillations in 0D and 1D, taken from [our another repo](https://github.com/YangLab-um/mitotic-waves-stochastic) (hash e6637b61c6a52d58ec81ecef394af7ae6924c6f2).
This supplements and cross-checks the simulation by the main codes.

To run the Gillespie's SSA to generate realizations of Cdk1 dynamics time series, use `cell_cycle_2ode` from `src/gillespie.jl`.
The following example generates stochastic trajectories shown in Fig. S6.
```bash
cd mitotic-waves-stochastic
mkdir data
julia --project=. scripts/ssa.jl
```

Mitotic waves can be simulated using `cell_cycle_2ode` from `src/waves.jl`. Refer to the following script file.
```bash
julia --project=. scripts/langevin.jl
```
