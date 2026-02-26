
# PRE
- Simulations (S): Lindblad, DM or Thermalize, KrylovSpectrum, Trajectory
  Note, KrylovSpectrum is replacing Lindblad in the main experiments for the thesis, but we need Lindblad still for crosschecking Krylov, since Lindblad is meant to be the exact variant, while Krylov brings in imprecision.
- Domains (D): BohrDomain, EnergyDomain, TimeDomain, TrotterDomain
- Configs {S, D}: Parametrized by Simulation type and Domain type. And can be varying dependent on these, but also has main variants: 
	  - __KMS__ - exact KMS DB with the coherent part B and corresponding $\gamma, \alpha, b_1, b_2$ etc.
	  - __GNS__ - approx GNS DB without the coherent part and has KMS condition fulfilling $\gamma$ and thus different $\alpha$ as well.
	  - KMS and GNS has with_linear_combination() where we always choose the right functions in precompute_data based on the config and this field being true or false, and having parameters, $a, b, \eta$
	  - A milestone soon will be adding: __DLL__ type config, which follows Ding et al construction, without energy integral, with different filter functions and Kossakowski matrices. Also has coherent term of the same principle. Has no liner combination variant and needs no energy labels at all.
- Usual jump set for all of them - keep the logic that for TrotterDomain, jump.in_eigenbasis is written in Trotter basis, and so is all objects along the simulations are in trotter basis so we never have to change from one basis to another in the loops. In the other domain cases they are in the H eigenbasis.

# MID
- For each simulation I want a main script under simulations/. That sets up everything via PRE that I can play with to check things out. Plus this is gonna be basically the scripts that would run on the cluster too for large scale simulations.
- Thus, I want also 4 run_* functions for each simulation (written for KMS and GNS cases atm, but later we will figure out how best to add the DLL variant.)
	- run_lindblad() - at this point for testing and crosschecks to have an exact Lindblad spectrum result even if not with full eigen()
	- run_thermalize() - runs the DM simulator, faithful to Chen's quantum algorithm, CPTP etc, sampling from the jumps $A^a$ every $\delta$ step too.
	- run_krylov_spectrum() - estimates the important part of the spectrum for both the exact Lindbladian or the CPTP one thats faithful to Chen's quantum algorithm again. This is applying the full Lindbladians, the one that the jump sampling from $A^a$ is averaging too for enough sampling.
	- run_trajectory() - Quantum Trajectory simulation, either with or without an observable set. If without  we are mainly interested in the DM and possibly can average to a DM every some steps to see how it's evolving. Or with observables where we don't want to reconstruct averaged DM's midrun but only track how the observable values are approaching the thermal average, and construct only the end averaged DM from the trajectories. We have started writing up a spectral gap estimator with trajectories, but it is unfinished (see unfinished phases) and we might not need it for the thesis, so that can be put on ice for a while. But generally a quantum trajectory simulator for observable convergence is useful. Only the spectral gap estimation seems to be difficult with it. 
# Post
I dont think we need a unifying Results struct, we can just keep 4 different ones for 4 different simulations.
- LindbladResults
- ThermalizeResults
- KrylovSpectrumResults
- TrajectoryResults
each with a boolean field to save the results struct under results/ or not. (We need a function for writing data that takes in either of the 4 results structs and writes a name descriptive name accordingly.)

General rules:
- Most important rule: keep the code as minimal as possible, don't have duplicates. If one function is very similar to another like currently krylov_oft and oft! then keep choose the one that's more efficient and use that for both cases. In this example I would just keep oft!() name, but maybe with the slightly more efficient logic from krylov_oft. 
  Of course, functions should generally be kept to more or less one function, but I really want you eliminate unnecessary extra functions and layers. The code feels bloated to me right now and I keep finding functions written in the past phases where you just created a new function instead of using an already existing one. Now the core logic is more or less there, so you can really think through how to best slim down the codebase.
- One tip right now specifically: It seems that Krylov functions like apply_delta_channel!() has a lot of overlaps with the DM thermalization simulator. Only that DM thermalization needs R^a, B^a, U_residual^a  per delta steps while Krylov can use the summed up ones R_total, B_total, U_residual_total, etc. Which means there are quite some overlaps in logics between the two simulators, and I think the inner functions can definitely be unified such that both simulators can call them but just slightly differently.
- This slimming down is also true for test functions. You keep creating new setups and specialized functions for create this and that, small trotter, small hamiltonian etc. Thats just unnecessary duplications... I want you to go through the test functions as well, and try to slim down the helper functions, eliminate duplications and partial duplications as much as possible. I also want tests to write out an "info" so I can see in the printout what has been really tested and what passed. As you go through the tests also check how much sense they make, always worth to second check, because I felt like some of your testing parameters and values, thresholds have been off time to time and made the tests sometimes more deceiving than useful.
- The necessary workspaces that we need for all the preallocations and infrastructure can be kept internal. I think there you kept attention how to make if most efficient, but maybe rethink again if anything can be elminiated because of overlaps in between workspaces, structs and try to reduce the number of them. Of course up to reason. I think these workspaces could be also homogenized in naming over the 4 simulations, because I guess each simulation needs a workspace. Though some domains need different things preallocated like Time and TrotterDomains always need NUFFT prefactor matrices, etc.


Some future plans still that might influence how you refactor the code now:
- We will add a DLL config that will use different filter functions, no energy labels, has BohrDomain, TimeDomain, TrotterDomain versions. Results of it will be directly compared to results of Chen's KMS DB variant.
- We will add a Hamiltonian simulation time accumulator based on config settings.
- We will add qiskit quantum circuit for the subroutines in the quantum algorithms, like OFT, quantum signal processing for the coherent evolution via B, etc. So we can compute the amount of gates that are needed for these quantum algorithms to converge to the Gibbs state. (NO quantum circuit simulations, but just gate estimations.)
- Compute quantum discriminant of the Lindbladians, compute the norm of their antihermitian part to tell what setups actually lead to less antihermitian part and thus be better approximating Gibbs state as fixed point.
- Quite separately but we will create a plot also the support of the jump operators via Lieb Robinson bounds.
- Also look at Kossakowski matrices in the different cases, i.e. the BohrDomain coefficient matrix $\alpha$
- Add an Ising Hamiltonian in 1D option, and also possibly a 2D Heisenberg one.
- 