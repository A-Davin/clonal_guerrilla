# bwh

Vegetation pattern formation in drylands based on a modified BWH-type model.

This repository contains the numerical codes used in the paper:

“Guerilla Clonal Growth Strategy Leads to Amorphous Pattern Formation in a Drylands Vegetation Model”  
SSRN preprint: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5392061

The model builds upon the vegetation–water framework introduced by Gilad et al. (2004, 2006, 2007), and extends it to investigate the effects of clonal growth strategies on spatial pattern formation.
The scripts for numerical simulations are derived from the repository https://github.com/jhardenberg/bwh/tree/main

---

## Basic usage

Edit `Params.jl` to set the model parameters.  
Unless otherwise specified, default parameter values are consistent with those used in the literature and in the accompanying paper.

### 2D simulations

To run the two–dimensional vegetation–water simulator, execute:

run_example.jl

The script is ready to be run in one of the two pre-configured setups provided in the file.  
It performs time integration of the full 2D model and produces spatial patterns of biomass and soil water.

---

## Linear stability analysis

The scripts used for the linear stability analysis shown in Figure 3 of the paper are also included.

### Dispersion relation λ(k) for fixed precipitation

To compute the dispersion relation as a function of the wavenumber k for selected fixed values of the precipitation parameter p, run:

num_stab_k.jl

This script computes the eigenvalues of the Jacobian matrix of the linearized model around homogeneous steady states and plots the real part of the growth rate λ as a function of k.

### Growth rate λ(p) for fixed wavenumber

To compute the growth rate as a function of precipitation p at a fixed wavenumber k, run:

num_stab_mymod.jl

In this case, the script uses the data contained in:

unif_low_b.dat

This file stores homogeneous equilibrium solutions obtained by numerically iterating the spatially uniform model initialized with very low biomass values. These equilibria are then used as base states for the linear stability analysis.

---

## References

Gilad, E., von Hardenberg, J., Provenzale, A., Shachak, M., & Meron, E. (2007).  
A mathematical model of plants as ecosystem engineers.  
Journal of Theoretical Biology, 244(4), 680–691.  
https://doi.org/10.1016/j.jtbi.2006.08.006

Gilad, E., & von Hardenberg, J. (2006).  
A fast algorithm for convolution integrals with space and time variant kernels.  
Journal of Computational Physics, 216(1), 326–336.  
https://doi.org/10.1016/j.jcp.2005.12.003

Gilad, E., von Hardenberg, J., Provenzale, A., Shachak, M., & Meron, E. (2004).  
Ecosystem Engineers: From Pattern Formation to Habitat Creation.  
Physical Review Letters, 93(9), 098105.  
https://doi.org/10.1103/PhysRevLett.93.098105
