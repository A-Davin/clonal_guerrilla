using Parameters
using Statistics
using DiffEqOperators, SparseArrays, LinearAlgebra
using BenchmarkTools
using Printf
using DelimitedFiles
using Plots
#using PlotlyJS

include("Params.jl")

export Params

# ============================================================
# Load model parameters
# ============================================================

P = Params()

# ============================================================
# Read homogeneous equilibrium data from file
# ============================================================

bio = readdlm("unif_low_b.dat")

b_0 = bio[:,1]      # equilibrium biomass
w_0 = bio[:,2]      # equilibrium water
p   = bio[:,3]      # precipitation parameter

# Cast equilibrium values to complex type
b_0 .+= 0im
w_0 .+= 0im

# ============================================================
# Definition of the wavenumber k
# ============================================================

k = Array{Complex{Float64}}(undef, 256)

zr = zeros(256)     # zero reference line for plotting

# Arrays to store eigenvalues
lam1 = Array{Complex{Float64}}(undef, 256)
lam2 = Array{Complex{Float64}}(undef, 256)
lam3 = Array{Complex{Float64}}(undef, 256)

# Uniform sampling of k
for j = 1:256
    k[j] = 0.0001 + ((j-1)/7.5) + 0im
end

# ============================================================
# Loop over precipitation values
# ============================================================

for i = 1:19

    # ========================================================
    # Fixed wavenumber index (single-k analysis)
    # ========================================================

    for j = 6:6

        @printf("p=%1.4f, k=%1.4f \n", p[i+10], k[j])
        # +(0.2+0im)*b_0[i]*k[j]*im

        # ====================================================
        # Jacobian matrix of the linearized model (2 × 2)
        # ====================================================

        a11 = (1 + P.η*b_0[i+10])*P.ν*w_0[i+10] *
	      (1 - 2*b_0[i+10] + P.η*b_0[i+10]*(3-4*b_0[i+10])) - 1 +
  	      (0.2+0im)*w_0[i+10]*k[j]*im -
      	      0.0005*k[j]^2

	a12 = P.ν*b_0[i+10]*(1-b_0[i+10]) *
	      ((1 + P.η*b_0[i+10])^2) *
	      exp(-(k[j]^2)*((1 + P.η*b_0[i+10])^2)/2) +
	      (0.2+0im)*b_0[i+10]*k[j]*im

	a21 = -P.γ*w_0[i+10]*(1 + P.η*b_0[i+10]) *
	      exp(-(k[j]^2)*((1 + P.η*b_0[i+10])^2)/2) *
	      (1 + P.η*b_0[i+10]*
	      (3 - (k[j]^2)*((1 + P.η*b_0[i+10])^2)))

	a22 = -P.ν -
	      P.γ*b_0[i+10]*((1 + P.η*b_0[i+10])^2) -
 	      P.dw*(k[j]^2)

	J = [
	    a11  a12;
	    a21  a22
	]

        # ====================================================
        # Eigenvalues of the Jacobian
        # ====================================================

        ei = eigvals(J)

        # Store the second eigenvalue as a function of p
        lam2[i] = ei[2]

        @printf("eigenvalues: \n %1.4f %1.4f \n\n",
                real(ei[1]), real(ei[2]))

    end
end

# ============================================================
# Plot: growth rate as a function of precipitation p
# (fixed wavenumber k)
# ============================================================

labels = ["p=2.5" "p=5.0" "p=7.5"]

pl = plot(p[11:29],
          real(lam2[1:19]),
          lw = 4,
          label = false,
          xlabel = "precipitation rate",
          ylabel = "Re[λ(k = 0.7, p)]",
          grid = false,
          legend = :bottomright,
          linecolor = :red)

# Zero growth-rate reference line
pl = plot!(p[5:29],
           zr[5:29],
           style = :dash,
           lw = 2,
           label = false,
           linecolor = :black)

display(pl)
