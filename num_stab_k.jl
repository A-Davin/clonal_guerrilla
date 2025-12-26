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

function main()

    # ============================================================
    # Load model parameters
    # ============================================================

    P = Params()

    # ============================================================
    # Read homogeneous equilibrium data from file
    # ============================================================

    bio = readdlm("unif_low_b.dat")

    b_0 = bio[:,1]    # equilibrium biomass
    w_0 = bio[:,2]    # equilibrium water
    p   = bio[:,3]    # precipitation parameter

    # Cast equilibrium values to complex type
    b_0 .+= 0im
    w_0 .+= 0im

    # ============================================================
    # Definition of the wavenumber k
    # ============================================================

    nk = 256          # number of k values

    k = Array{Complex{Float64}}(undef, nk)

    # Uniform sampling of k
    for j = 1:nk
        k[j] = 0.0001 + ((j-1)/7.5) + 0im
    end

    # ============================================================
    # Selection of precipitation values
    # ============================================================

    ip_list = 4:5:14          # indices of selected p values
    np = length(ip_list)

    # ============================================================
    # Arrays to store eigenvalues
    # ============================================================

    lam1 = Array{ComplexF64}(undef, nk, np)   # first eigenvalue
    lam2 = Array{ComplexF64}(undef, nk, np)   # second eigenvalue

    # ============================================================
    # Main loop over precipitation p and wavenumber k
    # ============================================================

    ip = 0
    for i in ip_list
        ip += 1

        for j = 1:nk

            @printf("p=%1.4f, k=%1.4f \n", p[i+10], k[j])

            # ====================================================
            # Jacobian matrix of the linearized model (2 × 2)
            # ====================================================

            # Vegetation–vegetation interaction term
            a11 =
                (1 + P.η*b_0[i+10])*P.ν*w_0[i+10] *
                (1 - 2*b_0[i+10] + P.η*b_0[i+10]*(3-4*b_0[i+10])) - 1 +
                (0.2+0im)*w_0[i+10]*k[j]*im -
                0.0005*k[j]^2

            # Vegetation–water coupling term
            a12 =
                P.ν*b_0[i+10]*(1-b_0[i+10]) *
                (1 + P.η*b_0[i+10])^2 *
                exp(-(k[j]^2)*(1 + P.η*b_0[i+10])^2/2) +
                (0.2+0im)*b_0[i+10]*k[j]*im

            # Water–vegetation coupling term
            a21 =
                -P.γ*w_0[i+10]*(1 + P.η*b_0[i+10]) *
                exp(-(k[j]^2)*(1 + P.η*b_0[i+10])^2/2) *
                (1 + P.η*b_0[i+10]*
                 (3 - k[j]^2*(1 + P.η*b_0[i+10])^2))

            # Water–water interaction term
            a22 =
                -P.ν -
                P.γ*b_0[i+10]*(1 + P.η*b_0[i+10])^2 -
                P.dw*k[j]^2

            # Jacobian matrix
            J = [
                a11  a12;
                a21  a22
            ]

            # Eigenvalues of the Jacobian
            ei = eigvals(J)

            # Store eigenvalues
            lam1[j,ip] = ei[1]
            lam2[j,ip] = ei[2]

            @printf("eigenvalues: %1.4f %1.4f \n\n",
                    real(ei[1]), real(ei[2]))
        end
    end

    # ============================================================
    # Plot: dispersion relation λ(k) for different p values
    # ============================================================

    pl = plot(xlabel = "wavenumber k",
              ylabel = "Re[λ(k)]",
              legend = :bottomright,
              grid = false)

    ip = 0
    for i in ip_list
        ip += 1
        plot!(real(k),
              real(lam2[:,ip]),
              label = @sprintf("p = %.1f", p[i+10]),
              lw = 3)
    end

    # Zero growth-rate reference line
    pl = hline!([0], lw = 1.5, style = :dash, label = false, c = :black)

    display(pl)

end

# ============================================================
# Script execution
# ============================================================

main()
