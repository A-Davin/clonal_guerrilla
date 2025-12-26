using Parameters
using Statistics
using Plots
using ImageFiltering
using DiffEqOperators, SparseArrays, LinearAlgebra
using DifferentialEquations
using FFTW
using BenchmarkTools
using Printf
using DelimitedFiles
using OrdinaryDiffEq

include("Params.jl")
include("laplacian.jl")
include("convapprox.jl")
include("rhs.jl")
include("integrate.jl")
include("plotbwh.jl")
include("main.jl")

export main, Params, integrate, plotbwh

P = Params(nx=64, ny=64, Lx=28, Ly=28, p=2, fplot=true, dt=0.1, nstep=200) #small setup for quick trials
#P = Params(nx=256, ny=256, Lx=28, Ly=28, p=3, fplot=false, dt=1, nstep=20, db=0.0005, bamp=0.05, wamp=0.05, freadinit=true, initfile="bwh.init.spot.dat", finalfile="clonal_exp.dat") # heavier setup more similar to the ones used for our simulations
b,w,h = main(P)
plotbwh(b, w, h, P, P.dt*100)
