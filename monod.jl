# -*- coding: utf-8 -*-
"""
Created on Monday 8 12:37:00 2024

@author: fiacac
"""
using DifferentialEquations
using BenchmarkTools
using Plots

# Define the parameters for the model
μ_max = 0.4     # Maximum growth rate (1/h)
Ks = 0.1        # Substrate half-saturation constant (g/L)
Yxs = 0.6       # Biomass yield coefficient (g/g)
P_max = 120.0   # Maximum product concentration (g/L)
μ = 1.0
r_s = 1.0
r_p = 1.0

p = (μ_max, Ks, Yxs, P_max, μ, r_s, r_p)
# Define the batch fermentation model
function fermentation_rate!(du, u, p, t)
    X, S, P = u
    μ_max, Ks, Yxs, P_max, μ, r_s, r_p = p

    μ = μ_max * S / (Ks + S)  # Monod equation for growth rate
    r_s = -μ / Yxs * X        # Substrate consumption rate
    r_p = μ * X               # Product formation rate

    # Rate of change
    du[1] = μ * X             # dX/dt
    du[2] = r_s               # dS/dt
    du[3] = r_p               # dP/dt
    return nothing 
end

# Initial conditions
X0 = 0.1        # Initial biomass concentration (g/L)
S0 = 20.0       # Initial substrate concentration (g/L)
P0 = 0.0        # Initial product concentration (g/L)
u0 = [X0, S0, P0]

# Time span for the simulation
tspan = (0.0, 50.0)
saveat = LinRange(0, 50, 101)
#print(tspan)

# Solve the differential equations
prob = ODEProblem(fermentation_rate!, u0, tspan, p)
#@btime sol = solve(prob, Tsit5(), saveat = saveat, reltol=1e-8, abstol=1e-8)
@btime sol = solve(prob, Tsit5(), saveat=saveat, reltol=1e-8, abstol=1e-8)
#runtime = @elapsed sol = solve(prob, Tsit5(), saveat=saveat, reltol=1e-8, abstol=1e-8)
# sol.u is the solution
# sol.t is the time 

# Output the solution
#println(sol)
plot(sol)