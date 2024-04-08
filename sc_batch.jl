# -*- coding: utf-8 -*-
"""
Created on Monday 8 12:37:00 2024

@author: fiacac
"""
using DifferentialEquations
using Plots
using BenchmarkTools

# Constants
Yox_XG = 0.8
Yred_XG = 0.05
Yox_XE = 0.72
Y_OG = 1.067
Y_EG = 0.5
Y_OE = 1.5

q_g = 3.5
q_o = 0.37
q_e = 0.32
t_lag = 4.66

# Rates
Kg = 0.17
Ke = 0.56
Ko = 0.0001
Ki = 0.31

# Oxygen saturation
O_sat = 0.00755

# kla
kla = 1004
#k = 1
rho = ones(Float64,4)

# Stoichiometric matrix
# s = [ -1  0  -Y_OG  Yox_XG;
# -1  Y_EG  0  Yred_XG;
#  0  -1  -Y_OE  Yox_XE;
#  0   0   1   0 ]

# p = (Y_OG, Y_OE, q_g, q_o, q_e, t_lag,
#      Kg, Ke, Ko, Ki, O_sat, kla)#, s)

p = (Yox_XG, Yred_XG, Yox_XE, Y_OG, Y_EG, Y_OE, q_g, q_o, q_e, t_lag, Kg, Ke, Ko, Ki, O_sat, kla, rho)

# Reaction function
function rxn!(dC, C, p, t)
    # G0, E0, O0, X0, V0, T0 = C

    Yox_XG, Yred_XG, Yox_XE, Y_OG, Y_EG, Y_OE, q_g, q_o, q_e, t_lag, Kg, Ke, Ko, Ki, O_sat, kla, rho = p

    # Rate vector
    # rho = [ ((1 / Y_OG) * min(q_o * (C[3] / (C[3] + Ko)), Y_OG * (q_g * (C[1] / (C[1] + Kg)))) * C[4]);
    #         ((1 - exp(-t / t_lag)) * ((q_g * (C[1] / (C[1] + Kg))) - (1 / Y_OG) * min(q_o * (C[3] / (C[3] + Ko)), Y_OG * (q_g * (C[1] / (C[1] + Kg)))) * C[4]));
    #         ((1 / Y_OE) * min(q_o * (C[3] / (C[3] + Ko)) - (1 / Y_OG) * min(q_o * (C[3] / (C[3] + Ko)), Y_OE * (q_e * (C[2] / (C[2] + Ke)) * (Ki / (C[1] + Ki)))), Y_OE * (q_e * (C[2] / (C[2] + Ke)) * (Ki / (C[1] + Ki)))) * C[4]);
    #          kla * (O_sat - C[3])]# * k

    rho[1] = ((1 / Y_OG) * min(q_o * (C[3] / (C[3] + Ko)), Y_OG * (q_g * (C[1] / (C[1] + Kg)))) * C[4])
    rho[2] = ((1 - exp(-t / t_lag)) * ((q_g * (C[1] / (C[1] + Kg))) - (1 / Y_OG) * min(q_o * (C[3] / (C[3] + Ko)), Y_OG * (q_g * (C[1] / (C[1] + Kg)))) * C[4]))
    rho[3]= ((1 / Y_OE) * min(q_o * (C[3] / (C[3] + Ko)) - (1 / Y_OG) * min(q_o * (C[3] / (C[3] + Ko)), Y_OE * (q_e * (C[2] / (C[2] + Ke)) * (Ki / (C[1] + Ki)))), Y_OE * (q_e * (C[2] / (C[2] + Ke)) * (Ki / (C[1] + Ki)))) * C[4])
    rho[4] = kla * (O_sat - C[3])

    # Overall conversion vector
    # r = s * rho

    # Mass balances
    dC[1] = (-1*rho[1])+(-1*rho[2])+(0*rho[3])+(0*rho[4])
    dC[2] = (0*rho[1])+(Y_EG*rho[2])+(-1*rho[3])+(0*rho[4])
    dC[3] = (-Y_OG*rho[1])+(0*rho[2])+(-Y_OE*rho[3])+(1*rho[4])
    dC[4] = (Yox_XG*rho[1])+(Yred_XG*rho[2])+(Yox_XE*rho[3])+(0*rho[4])
    dC[5] = 0
    dC[6] = 0
    return nothing
end

# Solve function
tspan = (0.0, 30.0)
trange = LinRange(0.0, 30.0, 30)

C0 = [18, 0.0, 0.00755, 0.1, 2, 30]

prob = ODEProblem(rxn!, C0, tspan, p)
@btime sol = solve(prob, Tsit5(), saveat=trange, reltol=1e-7, abstol=1e-9)
plot(sol)


# @btime rxn!(zeros(Float64,6), C0, p, 0.0)