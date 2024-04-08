"""
Created on Thu Sep 6 18:33:18 2018
Refactored on Monday 8 12:37:00 2024

@author: bjogut, simoca and fiacac
"""
from scipy.integrate import solve_ivp
import numpy as np
import time
import psutil
import matplotlib.pyplot as plt
import math

# Constants
Yox_XG = 0.8
Yred_XG = 0.05
Yox_XE = 0.72
Y_OG = 1.067
Y_EG = 0.5
Y_OE = 1.5

q_g = 4.68 #3.5
q_o = 0.37
q_e = 0.86 #0.32
t_lag = 4.66

# Rates
Kg = 0.17
Ke = 0.0001 #0.56
Ko = 0.56 #0.0001
Ki = 0.31

# Oxygen saturation
O_sat = 0.00755

# kla
kla = 1004

# Initial variables
G0 = 18
E0 = 0.34 #0.0
O0 = 0.00755
X0 = 0.1
V0 = 2
T0 = 30

# Time
t_end = 30
t_span = (0, t_end)

# fedbatch parameters
Cin = 100
F0 = 0.05 # flow can be modified
SFR = 0.15
t_expfb_start = 15
t_constfb_start = 18


# Reaction function
def rxn(t, C):
    G, E, O, X, V, T = C

    if (t >= t_expfb_start):
        if (t < t_constfb_start):
            Fin = F0 * math.exp(SFR * (t - t_expfb_start))
            Fout = 0
        else:
            Fin = F0 * math.exp(SFR * (t_constfb_start - t_expfb_start))
            Fout = 0
    else:
        Fin = 0
        Fout = 0

    F = Fin - Fout

    rho1 = ((1 / Y_OG) * min(q_o * (O / (O + Ko)), Y_OG * (q_g * (G / (G + Kg)))) * X)
    rho2 = ((1 - np.exp(-t / t_lag)) * ((q_g * (G / (G + Kg))) - (1 / Y_OG) * min(q_o * (O / (O + Ko)), Y_OG * (q_g * (G / (G + Kg)))) * X))
    rho3 = ((1 / Y_OE) * min(q_o * (O / (O + Ko)) - (1 / Y_OG) * min(q_o * (O / (O + Ko)), Y_OE * (q_e * (E / (E + Ke)) * (Ki / (G + Ki)))), Y_OE * (q_e * (E / (E + Ke)) * (Ki / (G + Ki)))) * X)
    rho4 = kla * (O_sat - O)

    r1 = (-1 * rho1) + (-1 * rho2) + (0 * rho3) + (0 * rho4)
    r2 = (0 * rho1) + (Y_EG * rho2) + (-1 * rho3) + (0 * rho4)
    r3 = (-Y_OG * rho1) + (0 * rho2) + (-Y_OE * rho3) + (1 * rho4)
    r4 = (Yox_XG * rho1) + (Yred_XG * rho2) + (Yox_XE * rho3) + (0 * rho4)


    dGdt = r1 - F / C[4] * C[0] + Fin / C[4] * Cin - Fout / C[4] * C[0]
    dEdt = r2 - F / C[4] * C[1] - Fout / C[4] * C[1]
    dOdt = r3
    dXdt = r4 - F / C[4] * C[3] - Fout / C[4] * C[3]
    dVdt = F
    dTdt = 0

    return [dGdt, dEdt, dOdt, dXdt, dVdt, dTdt]

# Solve function
def solve():
    sol = solve_ivp(rxn, t_span, [G0, E0, O0, X0, V0, T0], method='BDF', dense_output=True, max_step=1)
    return sol

# Main function
if __name__ == "__main__":
    process = psutil.Process()
    start_time = time.time()

    s = []
    for i in range(3):
        sol = solve()
        print("--- %s seconds ---" % (time.time() - start_time))
        s.append(time.time() - start_time)
        print(process.memory_info().rss * 0.001)

    print(np.mean(s))
    
    plt.plot(sol.t, sol.y[0, :], c='red', alpha=0.4)
    plt.plot(sol.t, sol.y[1, :], c='blue', alpha=0.4)
    plt.plot(sol.t, sol.y[2, :], c='green', alpha=0.4)
    plt.plot(sol.t, sol.y[3, :], c='orange', alpha=0.4)
    plt.legend()
    plt.grid()
    plt.show()
