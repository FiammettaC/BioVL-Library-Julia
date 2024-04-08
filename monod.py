# -*- coding: utf-8 -*-
"""
Created on Monday 8 12:37:00 2024

@author: fiacac
"""
from scipy.integrate import odeint, solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import time
import psutil
import os
process = psutil.Process()
start_time = time.time()

# Parameters
μ_max = 0.4  # Maximum growth rate (1/h)
Ks = 0.1     # Substrate half-saturation constant (g/L)
Yxs = 0.6    # Biomass yield coefficient (g/g)
P_max = 120.0  # Maximum product concentration (g/L)

# Batch fermentation model
def fermentation_rate(t, u, μ_max, Yxs):
    X, S, P = u
    μ = μ_max * S / (Ks + S)  # Monod equation for growth rate
    r_s = -μ / Yxs * X        # Substrate consumption rate
    r_p = μ * X               # Product formation rate

    return [μ * X, r_s, r_p]

# Initial conditions
X0 = 0.1  # Initial biomass concentration (g/L)
S0 = 20.0  # Initial substrate concentration (g/L)
P0 = 0.0  # Initial product concentration (g/L)
u0 = np.array([X0, S0, P0])

# Time span for the simulation
t = np.linspace(0, 50, 100)  # 100 time points between 0 and 50 hours

# Solve the differential equations
#sol = odeint(fermentation_rate, u0, t)
sol = solve_ivp(fermentation_rate, y0=u0, t_span=(0, 50), args=(μ_max, Yxs), method="LSODA")
print("--- %s seconds ---" % (time.time() - start_time))
print(process.memory_info().rss / 1000)
# print(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)

# Output the solution
# for i in range(len(t)):
#     #print(f'Time: {t[i]:.1f} h, Biomass: {sol[i][0]:.3f} g/L, Substrate: {sol[i][1]:.3f} g/L, Product: {sol[i][2]:.3f} g/L')

plt.plot(sol.t, sol.y[0,:])
plt.plot(sol.t, sol.y[1,:])
plt.plot(sol.t, sol.y[2,:])
plt.show()
