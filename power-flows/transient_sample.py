'''
    Use scipy.integrate.solve_ivp to integrate the function
'''
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import math


# Sample problem:
"""
    di/dt = sqrt(2)*V/L * sin(w*t+a) - R/L*i
"""
def current_in_circuit(t, i, V, L, w, a, R):
    return 2**0.5*V/L * np.sin(w*t+a) - (R/L)*i

def explicit_step(i_old, dt, t, V, L, w, a, R):
    return i_old + dt*current_in_circuit(t, i_old, V, L, w, a, R)

def current_in_circuit_exact(t, V, L, w, a, R):
    Z = (R**2 + (w*L)**2)**0.5
    theta = np.arctan(w*L/R)
    T = L/R
    return 2**0.5*V/Z*(np.sin(w*t+a-theta) - np.sin(a-theta)*np.exp(-t/T))

# Parameters
w = 0.5
V = 20.
a = -math.pi/4
R = 0.8
io = 0.
L = 8./w
t0= 0
tf = 30



t_exp = []
i_exp = []
i_exp.append(io)
t_exp.append(t0)
t_sim = 0
dt = 1

dt = 0.5*min(max((1/w)*(math.pi/2 - a),(1/w)*(math.pi/2 + a)), (1/w)*(math.pi/2))

i=0
while t_sim < tf:
    t_sim = dt+t_sim
    i+=1

    t_exp.append(t_sim)
    i_exp.append(explicit_step(i_exp[i-1], dt, t_exp[i], V, L, w, a, R))


sol = solve_ivp(current_in_circuit, t_span=[t0, tf], y0=[io],
                args=(V, L, w, a, R), dense_output=True, max_step=dt,
                method='RK45', vectorized=True, first_step=dt)


# Plot just the calculated points
t = sol.t
z = sol.sol(t)
plt.plot(t, z.T, marker = 'o')

# Plot 1000s of points for a smooth representation of the exact soln
t = np.linspace(t0, tf, 1000)
z_exact = current_in_circuit_exact(t, V, L, w, a, R)
plt.plot(t, z_exact)

#plt.plot(t_exp, i_exp, marker = '.')

plt.xlabel('t')
plt.legend(['i_numerical', 'i_exact', 'i_bad'], shadow=True)
#plt.legend(['i_exact', 'i_stable_explicit'], shadow=True)
plt.title('Current in Circuit')
plt.savefig('Result.png')
plt.show()
