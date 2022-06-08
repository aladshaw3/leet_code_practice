'''
    Use scipy.integrate.solve_ivp to integrate the function
'''
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import math

'''
def lotkavolterra(t, z, a, b, c, d):
    x, y = z
    return [a*x - b*x*y, -c*y + d*x*y]

sol = solve_ivp(lotkavolterra, [0, 15], [10, 5], args=(1.5, 1, 3, 1),
                dense_output=True)

t = np.linspace(0, 15, 300)
z = sol.sol(t)

plt.plot(t, z.T)
plt.xlabel('t')
plt.legend(['x', 'y'], shadow=True)
plt.title('Lotka-Volterra System')
plt.show()
'''

# Sample problem:
"""
    di/dt = sqrt(2)*V/L * sin(w*t+a) - R/L*i
"""
def current_in_circuit(t, i, V, L, w, a, R):
    return 2**0.5*V/L * np.sin(w*t+a) - (R/L)*i

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

dt = 0.5*min(max((1/w)*(math.pi/2 - a),(1/w)*(math.pi/2 + a)), (1/w)*(math.pi/2))

sol = solve_ivp(current_in_circuit, t_span=[t0, tf], y0=[io],
                args=(V, L, w, a, R), dense_output=True, max_step=dt,
                method='RK45', vectorized=True, first_step=dt)

print(sol)

# Plot just the calculated points
t = sol.t
z = sol.sol(t)
plt.plot(t, z.T, marker = 'o')

# Plot 1000s of points for a smooth representation of the exact soln
t = np.linspace(t0, tf, 1000)
z_exact = current_in_circuit_exact(t, V, L, w, a, R)
plt.plot(t, z_exact)


plt.xlabel('t')
plt.legend(['i_numerical', 'i_exact'], shadow=True)
plt.title('Current in Circuit')
plt.savefig('Result.png')
plt.show()
