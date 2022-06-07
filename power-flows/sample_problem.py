"""
    Sample problem for solving a non-linear system that arrises
    in power-flow simulations (at steady-state). The purpose of
    this excerise is just to demonstrate knowledge of how to setup
    and solve a basic power-flow problem. This will not contain
    specific objects and logic with determining the Y-bus matrix.
"""
from scipy.optimize import least_squares
import numpy as np

'''
def fun_rosenbrock(x, mapping):
    X1 = x[1]
    X0 = x[0]
    return np.array([mapping["a"] * (X1 - X0**2), (1 - X0)])


mapping = [{"a": 10}]
x0_rosenbrock = np.array([2, 2])
res_1 = least_squares(fun_rosenbrock, x0_rosenbrock, args=mapping)
print(res_1)
'''

'''
    Power system of equations follows this logic:
        Pk = Vk * sum(over n, Vn * [Gkn*cos(dk-dn) + Bkn*sin(dk-dn)])
        Qk = Vk * sum(over n, Vn * [Gkn*sin(dk-dn) - Bkn*cos(dk-dn)])

    NOTE: This will be put in manually for testing, however, we can
        automate the formation of the residuals later. For now,
        the difficult part of automation is the fact that at some
        buses (k) different vars are known or uknown, so we would need
        some form of mapping to make it more clear.

    Unknown Var Map:
        x[0] = P0
        x[1] = Q0
        x[2] = V1
        x[3] = d1
        x[4] = Q2
        x[5] = d2
        x[6] = V3
        x[7] = d3
        x[8] = V4
        x[9] = d4
'''
def power_sys_fun(x, mapping):
    # Knowns (Given from circuit diagram)
    V0 = 1
    d0 = 0
    P1 = -8
    Q1 = -2.8
    V2 = 1.05
    P2 = 4.4
    P3 = 0
    Q3 = 0
    P4 = 0
    Q4 = 0

    # Unknowns (vars we are solving for, and how they manually map)
    P0 = x[0]
    Q0 = x[1]
    V1 = x[2]
    d1 = x[3]
    Q2 = x[4]
    d2 = x[5]
    V3 = x[6]
    d3 = x[7]
    V4 = x[8]
    d4 = x[9]

    # Vectors (These come from the Y-bus matrix)
    #       Can optimize my only storing half this info
    #       because Y-bus matrix is symmetric 
    G0 = np.array([3.73, 0, 0, 0, -3.73])
    B0 = np.array([-49.72, 0, 0, 0, 49.72])
    G1 = np.array([0, 2.68, 0, -0.89, -1.79])
    B1 = np.array([0, -28.46, 0, 9.92, 19.84])
    G2 = np.array([0, 0, 7.46, -7.46, 0])
    B2 = np.array([0, 0, -99.44, 99.44, 0])
    G3 = np.array([0, -0.89, -7.46, 11.92, -3.57])
    B3 = np.array([0, 9.92, 99.44, -147.96, 39.68])
    G4 = np.array([-3.73, -1.79, 0, -3.57, 9.09])
    B4 = np.array([49.72, 19.84, 0, 39.68, -108.58])

    d = np.array([d0, d1, d2, d3, d4])
    V = np.array([V0, V1, V2, V3, V4])
    P = np.array([P0, P1, P2, P3, P4])
    Q = np.array([Q0, Q1, Q2, Q3, Q4])

    G = np.array([G0, G1, G2, G3, G4])
    B = np.array([B0, B1, B2, B3, B4])

    # Vectorized forms of the residuals
    def _Pres(P, V, G, B, d, k):
        return P[k] - (V[k] * sum( V * (G[k]*np.cos(d[k]-d) + B[k]*np.sin(d[k]-d)) ))

    def _Qres(Q, V, G, B, d, k):
        return Q[k] - (V[k] * sum( V * (G[k]*np.sin(d[k]-d) - B[k]*np.cos(d[k]-d)) ))


    # Residuals
    res = np.array([0,0,0,0,0,0,0,0,0,0])

    res[0] = _Pres(P, V, G, B, d, 0)
    res[1] = _Qres(Q, V, G, B, d, 0)

    res[2] = _Pres(P, V, G, B, d, 1)
    res[3] = _Qres(Q, V, G, B, d, 1)

    res[4] = _Pres(P, V, G, B, d, 2)
    res[5] = _Qres(Q, V, G, B, d, 2)

    res[6] = _Pres(P, V, G, B, d, 3)
    res[7] = _Qres(Q, V, G, B, d, 3)

    res[8] = _Pres(P, V, G, B, d, 4)
    res[9] = _Qres(Q, V, G, B, d, 4)

    # NOTE: Did not allow me to return res for some reason
    return np.array([_Pres(P, V, G, B, d, 0),
                    _Qres(Q, V, G, B, d, 0),
                    _Pres(P, V, G, B, d, 1),
                    _Qres(Q, V, G, B, d, 1),
                    _Pres(P, V, G, B, d, 2),
                    _Qres(Q, V, G, B, d, 2),
                    _Pres(P, V, G, B, d, 3),
                    _Qres(Q, V, G, B, d, 3),
                    _Pres(P, V, G, B, d, 4),
                    _Qres(Q, V, G, B, d, 4)])

# Example of mapping (to use later)
mapping = [{"G": {(0,0): 3.73, (0,1): 0,} }]

# Solution (from literature):
#x0 = np.array([3.948, 1.144, 0.834, -0.39, 3.376-0.4, -0.0104, 1.019, -0.049462631, 0.974, -0.07941248])

# NOTE: Important to give 'good' initial guess.
#       A bad initial guess results in a different solution
#   This initial guess was same as given in literature and
#   we get to same solution as literature solution
x0 = np.array([0, 0, 1, 0, 0, 0, 1, 0, 1, 0])
'''
# Unknowns (map for reference)
P0 = x[0]
Q0 = x[1]
V1 = x[2]
d1 = x[3]
Q2 = x[4]
d2 = x[5]
V3 = x[6]
d3 = x[7]
V4 = x[8]
d4 = x[9]
'''
# Call solver
res_1 = least_squares(power_sys_fun, x0, args=mapping, ftol=1e-10, gtol=1e-10, xtol=1e-10, verbose=2)
# Print results
print()
print(res_1)
print()

# Check results by manually calling res func
x0 = res_1.x
print(power_sys_fun(x0, mapping))
