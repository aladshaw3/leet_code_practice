from scipy.optimize import fsolve
import numpy as np
import math
import matplotlib.pyplot as plt

# Function to evaluate the next value
def evaluate_next(rate_fun,u_old,t,dt):
    u_new = u_old + rate_fun(t,u_old)*dt
    return u_new

def evaluate_next_alt(rate_fun,u_old,t,dt):

    # Define internal residual function
    def _res(u_new, rate_fun, u_old, t, dt):
        return (u_new-u_old)/dt - rate_fun(t, u_new)

    # Call a solver to find new based on residual analysis
    u_new = fsolve(_res, u_old, args=(rate_fun, u_old, t, dt))

    # Returns an array
    return u_new[0]

# Basic loop function for integration
def integration_loop(rate_fun,u0,t0,tf,dt):
    t = []
    u = []

    t.append(t0)
    u.append(u0)
    t_sim = t0

    i=0
    while t_sim < tf:
        t_sim = dt+t_sim
        i+=1

        t.append(t_sim)
        #u.append( evaluate_next(rate_fun,u[i-1],t_sim,dt) )
        u.append( evaluate_next_alt(rate_fun,u[i-1],t_sim,dt) )

    return t,u


def rate_fun_1(t, u):
    # Parameters
    w = 0.5           # Omega
    V = 20.           # Voltage
    a = -math.pi/4    # Alpha
    R = 0.8           # Resistance
    L = 8./w          # Impedence
    return math.sqrt(2)*(V/L)*math.sin(w*t+a) - (R/L)*u

if __name__=="__main__":

    t0= 0
    tf = 30
    io = 0.
    w = 0.5           # Omega
    a = -math.pi/4    # Alpha
    t_exp = []
    i_exp = []
    i_exp.append(io)
    t_exp.append(t0)
    t_sim = 0
    dt = 1

    dt = 0.5*min(max((1/w)*(math.pi/2 - a),(1/w)*(math.pi/2 + a)), (1/w)*(math.pi/2))

    t_exp, i_exp = integration_loop(rate_fun_1,io,t0,tf,dt)

    # Plot just the calculated points
    plt.plot(t_exp, i_exp, marker = '.')

    plt.xlabel('t')
    plt.legend(['i_numerical'], shadow=True)
    #plt.legend(['i_exact', 'i_stable_explicit'], shadow=True)
    plt.title('Current in Circuit')
    plt.savefig('Result.png')
    plt.show()
