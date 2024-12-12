import argparse
from plot_orbit2 import plot_orbit2, polar_to_cartesian2
from euler2 import euler_polar2,  euler_cartesian2
from rk22 import rk2_cartesian2, rk2_polar2
from drag3 import drag_rk2_polar3, drag_rk2_cartesian3, drag_euler_cartesian3, drag_euler_polar3
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from numpy import inf, exp
G = 1
M = 1
m = 1
k = G * m * M



def g(z):
    s = np.linspace(0,0.999999,1000)
    y= s/(1-s)*np.exp(-s/(1-s)*z)/((s/(1-s))**2+1)/(1-s)**2
    resultado = np.trapz(y,s)
    return resultado
def normalize_angle(angle):
    return angle % (2 * np.pi)


def analytical3(P, tf, dt):

    e,a, alpha = P
    vphi = np.sqrt(G*M/a*(1-e)/(1+e))
    vr= 0
    r0 = a*(1-e)
    v0 = vphi
    phi0 = 0
    h0=np.sqrt(k*a*(1-e**2)-alpha*0.5*v0**(-2)*0)
    #h0 =r0*vphi+alpha*phi0
    p = h0**2/k

    solution1, time = drag_rk2_polar3(P, tf, dt)
    phi_an = normalize_angle(solution1[:, 1])  
    z_values = h0/alpha - phi_an  

    g_vec = np.vectorize(g)
    g_results = g_vec(z_values) 
    if alpha== 0.0:
        r_an = p/(1+ e*np.cos(phi_an))
    else:
        r_an = p / (e * np.cos(phi_an) + (h0 / alpha)**2 * g_results)

    x_analytical = r_an * np.cos(phi_an)
    y_analytical = r_an * np.sin(phi_an)

    solution = np.column_stack((x_analytical, y_analytical))
    return solution, time
   