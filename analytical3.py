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
    s = np.linspace(0,0.9999999999999,1000)
    y= s/(1-s)*np.exp(-s/(1-s)*z)/((s/(1-s))**2+1)/(1-s)**2
    a1= 42.242855
    a2=302.757865
    a3=352.018498
    a4= 21.821899
    b1= 48.196927
    b2= 482.485984
    b3= 1114.978885
    b4=449.690326
    if z<1:
      resultado = np.trapz(y,s)
    else:
      resultado = (1/z**2)*(z**8+a1*z**6+a2*z**4+a3*z**2+a4)/(z**8+b1*z**6+b2*z**4+b3*z**2+b4)

    return resultado
def f(z):
    s = np.linspace(0,0.99999999999,1000)
    y= np.exp(-s/(1-s)*z)/((s/(1-s))**2+1)/(1-s)**2
    a1,a2,a3,a4 = 38.027264, 265.187033, 335.677320, 38.102495
    b1, b2, b3, b4= 40.021433, 322.624911, 570.236280, 157.105423
    if z<1:
      resultado = np.trapz(y,s)
    else:
      resultado = (1/z)*(z**8+a1*z**6+a2*z**4+a3*z**2+a4)/(z**8+b1*z**6+b2*z**4+b3*z**2+b4)
    return resultado

def i(z):
    resultado = -1/z + f(z)
    return resultado


 
def normalize_angle(angle):
    return angle % (2 * np.pi)


def analytical3(P, tf, dt):

    r0,vphi0,vr0, alpha = P
    dr= vr0
    dphi= vphi0/r0
    phi0= 0
    h0 =r0**2*dphi
    p = h0**2/k
    e= np.sqrt(h0**4/k**2*((dr/(r0**2*dphi)-k/alpha**2*i(h0/alpha-phi0))**2+(1/r0-k/alpha**2*g(h0/alpha-phi0))**2))
    theta0 = phi0-np.arctan2(dr/(r0**2*dphi)  -  k/alpha**2*i(h0/alpha-phi0),(1/r0 - k/alpha**2*g(h0/alpha-phi0))  )


    solution1, time = drag_rk2_polar3(P, tf, dt)
    phi_an = solution1[:, 1]
    z_values = h0/alpha - phi_an

    g_vec = np.vectorize(g)
    g_results = g_vec(z_values)

    if alpha== 0.0:
        r_an = p/(1+ e*np.cos(phi_an))
    else:
        r_an = p / (e * np.cos(phi_an-theta0) + (h0 / alpha)**2* g_results)

    x_analytical = r_an * np.cos(phi_an)
    y_analytical = r_an * np.sin(phi_an)

    solution = np.column_stack((x_analytical, y_analytical))
    return solution, time
   