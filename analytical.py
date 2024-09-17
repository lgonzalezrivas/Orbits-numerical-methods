#the analytical solution is given only for a circular orbit, that is, v = np.sqrt(G*M/R) (a_centrifugal = a_gravitation) 
import numpy as np
import matplotlib.pyplot as plt

G = 1
M = 1

 
def analytical(IC, tf, dt):
    x0= IC[0]
    y0= IC[1]
    vx0= IC[2]        
    vy0= IC[3]
    r0 = np.sqrt(x0**2 + y0**2)  
    v0 = np.sqrt(vx0**2 + vy0**2)
    phi0 = np.arctan2(y0, x0) 
    time = np.linspace(0, tf, 1000)
    
    x_analytical = r0 * np.cos(v0/r0 * time+phi0)
    y_analytical = r0 * np.sin(v0/r0 * time+phi0)
    
    solution = np.column_stack((x_analytical, y_analytical))
    return solution, time



