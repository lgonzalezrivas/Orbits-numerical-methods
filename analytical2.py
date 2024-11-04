import numpy as np
import matplotlib.pyplot as plt

G = 1
M = 1
m = 1
k = G * m * M
def normalize_angle(angle):
    return angle % (2 * np.pi)
def analytical2(P, tf, dt):
    e,a, alpha = P
    r = a*(1-e)
    phi = 0
    vphi = np.sqrt(G*M/a*(1-e)/(1+e))
    vr= 0
    num_t = int(tf / dt)
    time = np.linspace(0, tf, num_t)    
    r0 = a
    v0 = vphi
    phi0 = 0
    E = 0.5 * m*(v0)**2 - (G * M *m/ r0)  
    L = m*(r*vphi)
    if e == 0:
        x_analytical = r0 * np.cos(v0/r0  * time + phi0)
        y_analytical = r0 * np.sin(v0/r0 * time + phi0)

    elif 0 < e < 1:
        n = np.sqrt(G * (M) / a**3)  
        Mean = n * time  
        
        E = Mean  

        for _ in range(1000):  
            E_prev = E
            E = Mean + e * np.sin(E)
            if np.all(np.abs(E - E_prev) < 1e-10):  
                break
        phi_an = normalize_angle(2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))) 
        r_an = a * (1 - e**2) / (1 + e * np.cos(phi_an))  
        
        x_analytical = r_an * np.cos(phi_an)
        y_analytical = r_an * np.sin(phi_an)

    solution = np.column_stack((x_analytical, y_analytical))
    return solution, time
   