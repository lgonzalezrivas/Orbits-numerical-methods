import numpy as np
import matplotlib.pyplot as plt

G = 1
M = 1
m = 1
k = G * m * M
def normalize_angle(angle):
    return angle % (2 * np.pi)
def analytical_cartesian(IC, tf, dt):
    x0 = IC[0]
    y0 = IC[1]
    vx0 = IC[2]
    vy0= IC[3]
    num_t = int(tf / dt)
    time = np.linspace(0, tf, num_t)    
    r0 = np.sqrt(x0**2 + y0**2)
    v0 = np.sqrt(vx0**2 + vy0**2)
    phi0 = np.arctan2(y0, x0)
    E = 0.5 * m*(v0)**2 - (G * M *m/ r0)  
    L = m*(x0*vy0-y0*vx0)
    a = -G * M*m / (2 * E)
    e = np.sqrt(1 + (2 * E * L**2) / (k**2))

    if e == 0.0:
        x_analytical = r0 * np.cos(v0/r0  * time + phi0)
        y_analytical = r0 * np.sin(v0/r0 * time + phi0)

    elif 0.0 < e < 1.0:
        n = np.sqrt(G * (M) / a**3)  
        Mean = n * time  
        
        E = Mean  

        for _ in range(1000):  
            E = Mean + e * np.sin(E)
        phi_an = normalize_angle(2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))) 
        r_an = a * (1 - e**2) / (1 + e * np.cos(phi_an))  
        
        x_analytical = r_an * np.cos(phi_an+ phi0)
        y_analytical = r_an * np.sin(phi_an+ phi0)

    solution = np.column_stack((x_analytical, y_analytical))
    return solution, time
def analytical_polar(IC, tf, dt):
    r0 = IC[0]
    phi0 = IC[1]
    vr0 = IC[2]        
    vphi0 = IC[3]
    num_t = int(tf / dt)
    time = np.linspace(0, tf, num_t) 
    L = m*r0 * vphi0 
    E = 0.5 * m*(vr0**2 + vphi0**2) - (G * M*m / r0) 

    a = -G * M*m / (2 * E)
    e = np.sqrt(1 + (2 * E * L**2) / (k**2))  

    if e == 0.0:
        x_analytical = r0 * np.cos(vphi0 * time + phi0)
        y_analytical = r0 * np.sin(vphi0 * time + phi0)

    elif 0.0 < e < 1.0:
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
        
        x_analytical = r_an * np.cos(phi_an+ phi0)
        y_analytical = r_an * np.sin(phi_an+ phi0)

    solution = np.column_stack((x_analytical, y_analytical))
    return solution, time
