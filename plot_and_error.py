#the analytical solution is given only for a circular orbit, that is, v = np.sqrt(G*M/R) (a_centrifugal = a_gravitation) 
import numpy as np
import matplotlib.pyplot as plt

G = 1
M = 1

def plot_orbit_and_error(solution, time, x0, y0, vx0, vy0, tf, dt):
    x = solution[:, 0]
    y = solution[:, 1]
    vx = solution[:, 2]
    vy = solution[:, 3]

    r0 = np.sqrt(x0**2 + y0**2)
    v0 = np.sqrt(vx0**2 + vy0**2)
    phi0 = np.arctan2(y0, x0) 

    x_an = r0 * np.cos(v0/r0* time + phi0)
    y_an = r0 * np.sin(v0/r0 * time + phi0)

    plt.figure(figsize=(12, 6))
    ax = plt.gca() 
    ax.plot(x, y, label='RK2 Orbit', color= 'purple')
    ax.plot(x_an, y_an, label='Analitical Orbit', linestyle='--', color= 'red')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Orbit comparison')
    ax.grid(True)
    ax.legend()
    
    error_x = np.abs(x - x_an)
    error_y = np.abs(y - y_an)
    total_absolute_error = np.mean(np.sqrt(error_x**2 + error_y**2))

    ax.text(0.5, -0.1, 
            f'Total absolute error: {round(total_absolute_error * 100, 6)}%',
            transform=ax.transAxes,
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='center')
    
    plt.show()




