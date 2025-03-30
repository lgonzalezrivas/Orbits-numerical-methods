import numpy as np
import time

G, M = 1, 1
def rk2_cartesian2(P, tf, dt):
    e,a, alpha = P
    y=0
    x=a*(1-e)
    vx=0
    vy= np.sqrt(G*M/a*(1+e)/(1-e))

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [x, y, vx, vy]
    
    for i in range(1, len(time)):
        r = np.sqrt(x**2 + y**2)
        if r<0.1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with rk2_cartesian.")
            return solution[:i], time[:i]
        else:
            ax = -G*M*x/r**3
            ay = -G*M*y/r**3
            k1 = np.array([vx, vy, ax, ay])


            x_pred = x+dt*k1[0]
            y_pred = y+dt*k1[1]
            vx_pred = vx+dt*k1[2]
            vy_pred = vy+dt*k1[3]

            r_pred = np.sqrt(x_pred**2+y_pred**2)
            ax_pred= -G*M*x_pred/r_pred**3
            ay_pred = -G*M*y_pred/r_pred**3  
            k2 = np.array([vx_pred,vy_pred, ax_pred, ay_pred])   

            x += 0.5 * dt * (k1[0] + k2[0])
            y += 0.5 * dt * (k1[1] + k2[1])
            vx += 0.5*dt*(k1[2] + k2[2])
            vy += 0.5*dt*(k1[3] + k2[3])

        solution[i] = [x,y,vx,vy]

    return solution, time

def rk2_polar2(P, tf, dt):
    e,a, alpha = P
    phi=0
    r=a*(1-e)
    vr=0
    vphi= np.sqrt(G*M/a*(1+e)/(1-e))

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr, vphi]

    for i in range(1, len(time)):
        if r<0.04:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with rk2_polar.")
            return solution[:i], time[:i]
        else:
            dphi= vphi/r
            dr= vr
            d2r = r*dphi**2 - G*M/r**2
            d2phi= -2*dr * dphi/r
            k1 = np.array([dr, dphi, d2r, d2phi])

            r_pred = r + dt * k1[0]
            phi_pred = phi + dt * k1[1]
            dr_pred = dr + dt * k1[2]
            dphi_pred = dphi + dt * k1[3]
            d2r_pred = r_pred*dphi_pred**2 - G*M/r_pred**2
            d2phi_pred = -2*dr_pred * dphi_pred/r_pred
            k2 = np.array([dr_pred, dphi_pred, d2r_pred, d2phi_pred])

            r += 0.5 * dt *(k1[0] + k2[0])
            phi += 0.5 * dt * (k1[1] + k2[1])
            dr += 0.5 * dt * (k1[2] + k2[2])
            dphi+= 0.5 * dt* (k1[3] + k2[3]) 

        vr=dr
        vphi=r*dphi
        solution[i] = [r, phi, vr, vphi]

    return solution, time



def main():
    P = [0, 1, 0.0011]
    tf = 30
    dt = 0.001

    start_time = time.time()  
    solution, time_vals = rk2_cartesian2(P, tf, dt)
    end_time = time.time()  
    print(f"Tiempo de ejecución para rk2_cartesian2: {end_time - start_time:.6f} segundos")

    start_time = time.time()  
    solution, time_vals = rk2_polar2(P, tf, dt)
    end_time = time.time()  
    print(f"Tiempo de ejecución para rk2_polar2: {end_time - start_time:.6f} segundos")

if __name__ == "__main__":
    main()