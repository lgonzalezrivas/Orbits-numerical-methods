import numpy as np
import time

G, M = 1, 1
def euler_cartesian2(P, tf, dt):
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

        if r<1e-1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with euler_cartesian.")
            return solution[:i], time[:i]
        else:
            ax = -G*M*x/r**3
            ay = -G*M*y/r**3
            
            x += dt*vx
            y += dt*vy

            vx += dt*ax
            vy += dt*ay


        solution[i] = [x,y,vx,vy]

    return solution, time



def euler_polar2(P, tf, dt):
    e,a, alpha = P
    phi=0
    r=a*(1-e)
    vr=0
    vphi= np.sqrt(G*M/a*(1+e)/(1-e))
    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr, vphi]

    for i in range(1, len(time)):
        if r<1e-1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with euler_polar.")
            return solution[:i], time[:i]
        else:
            dphi= vphi/r
            dr= vr

            d2r = r*dphi**2 - G*M/r**2
            d2phi= -2*dr * dphi/r


            r += dt * vr
            phi += dt * dphi
            vr +=  dt * d2r
            dphi+= dt* d2phi


        vphi=r*dphi
        solution[i] = [r, phi, vr, vphi]

    return solution, time




def main():
    P = [0, 1, 0.0011]
    tf = 30
    dt = 0.001

    start_time = time.time()  
    solution, time_vals = euler_cartesian2(P, tf, dt)
    end_time = time.time()  
    print(f"Tiempo de ejecución para euler_cartesian2: {end_time - start_time:.6f} segundos")

    start_time = time.time()  
    solution, time_vals = euler_polar2(P, tf, dt)
    end_time = time.time()  
    print(f"Tiempo de ejecución para euler_polar2: {end_time - start_time:.6f} segundos")

if __name__ == "__main__":
    main()
