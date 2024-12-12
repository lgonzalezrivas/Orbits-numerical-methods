import numpy as np

G, M,m = 1, 1, 1  


def drag_rk2_polar3(P, tf, dt):
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
            d2r = r*dphi**2 - G*M/r**2 -alpha*dr/r**2
            d2phi= -2*dr * dphi/r -alpha*dphi/r**2
            k1 = np.array([dr, dphi, d2r, d2phi])

            r_pred = r + dt * k1[0]
            phi_pred = phi + dt * k1[1]
            dr_pred = dr + dt * k1[2]
            dphi_pred = dphi + dt * k1[3]
            d2r_pred = r_pred*dphi_pred**2 - G*M/r_pred**2-alpha*dr_pred/r_pred**2
            d2phi_pred = -2*dr_pred * dphi_pred/r_pred-alpha*dphi_pred/r_pred**2  
            k2 = np.array([dr_pred, dphi_pred, d2r_pred, d2phi_pred])

            r += 0.5 * dt *(k1[0] + k2[0])
            phi += 0.5 * dt * (k1[1] + k2[1])
            dr += 0.5 * dt * (k1[2] + k2[2])
            dphi+= 0.5 * dt* (k1[3] + k2[3]) 

        vr=dr
        vphi=r*dphi
        solution[i] = [r, phi, vr, vphi]

    return solution, time

def drag_rk2_cartesian3(P, tf, dt):
    e,a, alpha = P
    x=a*(1-e)
    y=0
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
            ax = -G*M*x/r**3 -alpha*vx/r**2
            ay = -G*M*y/r**3 -alpha*vy/r**2
            k1 = np.array([vx, vy, ax, ay])


            x_pred = x+dt*k1[0]
            y_pred = y+dt*k1[1]
            vx_pred = vx+dt*k1[2]
            vy_pred = vy+dt*k1[3]

            r_pred = np.sqrt(x_pred**2+y_pred**2)
            ax_pred= -G*M*x_pred/r_pred**3 -alpha*vx_pred/r_pred**2
            ay_pred = -G*M*y_pred/r_pred**3  -alpha*vy_pred/r_pred**2
            k2 = np.array([vx_pred,vy_pred, ax_pred, ay_pred])   

            x += 0.5 * dt * (k1[0] + k2[0])
            y += 0.5 * dt * (k1[1] + k2[1])
            vx += 0.5*dt*(k1[2] + k2[2])
            vy += 0.5*dt*(k1[3] + k2[3])

        solution[i] = [x,y,vx,vy]

    return solution, time

def drag_euler_cartesian3(P, tf, dt):
    e,a, alpha = P
    x=a*(1-e)
    y=0
    vy=np.sqrt(G*M/a*(1+e)/(1-e))
    vx=0
    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [x, y, vx, vy]
    
    for i in range(1, len(time)):
        r = np.sqrt(x**2 + y**2)
        if r<1e-1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with euler_cartesian.")

            return solution[:i], time[:i]
        else:
            ax = -G*M*x/r**3 -alpha*vx/r**2
            ay = -G*M*y/r**3 -alpha*vy/r**2
            
            x += dt*vx
            y += dt*vy

            vx += dt*ax
            vy += dt*ay


        solution[i] = [x,y,vx,vy]

    return solution, time



def drag_euler_polar3(P, tf, dt):
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

            d2r = r*dphi**2 - G*M/r**2 -alpha*dr/r**2
            d2phi= -2*dr * dphi/r -alpha*dphi/r**2


            r += dt * vr
            phi += dt * dphi
            vr +=  dt * d2r
            dphi+= dt* d2phi


        vphi=r*dphi
        solution[i] = [r, phi, vr, vphi]

    return solution, time


