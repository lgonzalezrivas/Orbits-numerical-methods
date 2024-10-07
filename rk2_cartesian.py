import numpy as np
G, M = 1, 1
def rk2_cartesian(IC, tf, dt):
    x, y , vx, vy = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [x, y, vx, vy]
    
    for i in range(1, len(time)):
        r = np.sqrt(x**2 + y**2)
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
