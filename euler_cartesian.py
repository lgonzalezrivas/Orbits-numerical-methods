import numpy as np
G, M = 1, 1
def euler_cartesian(IC, tf, dt):
    x, y , vx, vy = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [x, y, vx, vy]
    
    for i in range(1, len(time)):
        r = np.sqrt(x**2 + y**2)
        ax = -G*M*x/r**3
        ay = -G*M*y/r**3


        vx += dt*ax
        vy += dt*ay
        x += dt*vx
        y += dt*vy

        solution[i] = [x,y,vx,vy]

    return solution, time
