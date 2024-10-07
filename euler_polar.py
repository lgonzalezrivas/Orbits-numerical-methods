import numpy as np

G, M,m = 1, 1, 1  

def euler_polar(IC, tf, dt):
    r, phi, vr, vphi = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr, vphi]

    for i in range(1, len(time)):
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


