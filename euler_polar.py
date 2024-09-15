import numpy as np

G, M, m = 1, 1, 1  

def euler_polar(IC, tf, dt):
    r, phi, vr, vphi = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr, vphi]
    dphi = vphi / r
    l = m * r**2 * dphi

    for i in range(1, len(time)):
        dphi = l / (m * r**2)
        phi += dt * dphi
        vphi = dphi * r

        d2r = l**2 / (m**2 * r**3) - G * M / r**2
        vr += dt * d2r

        r += dt * vr

        solution[i] = [r, phi, vr, vphi]

    return solution, time


