import numpy as np

G, M,m = 1, 1, 1  

def rk2_polar(IC, tf, dt):
    r, phi, vr, vphi = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr, vphi]
    l=m*r*vphi
    for i in range(1, len(time)):
        d2r = l**2/(m**2 *r**3) - G*M/r**2
        dphi= vphi/r
        k1 = np.array([vr, dphi, d2r])

        r_pred = r + dt * k1[0]
        phi_pred = phi + dt * k1[1]
        vr_pred = vr + dt * k1[2]
        d2r_pred = l**2/(m**2 *r_pred**3) - G*M/r_pred**2
        dphi_pred = l/(m*r_pred**2)
        k2 = np.array([vr_pred, dphi_pred, d2r_pred])

        r += 0.5 * dt * (k1[0] + k2[0])
        phi += 0.5 * dt * (k1[1] + k2[1])
        vr += 0.5 * dt * (k1[2] + k2[2])

        
        vphi=l/(m*r)
        solution[i] = [r, phi, vr, vphi]

    return solution, time


