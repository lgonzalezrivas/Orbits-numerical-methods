import numpy as np

G, M,m = 1, 1, 1  
alpha = 0.002
def drag(IC, tf, dt):
    r, phi, vr, vphi = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr, vphi]

    for i in range(1, len(time)):
        dphi= vphi/r
        dr= vr
        d2r = r*dphi**2 - G*M/r**2 -alpha*dr
        d2phi= -2*dr * dphi/r -alpha*dphi
        k1 = np.array([dr, dphi, d2r, d2phi])

        r_pred = r + dt * k1[0]
        phi_pred = phi + dt * k1[1]
        dr_pred = dr + dt * k1[2]
        dphi_pred = dphi + dt * k1[3]
        d2r_pred = r_pred*dphi_pred**2 - G*M/r_pred**2-alpha*dr_pred
        d2phi_pred = -2*dr_pred * dphi_pred/r_pred-alpha*dphi_pred  
        k2 = np.array([dr_pred, dphi_pred, d2r_pred, d2phi_pred])

        r += 0.5 * dt *(k1[0] + k2[0])
        phi += 0.5 * dt * (k1[1] + k2[1])
        dr += 0.5 * dt * (k1[2] + k2[2])
        dphi+= 0.5 * dt* (k1[3] + k2[3]) 

        vr=dr
        vphi=r*dphi
        solution[i] = [r, phi, vr, vphi]

    return solution, time


