import argparse
from plot_orbit2 import plot_orbit2, polar_to_cartesian2
from euler2 import euler_polar2,  euler_cartesian2
from rk22 import rk2_cartesian2, rk2_polar2
from analytical2 import analytical2
from drag2 import drag_rk2_polar2, drag_rk2_cartesian2, drag_euler_cartesian2, drag_euler_polar2
import matplotlib.pyplot as plt
import numpy as np

G = 1
M = 1
m=1
k = G*M*m


def normalize_angle(angle):
    return angle % (2 * np.pi)

def totalerrors_polar2(method, alpha_values, a, e, tf, dt):
    errors = []
    for alpha in alpha_values:
        P = [e, a, alpha]
        if method == 'drag_rk2_polar2':
            solution, time = drag_rk2_polar2(P, tf, dt)
        elif method == 'drag_euler_polar2':
            solution, time = drag_euler_polar2(P, tf, dt)

        r = solution[:, 0]
        phi = normalize_angle(solution[:, 1])  
        e,a, alpha = P
        phi0=0
        r0=a*(1-e)
        vr0=0
        vphi0= np.sqrt(G*M/a*(1+e)/(1-e))
        L = m*r0 * vphi0 
        E = 0.5 * m*(vr0**2 + vphi0**2) - (k / r0)  
        error_L = None
        r = solution[:,0]
        vphi= solution[:,3]
        L = m*r*vphi
        L0 = m*r0*vphi0
        L_an = L0*np.exp(-alpha*time)
        error_L = np.abs(L-L_an)
        total_absolute_error = np.mean(error_L)
        errors.append(total_absolute_error)
    return errors



def totalerrors_cartesian2(method, alpha_values, a, e, tf, dt):
    errors = []
    for alpha in alpha_values:
        P = [e, a, alpha]
        if method == 'drag_rk2_cartesian2':
            solution, time = drag_rk2_cartesian2(P, tf, dt)
        elif method == 'drag_euler_cartesian2':
            solution, time = drag_euler_cartesian2(P, tf, dt)
        x = solution[:, 0]
        y = solution[:, 1]  
        e,a, alpha = P
        x0=a*(1-e)
        y0=0
        vx0=0
        vy0= np.sqrt(G*M/a*(1+e)/(1-e))
        error_L = None
        r0 = np.sqrt(x0**2 + y0**2)
        v0 = np.sqrt(vx0**2 + vy0**2)
        phi0 = np.arctan2(y0, x0)
        E = 0.5 * m*(v0)**2 - (k / r0)  
        L = m*(x0*vy0-y0*vx0)

        vx = solution[:,2]
        vy= solution[:,3]
        L = m*(x*vy-y*vx)
        L0 = m*(x0*vy0-y0*vx0)
        L_an = L0*np.exp(-alpha*time)
        error_L = np.abs(L-L_an)
        total_absolute_error = np.mean(error_L)
        errors.append(total_absolute_error)
    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--e', type=float, required=True)
    parser.add_argument('--a', type=float, required=True)

    parser.add_argument('--alpha', type=float, required=False)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--method1', choices=['drag_euler_cartesian2', 'drag_euler_polar2', 'drag_rk2_cartesian2', 'drag_rk2_polar2'])
    parser.add_argument('--method2', choices=['drag_euler_cartesian2', 'drag_euler_polar2', 'drag_rk2_cartesian2', 'drag_rk2_polar2'])
    args = parser.parse_args()

    alpha_values = np.logspace(-5, 1, 7)
    e= args.e
    a = args.a
    tf = args.tf
    dt = args.dt

    plt.figure(figsize=(12, 6))
    if 'polar' in args.method1:
        errors1 = totalerrors_polar2(args.method1,alpha_values, a, e, tf, dt)
    else:
        errors1 = totalerrors_cartesian2(args.method1, alpha_values, a, e, tf, dt)
    plt.loglog(alpha_values, errors1, marker='o', linestyle='--', color='blue', label=f'{args.method1}')

    if args.method2:
        if 'polar' in args.method2:
            errors2 = totalerrors_polar2(args.method2, alpha_values, a, e, tf, dt)
        else:
            errors2 = totalerrors_cartesian2(args.method2, alpha_values, a, e, tf, dt)
        plt.loglog(alpha_values, errors2, marker='o', linestyle='--', color='red', label=f'{args.method2}')

    plt.xlabel('Alpha')
    plt.ylabel('Absolute error')
    plt.title('Relation between error and alpha')
    plt.grid(True)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()

#python3 error2.py --e 0 --a 1 --alpha 0.02 --dt 0.1 --tf 10.0 --method1 rk2_polar2 --method2 euler_polar2