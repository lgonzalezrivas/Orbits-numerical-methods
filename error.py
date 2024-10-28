import argparse
from plot_orbit import plot_orbit, polar_to_cartesian
from euler import euler_polar, euler_cartesian
from rk2 import rk2_cartesian, rk2_polar
from analytical import analytical_polar, analytical_cartesian
from drag import drag_rk2_polar, drag_rk2_cartesian, drag_euler_cartesian, drag_euler_polar
import matplotlib.pyplot as plt
import numpy as np


G = 1
M = 1
m=1
k = G*M*m


def normalize_angle(angle):
    return angle % (2 * np.pi)

def totalerrors_polar(method, dt_values, IC, tf):
    errors = []

    for dt in dt_values:
        if method == 'rk2_polar':
            solution, time = rk2_polar(IC, tf, dt)
        elif method == 'drag_rk2_polar':
            solution, time = drag_rk2_polar(IC, tf, dt)
        elif method == 'drag_euler_polar':
            solution, time = drag_euler_polar(IC, tf, dt)
        elif method == 'euler_polar':
            solution, time = euler_polar(IC, tf, dt)

        r = solution[:, 0]
        phi = normalize_angle(solution[:, 1])  
        r0 = IC[0]
        phi0 =  normalize_angle(IC[1])
        vr0 = IC[2]
        vphi0= IC[3]
        L = m*r0 * vphi0 
        E = 0.5 * m*(vr0**2 + vphi0**2) - (G * M*m / r0)  
        error_L = None
        error_r = None
        error_phi = None
        a = -G * M*m / (2 * E)
        e = np.sqrt(1 + (2 * E * L**2) / (k**2))
        if method == 'drag_rk2_polar' or method == 'drag_euler_polar':
            r = solution[:,0]
            vphi= solution[:,3]
            L = m*r*vphi
            L0 = m*r0*vphi0
            alpha = 0.002
            L_an = L0*np.exp(-alpha*time)
            error_L = np.abs(L-L_an)
        elif e == 0.0 and 'drag' not in method:
            phi_an =  normalize_angle(np.sqrt(G*M/r0**3)*time +phi0)
            r_an = r0

            error_phi = np.abs(phi - phi_an)
            error_r = np.abs(r - r_an)

        elif 0.0 < e < 1.0 and 'drag' not in method:
            n = np.sqrt(G * (M) / a**3)  
            Mean = n * time  
            E = Mean  
            for _ in range(1000):  
                E_prev = E
                E = Mean + e * np.sin(E)
                if np.all(np.abs(E - E_prev) < 1e-10):  
                    break
            phi_an = normalize_angle(2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2)))
            r_an = a * (1 - e**2) / (1 + e * np.cos(phi_an))  
            error_phi = np.abs(phi-phi_an-phi0)
            error_r = np.abs(r - r_an)
        total_absolute_error = np.mean(np.sqrt(
            (error_phi if error_phi is not None else 0)**2 + 
            (error_r if error_r is not None else 0)**2 + 
            (error_L if error_L is not None else 0)**2
        ))
        errors.append(total_absolute_error)
    return errors



def totalerrors_cartesian(method, dt_values, IC, tf):
    errors = []

    for dt in dt_values:
        if method == 'rk2_cartesian':
            solution, time = rk2_cartesian(IC, tf, dt)
        if method == 'drag_rk2_cartesian':
            solution, time = drag_rk2_cartesian(IC, tf, dt)
        if method == 'drag_euler_cartesian':
            solution, time = drag_euler_cartesian(IC, tf, dt)
        elif method == 'euler_cartesian':
            solution, time = euler_cartesian(IC, tf, dt)

        x = solution[:, 0]
        y = solution[:, 1]  
        x0 = IC[0]
        y0 = IC[1]
        vx0 = IC[2]
        vy0= IC[3]
        error_L = None
        error_x = None
        error_y = None
        r0 = np.sqrt(x0**2 + y0**2)
        v0 = np.sqrt(vx0**2 + vy0**2)
        phi0 = np.arctan2(y0, x0)
        E = 0.5 * m*(v0)**2 - (G * M*m / r0)  
        L = m*(x0*vy0-y0*vx0)
        a = -G * M*m / (2 * E)
        e = np.sqrt(1 + (2 * E * L**2) / (m*k**2))
        if method == 'drag_rk2_cartesian' or method == 'drag_euler_cartesian':
            vx = solution[:,2]
            vy= solution[:,3]
            L = m*(x*vy-y*vx)
            L0 = m*(x0*vy0-y0*vx0)
            alpha = 0.002
            L_an = L0*np.exp(-alpha*time)
            error_L = np.abs(L-L_an)
        elif e == 0.0 and 'drag' not in method:
            x_an = r0 * np.cos(v0/r0 * time + phi0)
            y_an = r0 * np.sin(v0/r0 * time + phi0)
            error_x = np.abs(x - x_an)
            error_y = np.abs(y - y_an)
    
        elif 0.0 < e < 1.0 and 'drag' not in method:
            n = np.sqrt(G * (M) / a**3)  
            Mean = n * time  
            E = Mean  
            for _ in range(1000):  
                E_prev = E
                E = Mean + e * np.sin(E)
                if np.all(np.abs(E - E_prev) < 1e-10):  
                    break
            phi_an =  normalize_angle(2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))) 
            r_an = a * (1 - e**2) / (1 + e * np.cos(phi_an))
            x_an = r_an * np.cos(phi_an+ phi0)
            y_an = r_an * np.sin(phi_an+ phi0)  
            error_x = np.abs(x - x_an)
            error_y = np.abs(y - y_an)
        total_absolute_error = np.mean(np.sqrt(
            (error_x if error_x is not None else 0)**2 + 
            (error_y if error_y is not None else 0)**2 + 
            (error_L if error_L is not None else 0)**2
        ))
        errors.append(total_absolute_error)

    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=float, required=True)
    parser.add_argument('--c20', type=float, required=True)
    parser.add_argument('--vc10', type=float, required=True)
    parser.add_argument('--vc20', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--method1', choices=['euler_cartesian', 'euler_polar', 'rk2_cartesian', 'rk2_polar', 'drag_euler_cartesian', 'drag_euler_polar', 'drag_rk2_cartesian', 'drag_rk2_polar'])
    parser.add_argument('--method2', choices=['euler_cartesian', 'euler_polar', 'rk2_cartesian', 'rk2_polar',  'drag_euler_cartesian', 'drag_euler_polar', 'drag_rk2_cartesian', 'drag_rk2_polar'])
    args = parser.parse_args()

    dt_values = np.logspace(-4, -1, 10)
    IC = [args.c10, args.c20, args.vc10, args.vc20]

    plt.figure(figsize=(12, 6))
    if 'polar' in args.method1:
        errors1 = totalerrors_polar(args.method1, dt_values, IC, args.tf)
    else:
        errors1 = totalerrors_cartesian(args.method1,dt_values, IC, args.tf )
    plt.loglog(dt_values, errors1, marker='o', linestyle='--',color='blue', label=f'{args.method1}')
    errors2 = None
    if args.method2:
        if 'polar' in args.method2:
            errors2 = totalerrors_polar(args.method2, dt_values, IC, args.tf)
        else:
            errors2 = totalerrors_cartesian(args.method2,dt_values, IC, args.tf )    
        plt.loglog(dt_values, errors2, marker='o', linestyle='--', color='red',label=f'{args.method2}')



    plt.xlabel('dt')
    plt.ylabel('Absolute error')
    plt.title('Relation between error and dt')
    plt.grid(True)
    plt.legend()

    plt.show()

if __name__ == "__main__":
    main()

#python3 error_polar.py --c10 1.0 --c20 0.0 --vc10 0.0 --vc20 1.0 --tf 10.0 --method1 rk2_polar --method2 euler_polar