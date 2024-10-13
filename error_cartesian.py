import numpy as np
import matplotlib.pyplot as plt
import argparse
from euler_cartesian import euler_cartesian
from rk2_cartesian import rk2_cartesian

G = 1
M = 1
m=1
k = G*M*m


def totalerrors(method, dt_values, IC, tf):
    errors = []

    for dt in dt_values:
        if method == 'rk2_cartesian':
            solution, time = rk2_cartesian(IC, tf, dt)
        elif method == 'euler_cartesian':
            solution, time = euler_cartesian(IC, tf, dt)

        x = solution[:, 0]
        y = solution[:, 1]  
        x0 = IC[0]
        y0 = IC[1]
        vx0 = IC[2]
        vy0= IC[3]
    
        r0 = np.sqrt(x0**2 + y0**2)
        v0 = np.sqrt(vx0**2 + vy0**2)
        phi0 = np.arctan2(y0, x0)
        E = 0.5 * m*(v0)**2 - (G * M*m / r0)  
        L = m*(x0*vy0-y0*vx0)
        a = -G * M*m / (2 * E)
        e = np.sqrt(1 + (2 * E * L**2) / (m*k**2))
        if e == 0.0:
            x_an = r0 * np.cos(v0/r0 * time + phi0)
            y_an = r0 * np.sin(v0/r0 * time + phi0)
            error_x = np.abs(x - x_an)
            error_y = np.abs(y - y_an)
    
        elif 0.0 < e < 1.0:
            n = np.sqrt(G * (M) / a**3)  
            Mean = n * time  
            E = Mean  
            for _ in range(1000):  
                E = Mean + e * np.sin(E)

            phi_an = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2)) 
            r_an = a * (1 - e**2) / (1 + e * np.cos(phi_an))
            x_an = r_an * np.cos(phi_an+ phi0)
            y_an = r_an * np.sin(phi_an+ phi0)  
            error_x = np.abs(x - x_an)
            error_y = np.abs(y - y_an)
    
        total_absolute_error = np.mean(np.sqrt(error_x**2 + error_y**2))
        errors.append(total_absolute_error)

    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=float, required=True)
    parser.add_argument('--c20', type=float, required=True)
    parser.add_argument('--vc10', type=float, required=True)
    parser.add_argument('--vc20', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--method1', type=str, choices=['rk2_cartesian', 'euler_cartesian'], required=True)
    parser.add_argument('--method2', type=str, choices=['rk2_cartesian', 'euler_cartesian'], default=None)
    args = parser.parse_args()

    dt_values = np.logspace(-4, 0, 20)
    IC = [args.c10, args.c20, args.vc10, args.vc20]

    plt.figure(figsize=(12, 6))

    errors1 = totalerrors(args.method1, dt_values, IC, args.tf)
    plt.loglog(dt_values, errors1, marker='o', linestyle='--',color='blue', label=f'{args.method1}')

    if args.method2:
        errors2 = totalerrors(args.method2, dt_values, IC, args.tf)
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