import numpy as np
import matplotlib.pyplot as plt
from rk2_cartesian import rk2_cartesian
from euler_cartesian import euler_cartesian
import argparse


G = 1
M = 1

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
        vy0 = IC[3]
        phi0 = np.arctan2(y0, x0) 

        r0 = np.sqrt(x0**2 + y0**2)
        v0 = np.sqrt(vx0**2 + vy0**2)
        x_an = r0 * np.cos(v0/r0 * time + phi0)
        y_an = r0 * np.sin(v0/r0 * time + phi0)

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

    dt_values = np.logspace(-3, 1, 20)
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
    plt.xlim(right=np.max(dt_values) * 0.6)

    plt.show()

if __name__ == "__main__":
    main()

#python3 error_cartesian.py --c10 1.0 --c20 0.0 --vc10 0.0 --vc20 1.0 --tf 10.0 --method1 rk2_cartesian --method2 euler_cartesian