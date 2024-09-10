import numpy as np
import matplotlib.pyplot as plt
from runge_kutta import runge_kutta_2_cartesian
import argparse


G = 1
M = 1

def totalerrors(dt_values, IC, tf):
    errors = []
    
    for dt in dt_values:
        solution, time = runge_kutta_2_cartesian(IC, tf, dt)
        x = solution[:, 0]
        y = solution[:, 1]  
        x0 = IC[0]
        y0 = IC[1]
        vx0= IC[2]
        vy0= IC[3]
        phi0 = np.arctan2(y0, x0) 

        r0 = np.sqrt(x0**2 + y0**2)
        v0 = np.sqrt(vx0**2 + vy0**2)
        x_an = r0 * np.cos(v0/r0* time + phi0)
        y_an = r0 * np.sin(v0/r0 * time + phi0)


        error_x = np.abs(x - x_an)
        error_y = np.abs(y - y_an)
        total_absolute_error = np.mean(np.sqrt(error_x**2 + error_y**2))
        errors.append(total_absolute_error)

    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--x0', type=float, required=True)
    parser.add_argument('--y0', type=float, required=True)
    parser.add_argument('--vx0', type=float, required=True)
    parser.add_argument('--vy0', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    args = parser.parse_args()

    dt_values= np.logspace(-3,1,20)
    IC = [args.x0, args.y0, args.vx0, args.vy0]
    errors = totalerrors(dt_values, IC, args.tf)

    plt.figure(figsize=(12, 6))
    plt.loglog(dt_values, errors, marker='o', linestyle='--')
    plt.xlabel('dt')
    plt.ylabel('Absolute error')
    plt.title('Relation between error and dt')
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()

#python3 evaluate_errors.py --x0 1.0 --y0 0.0 --vx0 0.0 --vy0 1.0 --tf 10.0
