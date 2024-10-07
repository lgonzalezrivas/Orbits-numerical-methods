import numpy as np
import matplotlib.pyplot as plt
import argparse
from euler_polar import euler_polar
from euler_cartesian import euler_cartesian
from rk2_cartesian import rk2_cartesian
from rk2_polar import rk2_polar

m = 1 

def plot_angular_momentum(time1, solution1, method1, time2=None, solution2=None, method2=None):
    if 'polar' in method1:
        r1 = solution1[:, 0]
        vphi1 = solution1[:, 3]
        L1 = m * r1 * vphi1  
    else:
        x1 = solution1[:, 0]
        y1 = solution1[:, 1]
        vx1 = solution1[:, 2]
        vy1 = solution1[:, 3]
        L1 = m * (x1 * vy1 - y1 * vx1)  

    plt.figure(figsize=(10, 6))
    plt.plot(time1, L1, label=f'Angular momentum ({method1})', color='blue')

    if time2 is not None and solution2 is not None:
        if 'polar' in method2:
            r2 = solution2[:, 0]
            vphi2 = solution2[:, 3]
            L2 = m * r2 * vphi2  
        else:
            x2 = solution2[:, 0]
            y2 = solution2[:, 1]
            vx2 = solution2[:, 2]
            vy2 = solution2[:, 3]
            L2 = m * (x2 * vy2 - y2 * vx2) 

        plt.plot(time2, L2, label=f'Angular momentum ({method2})', color='red')

    plt.xlabel('Time')
    plt.ylabel('Angular Momentum')
    plt.title('Angular Momentum Evolution')
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=float, required=True)
    parser.add_argument('--c20', type=float, required=True)
    parser.add_argument('--vc10', type=float, required=True)
    parser.add_argument('--vc20', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--method1', type=str, choices=['euler_cartesian', 'rk2_cartesian', 'euler_polar', 'rk2_polar'], required=True)
    parser.add_argument('--method2', type=str, choices=['euler_cartesian', 'rk2_cartesian', 'euler_polar', 'rk2_polar'], required=True)
    args = parser.parse_args()

    IC = [args.c10, args.c20, args.vc10, args.vc20]
    tf = args.tf
    dt = args.dt  

    # Método 1
    if args.method1 == 'euler_cartesian':
        solution1, time1 = euler_cartesian(IC, tf, dt)
    elif args.method1 == 'rk2_cartesian':
        solution1, time1 = rk2_cartesian(IC, tf, dt)
    elif args.method1 == 'euler_polar':
        solution1, time1 = euler_polar(IC, tf, dt)
    elif args.method1 == 'rk2_polar':
        solution1, time1 = rk2_polar(IC, tf, dt)

    # Método 2
    if args.method2 == 'euler_cartesian':
        solution2, time2 = euler_cartesian(IC, tf, dt)
    elif args.method2 == 'rk2_cartesian':
        solution2, time2 = rk2_cartesian(IC, tf, dt)
    elif args.method2 == 'euler_polar':
        solution2, time2 = euler_polar(IC, tf, dt)
    elif args.method2 == 'rk2_polar':
        solution2, time2 = rk2_polar(IC, tf, dt)

    plot_angular_momentum(time1, solution1, args.method1, time2, solution2, args.method2)

if __name__ == "__main__":
    main()
#python3 momentum.py --c10 1.0 --c20 0.0 --vc10 0.0 --vc20 1.0 --tf 1000 --dt 0.1 --method1 rk2_cartesian --method2 rk2_polar