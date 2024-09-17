import numpy as np
import matplotlib.pyplot as plt
import argparse
from euler_polar import euler_polar
from euler_cartesian import euler_cartesian
from rk2_cartesian import rk2_cartesian
from rk2_polar import rk2_polar

m = 1  

def plot_energy(time1, solution1, method1, time2=None, solution2=None, method2=None):
    vc1 = solution1[:, 2] 
    vc2 = solution1[:, 3]

    Ek1 = 0.5 * m * (vc1**2 + vc2**2)
    
    plt.figure(figsize=(10, 6))
    plt.plot(time1, Ek1, label=f'Kinetic energy ({method1})', color='purple')

    if time2 is not None and solution2 is not None:
        vc3 = solution2[:, 2] 
        vc4 = solution2[:, 3]
        Ek2 = 0.5 * m * (vc3**2 + vc4**2)
        plt.plot(time2, Ek2, label=f'Kinetic energy ({method2})', color='red')

    plt.xlabel('Time')
    plt.ylabel('Kinetic energy')
    plt.title('Kinetic energy evolution')
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

    if args.method1 == 'euler_cartesian':
        solution1, time1 = euler_cartesian(IC, tf, dt)
    elif args.method1 == 'rk2_cartesian':
        solution1, time1 = rk2_cartesian(IC, tf, dt)
    elif args.method1 == 'euler_polar':
        solution1, time1 = euler_polar(IC, tf, dt)
    elif args.method1 == 'rk2_polar':
        solution1, time1 = rk2_polar(IC, tf, dt)

    if args.method2 == 'euler_cartesian':
        solution2, time2 = euler_cartesian(IC, tf, dt)
    elif args.method2 == 'rk2_cartesian':
        solution2, time2 = rk2_cartesian(IC, tf, dt)
    elif args.method2 == 'euler_polar':
        solution2, time2 = euler_polar(IC, tf, dt)
    elif args.method2 == 'rk2_polar':
        solution2, time2 = rk2_polar(IC, tf, dt)


    plot_energy(time1, solution1, args.method1, time2, solution2, args.method2)

if __name__ == "__main__":
    main()
#python3 plot_energy.py --c10 1.0 --c20 0.0 --vc10 0.0 --vc20 1.0 --tf 100.0 --dt 0.0001 --method1 euler_cartesian --method2 rk2_cartesian
