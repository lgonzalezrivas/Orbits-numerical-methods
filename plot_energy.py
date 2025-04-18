import argparse
from plot_orbit import plot_orbit, polar_to_cartesian
from euler import euler_polar, euler_cartesian
from rk2 import rk2_cartesian, rk2_polar
from analytical import analytical_polar, analytical_cartesian
from drag import drag_rk2_polar, drag_rk2_cartesian, drag_euler_cartesian, drag_euler_polar
import matplotlib.pyplot as plt
import numpy as np

m = 1  

def initial_energy(IC, tf, dt):
    vc10 = IC[2]
    vc20 = IC[3]
    Ei = 0.5 * m * (vc10**2 + vc20**2)
    
    time = np.arange(0, tf, dt)
    energy_solution = np.full((len(time), 4), Ei)
    
    return energy_solution, time


def plot_energy(time1, solution1, method1, time2=None, solution2=None, method2=None):
    plt.figure(figsize=(10, 6))

    if method1 == 'initial_energy':
        Ek1 = solution1[:, 0]  
        plt.plot(time1, Ek1, label=f'Kinetic energy ({method1})', color='purple')

    else:
        vc1 = solution1[:, 2]
        vc2 = solution1[:, 3]
        Ek1 = 0.5 * m * (vc1**2 + vc2**2)
        plt.plot(time1, Ek1, label=f'Kinetic energy ({method1})', color='purple')

    if time2 is not None and solution2 is not None:
        if method2 == 'initial_energy':
            Ek2 = solution2[:, 0]  
            plt.plot(time2, Ek2, label=f'Kinetic energy ({method2})', color='red')
        else:
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
    parser.add_argument('--method1', type=str, choices=['euler_cartesian', 'rk2_cartesian', 'euler_polar', 'rk2_polar', 'initial_energy','drag_euler_cartesian', 'drag_euler_polar', 'drag_rk2_cartesian', 'drag_rk2_polar'])
    parser.add_argument('--method2', type=str, choices=['euler_cartesian', 'rk2_cartesian', 'euler_polar', 'rk2_polar', 'initial_energy','drag_euler_cartesian', 'drag_euler_polar', 'drag_rk2_cartesian', 'drag_rk2_polar'])
    args = parser.parse_args()

    IC = [args.c10, args.c20, args.vc10, args.vc20]
    tf = args.tf
    dt = args.dt  

    solution1, time1 = globals()[args.method1](IC, args.tf, args.dt)

    solution2, time2 = None, None
    if args.method2:
            solution2, time2 = globals()[args.method2](IC, args.tf, args.dt)
    plot_energy(time1, solution1, args.method1, time2, solution2, args.method2)

if __name__ == "__main__":
    main()
#python3 plot_energy.py --c10 1.000000 --c20 0 --vc10 0 --vc20 1.0 --tf 100 --dt 0.1 --method1 rk2_polar --method2 rk2_cartesian