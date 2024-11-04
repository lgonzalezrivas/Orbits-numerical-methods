import argparse
from plot_orbit2 import plot_orbit2, polar_to_cartesian2
from euler2 import euler_polar2,  euler_cartesian2
from rk22 import rk2_cartesian2, rk2_polar2
from analytical2 import analytical2
from drag2 import drag_rk2_polar2, drag_rk2_cartesian2, drag_euler_cartesian2, drag_euler_polar2
import matplotlib.pyplot as plt
import numpy as np

m = 1  
G, M=1,1

def initial_energy2(P, tf, dt):
    e,a, alpha = P
    vc10=0
    vc20= np.sqrt(G*M/a*(1+e)/(1-e))
    Ei = 0.5 * m * (vc10**2 + vc20**2)
    
    time = np.arange(0, tf, dt)
    energy_solution = np.full((len(time), 4), Ei)
    
    return energy_solution, time


def plot_energy2(time1, solution1, method1, time2=None, solution2=None, method2=None):
    plt.figure(figsize=(10, 6))

    if method1 == 'initial_energy2':
        Ek1 = solution1[:, 0]  
        plt.plot(time1, Ek1, label=f'Kinetic energy ({method1})', color='purple')

    else:
        vc1 = solution1[:, 2]
        vc2 = solution1[:, 3]
        Ek1 = 0.5 * m * (vc1**2 + vc2**2)
        plt.plot(time1, Ek1, label=f'Kinetic energy ({method1})', color='purple')

    if time2 is not None and solution2 is not None:
        if method2 == 'initial_energy2':
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
    parser.add_argument('--e', type=float, required=True)
    parser.add_argument('--a', type=float, required=True) 
    parser.add_argument('--alpha', type=float, required=True) 

    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--method1', type=str, choices=['euler_cartesian2', 'rk2_cartesian2', 'euler_polar2', 'rk2_polar2', 'initial_energy2','drag_euler_cartesian2', 'drag_euler_polar2', 'drag_rk2_cartesian2', 'drag_rk2_polar2'])
    parser.add_argument('--method2', type=str, choices=['euler_cartesian2', 'rk2_cartesian2', 'euler_polar2', 'rk2_polar2', 'initial_energy2','drag_euler_cartesian2', 'drag_euler_polar2', 'drag_rk2_cartesian2', 'drag_rk2_polar2'])
    args = parser.parse_args()

    P = [args.e, args.a,args.alpha]   
    tf = args.tf
    dt = args.dt  

    solution1, time1 = globals()[args.method1](P, args.tf, args.dt)

    solution2, time2 = None, None
    if args.method2:
            solution2, time2 = globals()[args.method2](P, args.tf, args.dt)
    plot_energy2(time1, solution1, args.method1, time2, solution2, args.method2)

if __name__ == "__main__":
    main()
#python3 plot_energ2y.py --e 0 --a 1.0 --tf 10 --dt 0.1--tf 100 --dt 0.1 --method1 rk2_polar2 --method2 rk2_cartesian2