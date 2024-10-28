import argparse
from plot_orbit import plot_orbit, polar_to_cartesian
from euler import euler_polar, euler_cartesian
from rk2 import rk2_cartesian, rk2_polar
from analytical import analytical_polar, analytical_cartesian
from drag import drag_rk2_polar, drag_rk2_cartesian, drag_euler_cartesian, drag_euler_polar
import matplotlib.pyplot as plt
import numpy as np

m = 1 
def drag_analytical(IC,tf,dt):
    alpha= 0.002
    r0 = IC[0]
    vphi0 = IC[3]
    L0 =  m * r0*vphi0
    time = np.arange(0, tf, dt)
    L = L0*np.exp(-alpha*time)
    return L, time

def plot_angular_momentum(time1, solution1, method1, time2=None, solution2=None, method2=None):
    if method1 == 'drag_analytical':
        L1 = solution1

    elif 'polar' in method1 and 'analytical' not in method1:
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
        if method2 == 'drag_analytical':
            L2=solution2
        elif 'polar' in method2 and 'analytical' not in method2:
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
    parser.add_argument('--method1', choices=['euler_cartesian', 'euler_polar', 'rk2_cartesian', 'rk2_polar', 'drag_analytical', 'drag_euler_cartesian', 'drag_euler_polar', 'drag_rk2_cartesian', 'drag_rk2_polar'])
    parser.add_argument('--method2', choices=['euler_cartesian', 'euler_polar', 'rk2_cartesian', 'rk2_polar', 'drag_analytical', 'drag_euler_cartesian', 'drag_euler_polar', 'drag_rk2_cartesian', 'drag_rk2_polar'])
    args = parser.parse_args()

    IC = [args.c10, args.c20, args.vc10, args.vc20]
    tf = args.tf
    dt = args.dt  
    if args.method1 == 'drag_analytical':
        solution1, time1 = drag_analytical(IC, args.tf, args.dt)
    else:
        solution1, time1 = globals()[args.method1](IC, args.tf, args.dt)

    solution2, time2 = None, None
    if args.method2:
        if args.method2 == 'drag_analytical':
            solution2, time2 = drag_analytical(IC, args.tf, args.dt)
        else:
            solution2, time2 = globals()[args.method2](IC, args.tf, args.dt)

    plot_angular_momentum(time1, solution1, args.method1, time2, solution2, args.method2)

if __name__ == "__main__":
    main()
#python3 angular_momentum.py --c10 1 --c20 0 --vc10 0 --vc20 1.3 --tf 100 --dt 0.01 --method1 drag_euler_polar --method2 drag_analytical