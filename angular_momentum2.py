import argparse
from plot_orbit2 import plot_orbit2, polar_to_cartesian2
from euler2 import euler_polar2,  euler_cartesian2
from rk22 import rk2_cartesian2, rk2_polar2
from analytical2 import analytical2
from drag2 import drag_rk2_polar2, drag_rk2_cartesian2, drag_euler_cartesian2, drag_euler_polar2
import matplotlib.pyplot as plt
import numpy as np
G, M = 1,1
m = 1 
def drag_analytical2(P,tf,dt):
    e,a, alpha= P
    r0 = a*(1-e)
    vphi0 = np.sqrt(G*M/a*(1+e)/(1-e))
    L0 =  m * r0*vphi0
    time = np.arange(0, tf, dt)
    L = L0*np.exp(-alpha*time)
    return L, time

def plot_angular_momentum2(time1, solution1, method1, time2=None, solution2=None, method2=None):
    if method1 == 'drag_analytical2':
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
        if method2 == 'drag_analytical2':
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
    parser.add_argument('--e', type=float, required=True)
    parser.add_argument('--a', type=float, required=True)
    parser.add_argument('--alpha', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--method1', choices=['euler_cartesian2', 'euler_polar2', 'rk2_cartesian2', 'rk2_polar2', 'drag_analytical2', 'drag_euler_cartesian2', 'drag_euler_polar2', 'drag_rk2_cartesian2', 'drag_rk2_polar2'])
    parser.add_argument('--method2', choices=['euler_cartesian2', 'euler_polar2', 'rk2_cartesian2', 'rk2_polar2', 'drag_analytical2', 'drag_euler_cartesian2', 'drag_euler_polar2', 'drag_rk2_cartesian2', 'drag_rk2_polar2'])
    args = parser.parse_args()

    P = [args.e, args.a, args.alpha]
    tf = args.tf
    dt = args.dt  
    if args.method1 == 'drag_analytical2':
        solution1, time1 = drag_analytical2(P, args.tf, args.dt)
    else:
        solution1, time1 = globals()[args.method1](P, args.tf, args.dt)

    solution2, time2 = None, None
    if args.method2:
        if args.method2 == 'drag_analytical2':
            solution2, time2 = drag_analytical2(P, args.tf, args.dt)
        else:
            solution2, time2 = globals()[args.method2](P, args.tf, args.dt)

    plot_angular_momentum2(time1, solution1, args.method1, time2, solution2, args.method2)

if __name__ == "__main__":
    main()
#python3 angular_momentum2.py --e 0 --a 1 --alpha 1 --tf 100 --dt 0.01 --method1 drag_euler_polar2 --method2 drag_analytical2