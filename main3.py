import argparse
from plot_orbit2 import plot_orbit2, polar_to_cartesian2
from euler2 import euler_polar2,  euler_cartesian2
from rk22 import rk2_cartesian2, rk2_polar2
from analytical3 import analytical3
from drag3 import drag_rk2_polar3, drag_rk2_cartesian3, drag_euler_cartesian3, drag_euler_polar3
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--e', type=float, required=True)
    parser.add_argument('--a', type=float, required=True) 
    parser.add_argument('--alpha', type=float, required=True) 

    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--method1', choices=['euler_cartesian2', 'euler_polar2', 'rk2_cartesian2', 'rk2_polar2', 'analytical3', 'drag_euler_cartesian3', 'drag_euler_polar3', 'drag_rk2_cartesian3', 'drag_rk2_polar3'])
    parser.add_argument('--method2', choices=['euler_cartesian2', 'euler_polar2', 'rk2_cartesian2', 'rk2_polar2', 'analytical3', 'drag_euler_cartesian3', 'drag_euler_polar3', 'drag_rk2_cartesian3', 'drag_rk2_polar3'])


    args = parser.parse_args()

    P = [args.e, args.a,args.alpha]

    if 'polar' in args.method1:
        solution1, time1 = globals()[args.method1](P, args.tf, args.dt)
    else:
        solution1, time1 = globals()[args.method1](P, args.tf, args.dt)

    solution2 = None
    if args.method2:
        if 'polar' in args.method2:
            solution2, time2 = globals()[args.method2](P, args.tf, args.dt)
        else:
            solution2, time2 = globals()[args.method2](P, args.tf, args.dt)

    plot_orbit2(time1, solution1, args.method1, solution2, args.method2)


if __name__ == "__main__":
    main()

#python3 main2.py --e 0 --a 1.0 --tf 10 --dt 0.1 --method1 euler_polar2 --method2 analytical2