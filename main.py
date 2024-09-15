import argparse
from plot_orbit import plot_orbit 
from plot_orbit import polar_to_cartesian
from euler_polar import euler_polar
from euler_cartesian import euler_cartesian
from rk2_cartesian import rk2_cartesian
from rk2_polar import rk2_polar
from analytical import analytical
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=float, required=True)
    parser.add_argument('--c20', type=float, required=True) 
    parser.add_argument('--vc10', type=float, required=True)    
    parser.add_argument('--vc20', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--method1', choices=['euler_cartesian', 'euler_polar', 'rk2_cartesian', 'rk2_polar', 'analytical'])
    parser.add_argument('--method2', choices=['euler_cartesian', 'euler_polar', 'rk2_cartesian', 'rk2_polar', 'analytical'])


    args = parser.parse_args()

    IC = [args.c10, args.c20, args.vc10, args.vc20]

    if 'polar' in args.method1:
        solution1, time1 = globals()[args.method1](IC, args.tf, args.dt)
        method_type1 = 'polar'
    else:
        solution1, time1 = globals()[args.method1](IC, args.tf, args.dt)
        method_type1 = 'cartesian'

    solution2, method_type2 = None, None
    if args.method2:
        if 'polar' in args.method2:
            solution2, time2 = globals()[args.method2](IC, args.tf, args.dt)
            method_type2 = 'polar'
        else:
            solution2, time2 = globals()[args.method2](IC, args.tf, args.dt)
            method_type2 = 'cartesian'

    plot_orbit(time1, solution1, args.method1, solution2, args.method2)



if __name__ == "__main__":
    main()

#python3 main.py --c10 1.0 --c20 0.0 --vc10 0.0 --vc20 1.0 --tf 1 --dt 0.1 --method1 euler_polar --method2 analytical