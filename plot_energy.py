import numpy as np
import matplotlib.pyplot as plt
import argparse
from euler_polar import euler_polar
from euler_cartesian import euler_cartesian
from rk2_cartesian import rk2_cartesian
from rk2_polar import rk2_polar

m = 1  

def plot_energy(time, solution):
    vc1 = solution[:, 2] 
    vc2 = solution[:, 3]

    Ek = 0.5 * m * (vc1**2 + vc2**2)

    plt.figure(figsize=(10, 6))
    plt.plot(time, Ek, color='purple')
    plt.xlabel('Time')
    plt.ylabel('Kinetic energy')
    plt.title('Kinetic energy evolution')
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
    parser.add_argument('--method', type=str, choices=['euler_cartesian', 'rk2_cartesian', 'euler_polar', 'rk2_polar'], required=True)
    args = parser.parse_args()

    IC = [args.c10, args.c20, args.vc10, args.vc20]
    tf = args.tf
    dt = args.dt  

    if args.method == 'euler_cartesian':
        solution, time = euler_cartesian(IC, tf, dt)
    elif args.method == 'rk2_cartesian':
        solution, time = rk2_cartesian(IC, tf, dt)
    elif args.method == 'euler_polar':
        solution, time = euler_polar(IC, tf, dt)
    elif args.method == 'rk2_polar':
        solution, time = rk2_polar(IC, tf, dt)
    else:
        raise ValueError(f"Unknown method: {args.method}")

    plot_energy(time, solution)

if __name__ == "__main__":
    main()
#python3 plot_energy.py --c10 1.0 --c20 0.0 --vc10 0.0 --vc20 1.0 --tf 100.0 --dt 0.0001 --method rk2_polar