import argparse
from runge_kutta import runge_kutta_2_cartesian
from plot_and_error import plot_orbit_and_error
from plot import plot


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--x0', type=float)
    parser.add_argument('--y0', type=float)
    parser.add_argument('--vx0', type=float)
    parser.add_argument('--vy0', type=float)
    parser.add_argument('--tf', type=float)
    parser.add_argument('--dt', type=float)
    
    args = parser.parse_args()

    IC = [args.x0, args.y0, args.vx0, args.vy0]
    solution, time = runge_kutta_2_cartesian(IC, args.tf, args.dt)

    plot(time, solution)

    #plot_orbit_and_error(solution, time, args.x0, args.y0, args.vx0, args.vy0, args.tf, args.dt) #circular  orbit


if __name__ == "__main__":
    main()

#python3 main.py --x0 1.0 --y0 0.0 --vx0 0.0 --vy0 1.0 --tf 3 --dt 0.01