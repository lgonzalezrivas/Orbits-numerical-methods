import numpy as np
import matplotlib.pyplot as plt
import argparse
from euler_polar import euler_polar
from rk2_polar import rk2_polar

G, M, m = 1, 1, 1  

parser = argparse.ArgumentParser()
parser.add_argument('--c10', type=float, required=True)
parser.add_argument('--c20', type=float, required=True)
parser.add_argument('--vc10', type=float, required=True)
parser.add_argument('--vc20', type=float, required=True)
parser.add_argument('--tf', type=float, required=True)
parser.add_argument('--method1', type=str, choices=['rk2_polar', 'euler_polar'], required=True)
parser.add_argument('--method2', type=str, choices=['rk2_polar', 'euler_polar'], default=None)
args = parser.parse_args()


IC = [args.c10, args.c20, args.vc10, args.vc20]

dt_values = np.logspace(-3, 1, 20)

if args.method1 == 'euler_polar':
    run_simulation1 = euler_polar
elif args.method1 == 'rk2_polar':
    run_simulation1 = rk2_polar

results1 = {}
for dt in dt_values:
    solution, time = run_simulation1(IC, args.tf, dt)
    results1[dt] = (solution[:, 0], solution[:, 1]) 

if args.method2:
    if args.method2 == 'euler_polar':
        run_simulation2 = euler_polar
    elif args.method2 == 'rk2_polar':
        run_simulation2 = rk2_polar

    results2 = {}
    for dt in dt_values:
        solution, time = run_simulation2(IC, args.tf, dt)
        results2[dt] = (solution[:, 0], solution[:, 1])  

def calculate_error(results, dt_values, tf, c10):
    total_errors = []
    time_an = np.linspace(0, tf, 1000)
    phi_an = np.sqrt(G * M / c10**3) * time_an
    r_an = c10

    for dt, (r_values, phi_values) in results.items():
        phi_an_interp = np.interp(np.linspace(0, tf, len(phi_values)), time_an, phi_an)
        error_r = np.abs(r_values - r_an)
        error_phi = np.abs(phi_values - phi_an_interp)
        total_absolute_error = np.mean(error_r + error_phi)
        total_errors.append(total_absolute_error)

    return total_errors

errors_total1 = calculate_error(results1, dt_values, args.tf, args.c10)

if args.method2:
    errors_total2 = calculate_error(results2, dt_values, args.tf, args.c10)

plt.figure(figsize=(10, 6))
plt.loglog(dt_values, errors_total1, marker='o', linestyle='-', color='red', label=f'Method 1: {args.method1}')
if args.method2:
    plt.loglog(dt_values, errors_total2, marker='x', linestyle='--', color='blue', label=f'Method 2: {args.method2}')

plt.xlabel('dt (s)')
plt.ylabel('Total absolute error')
plt.title('Absolute error in polar coordinates')
plt.legend()
plt.grid(True)
plt.show()
#python3 error_polar.py --c10 1.0 --c20 0.0 --vc10 0.0 --vc20 1.0 --tf 100 --method1 rk2_polar --method2 euler_polar