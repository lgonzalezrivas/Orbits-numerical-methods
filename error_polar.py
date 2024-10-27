import numpy as np
import matplotlib.pyplot as plt
import argparse
from euler_polar import euler_polar
from rk2_polar import rk2_polar

G = 1
M = 1
m=1
k = G*M*m

def normalize_angle(angle):
    return angle % (2 * np.pi)

def totalerrors(method, dt_values, IC, tf):
    errors = []

    for dt in dt_values:
        if method == 'rk2_polar':
            solution, time = rk2_polar(IC, tf, dt)
        elif method == 'euler_polar':
            solution, time = euler_polar(IC, tf, dt)

        r = solution[:, 0]
        phi = normalize_angle(solution[:, 1])  
        r0 = IC[0]
        phi0 =  normalize_angle(IC[1])
        vr0 = IC[2]
        vphi0= IC[3]
        L = m*r0 * vphi0 
        E = 0.5 * m*(vr0**2 + vphi0**2) - (G * M*m / r0)  

        a = -G * M*m / (2 * E)
        e = np.sqrt(1 + (2 * E * L**2) / (k**2))
        if e == 0.0:
            phi_an =  normalize_angle(np.sqrt(G*M/r0**3)*time +phi0)
            r_an = r0

            error_phi = np.abs(phi - phi_an)
            error_r = np.abs(r - r_an)

        elif 0.0 < e < 1.0:
            n = np.sqrt(G * (M) / a**3)  
            Mean = n * time  
            E = Mean  
            for _ in range(1000):  
                E_prev = E
                E = Mean + e * np.sin(E)
                if np.all(np.abs(E - E_prev) < 1e-10):  
                    break
            phi_an = normalize_angle(2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2)))
            r_an = a * (1 - e**2) / (1 + e * np.cos(phi_an))  
            error_phi = np.abs(phi-phi_an-phi0)
            error_r = np.abs(r - r_an)
        
        total_absolute_error = np.mean(np.sqrt(error_phi**2+error_r**2))
        errors.append(total_absolute_error)
    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=float, required=True)
    parser.add_argument('--c20', type=float, required=True)
    parser.add_argument('--vc10', type=float, required=True)
    parser.add_argument('--vc20', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--method1', type=str, choices=['rk2_polar', 'euler_polar'], required=True)
    parser.add_argument('--method2', type=str, choices=['rk2_polar', 'euler_polar'], default=None)
    args = parser.parse_args()

    dt_values = np.logspace(-4, -1, 10)
    IC = [args.c10, args.c20, args.vc10, args.vc20]

    plt.figure(figsize=(12, 6))

    errors1 = totalerrors(args.method1, dt_values, IC, args.tf)
    plt.loglog(dt_values, errors1, marker='o', linestyle='--',color='blue', label=f'{args.method1}')

    if args.method2:
        errors2 = totalerrors(args.method2, dt_values, IC, args.tf)
        plt.loglog(dt_values, errors2, marker='o', linestyle='--', color='red',label=f'{args.method2}')



    plt.xlabel('dt')
    plt.ylabel('Absolute error')
    plt.title('Relation between error and dt')
    plt.grid(True)
    plt.legend()

    plt.show()

if __name__ == "__main__":
    main()

#python3 error_polar.py --c10 1.0 --c20 0.0 --vc10 0.0 --vc20 1.0 --tf 10.0 --method1 rk2_polar --method2 euler_polar