import argparse
from plot_orbit2 import plot_orbit2, polar_to_cartesian2
from euler2 import euler_polar2,  euler_cartesian2
from rk22 import rk2_cartesian2, rk2_polar2
from analytical3 import analytical3
from drag3 import drag_rk2_polar3, drag_rk2_cartesian3, drag_euler_cartesian3, drag_euler_polar3
import matplotlib.pyplot as plt
import numpy as np

G = 1
M = 1
m=1
k = G*M*m
def regla_trapecio(f,a,b,n,z):
    h = (b - a)/n
    resultado = 1/2*(f(a,z) + f(b,z))
    for i in range(1, n):
        resultado += f(a + i * h, z)
    resultado *= h
    return resultado

def fun(s, z):
    return (s/ (s - 1)) *np.exp(-z * (s / (s - 1))) / ((s / (s - 1))**2 + 1)

def g(z, a, b, n):
    return regla_trapecio(fun, a, b, n, z)


def normalize_angle(angle):
    return angle % (2 * np.pi)

def totalerrors_polar3(method, dt_values, P, tf):
    errors = []

    for dt in dt_values:
        if method == 'rk2_polar2':
            solution, time = rk2_polar2(P, tf, dt)
        elif method == 'drag_rk2_polar3':
            solution, time = drag_rk2_polar3(P, tf, dt)
        elif method == 'drag_euler_polar3':
            solution, time = drag_euler_polar3(P, tf, dt)
        elif method == 'euler_polar2':
            solution, time = euler_polar2(P, tf, dt)

        r = solution[:, 0]
        phi = normalize_angle(solution[:, 1])  
        e,a, alpha = P
        phi0=0
        r0=a*(1-e)
        vr0=0
        vphi0= np.sqrt(G*M/a*(1+e)/(1-e))
        L = m*r0 * vphi0 
        E = 0.5 * m*(vr0**2 + vphi0**2) - (k / r0)  
        error_L = None
        error_r = None
        error_phi = None
        if method == 'drag_rk2_polar3' or method == 'drag_euler_polar3':
            e,a, alpha = P
            r = a*(1-e)
            phi = 0
            vphi = np.sqrt(G*M/a*(1-e)/(1+e))
            vr= 0
            num_t = int(tf / dt)
            time = np.linspace(0, tf, num_t)    
            r0 = a
            v0 = vphi
            phi0 = 0
            E = 0.5 * m*(v0)**2 - (G * M *m/ r0)  
            L = m*(r*vphi)
            h0=r0**2*v0/r0+alpha*phi0
            p = h0**2/(G*M*m)

            solution, time = rk2_polar2(P, tf, dt)
            r_an = p/(e*np.cos(phi-phi0)+(h0/alpha)**2*g(h0/alpha-phi,0,0.9,1000))
            error_r = np.abs(r-r_an)
        elif e == 0 and 'drag' not in method:
            phi_an =  normalize_angle(np.sqrt(G*M/r0**3)*time +phi0)
            r_an = r0

            error_phi = np.abs(phi - phi_an)
            error_r = np.abs(r - r_an)

        elif 0 < e < 1 and 'drag' not in method:
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
        total_absolute_error = np.mean(np.sqrt(
            (error_phi if error_phi is not None else 0)**2 + 
            (error_r if error_r is not None else 0)**2 + 
            (error_L if error_L is not None else 0)**2
        ))
        errors.append(total_absolute_error)
    return errors



def totalerrors_cartesian3(method, dt_values, P, tf):
    errors = []

    for dt in dt_values:
        if method == 'rk2_cartesian2':
            solution, time = rk2_cartesian2(P, tf, dt)
        if method == 'drag_rk2_cartesian3':
            solution, time = drag_rk2_cartesian3(P, tf, dt)
        if method == 'drag_euler_cartesian3':
            solution, time = drag_euler_cartesian3(P, tf, dt)
        elif method == 'euler_cartesian2':
            solution, time = euler_cartesian2(P, tf, dt)

        x = solution[:, 0]
        y = solution[:, 1]  
        phi_an= np.arctan2(y/x)
        e,a, alpha = P
        x0=a*(1-e)
        y0=0
        vx0=0
        vy0= np.sqrt(G*M/a*(1+e)/(1-e))
        error_L = None
        error_x = None
        error_y = None
        r0 = np.sqrt(x0**2 + y0**2)
        v0 = np.sqrt(vx0**2 + vy0**2)
        phi0 = np.arctan2(y0, x0)
        E = 0.5 * m*(v0)**2 - (k / r0)  
        L = m*(x0*vy0-y0*vx0)

        if method == 'drag_rk2_cartesian3' or method == 'drag_euler_cartesian3':
            e,a, alpha = P
            r = a*(1-e)
            phi = 0
            vphi = np.sqrt(G*M/a*(1-e)/(1+e))
            vr= 0
            num_t = int(tf / dt)
            time = np.linspace(0, tf, num_t)    
            r0 = a
            v0 = vphi
            phi0 = 0
            E = 0.5 * m*(v0)**2 - (G * M *m/ r0)  
            L = m*(r*vphi)
            h0=r0**2*v0/r0+alpha*phi0
            p = h0**2/(G*M*m)

            r_an = p/(e*np.cos(phi_an-phi0)+(h0/alpha)**2*g(h0/alpha-phi_an,0,0.9,1000))
            x_an = r_an*np.cos(phi_an+phi0)
            y_an=r_an*np.sin(phi_an+phi0)
            error_x = np.abs(x - x_an)
            error_y = np.abs(y - y_an)
    
        elif e == 0 and 'drag' not in method:
            x_an = r0 * np.cos(v0/r0 * time + phi0)
            y_an = r0 * np.sin(v0/r0 * time + phi0)
            error_x = np.abs(x - x_an)
            error_y = np.abs(y - y_an)
    
        elif 0 < e < 1 and 'drag' not in method:
            n = np.sqrt(G * (M) / a**3)  
            Mean = n * time  
            E = Mean  
            for _ in range(1000):  
                E_prev = E
                E = Mean + e * np.sin(E)
                if np.all(np.abs(E - E_prev) < 1e-10):  
                    break
            phi_an =  normalize_angle(2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))) 
            r_an = a * (1 - e**2) / (1 + e * np.cos(phi_an))
            x_an = r_an * np.cos(phi_an+ phi0)
            y_an = r_an * np.sin(phi_an+ phi0)  
            error_x = np.abs(x - x_an)
            error_y = np.abs(y - y_an)
        total_absolute_error = np.mean(np.sqrt(
            (error_x if error_x is not None else 0)**2 + 
            (error_y if error_y is not None else 0)**2 + 
            (error_L if error_L is not None else 0)**2
        ))
        errors.append(total_absolute_error)

    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--e', type=float, required=True)
    parser.add_argument('--a', type=float, required=True)
    parser.add_argument('--alpha', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=False)
    parser.add_argument('--method1', choices=['euler_cartesian2', 'euler_polar2', 'rk2_cartesian2', 'rk2_polar2', 'drag_euler_cartesian3', 'drag_euler_polar3', 'drag_rk2_cartesian3', 'drag_rk2_polar3'])
    parser.add_argument('--method2', choices=['euler_cartesian2', 'euler_polar2', 'rk2_cartesian2', 'rk2_polar2',  'drag_euler_cartesian3', 'drag_euler_polar3', 'drag_rk2_cartesian3', 'drag_rk2_polar3'])
    args = parser.parse_args()

    dt_values = np.logspace(-4, -1, 10)
    P = [args.e, args.a,args.alpha]

    plt.figure(figsize=(12, 6))
    if 'polar' in args.method1:
        errors1 = totalerrors_polar3(args.method1, dt_values, P, args.tf)
    else:
        errors1 = totalerrors_cartesian3(args.method1,dt_values, P, args.tf )
    plt.loglog(dt_values, errors1, marker='o', linestyle='--',color='blue', label=f'{args.method1}')
    errors2 = None
    if args.method2:
        if 'polar' in args.method2:
            errors2 = totalerrors_polar3(args.method2, dt_values, P, args.tf)
        else:
            errors2 = totalerrors_cartesian3(args.method2,dt_values, P, args.tf )    
        plt.loglog(dt_values, errors2, marker='o', linestyle='--', color='red',label=f'{args.method2}')
    


    if 'euler' in args.method1:
        min_error1 = min(errors1)
        slope1 = min_error1*(dt_values/dt_values[0])
        plt.loglog(dt_values, slope1, 'k--', label="Slope equal to 1")
    else: 
        min_error1 = min(errors1)
        slope1 = min_error1*(dt_values/dt_values[0])**2
        plt.loglog(dt_values, slope1, 'k--', label="Slope equal to 2")     

    if 'euler' in args.method2:
        min_error2 = min(errors2)
        slope2 = min_error2*(dt_values/dt_values[0])
        plt.loglog(dt_values, slope2, color='gray', linestyle='--', label="Slope equal to 1")
    else: 
        min_error2 = min(errors2)
        slope2 = min_error2*(dt_values/dt_values[0])**2
        plt.loglog(dt_values, slope2, color='gray', linestyle='--', label="Slope equal to 2")     

    plt.xlabel('dt')
    plt.ylabel('Absolute error')
    plt.title('Relation between error and dt')
    plt.grid(True)
    plt.legend()

    plt.show()

if __name__ == "__main__":
    main()

#python3 error2.py --e 0 --a 1 --alpha 0.02 --tf 10.0 --method1 rk2_polar2 --method2 euler_polar2--vc20 1.0