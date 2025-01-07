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

def g(z):
    s = np.linspace(0,0.999999,1000)
    y= s/(1-s)*np.exp(-s/(1-s)*z)/((s/(1-s))**2+1)/(1-s)**2
    a1= 42.242855
    a2=302.757865
    a3=352.018498
    a4= 21.821899
    b1= 48.196927
    b2= 482.485984
    b3= 1114.978885
    b4=449.690326
    if z<1:
      resultado = np.trapz(y,s)
    else:
      resultado = (1/z**2)*(z**8+a1*z**6+a2*z**4+a3*z**2+a4)/(z**8+b1*z**6+b2*z**4+b3*z**2+b4)

    return resultado
def f(z):
    s = np.linspace(0,0.999999,1000)
    y= np.exp(-s/(1-s)*z)/((s/(1-s))**2+1)/(1-s)**2
    a1,a2,a3,a4 = 38.027264, 265.187033, 335.677320, 38.102495
    b1, b2, b3, b4= 40.021433, 322.624911, 570.236280, 157.105423
    if z<1:
      resultado = np.trapz(y,s)
    else:
      resultado = (1/z)*(z**8+a1*z**6+a2*z**4+a3*z**2+a4)/(z**8+b1*z**6+b2*z**4+b3*z**2+b4)
    return resultado

def i(z):
    resultado = -1/z + f(z)
    return resultado



def normalize_angle(angle):
    return angle % (2 * np.pi)

def totalerrors_polar3(method, dt_values, P, tf):
    errors = []

    for dt in dt_values:
        if method == 'drag_rk2_polar3':
            solution, time = drag_rk2_polar3(P, tf, dt)
        elif method == 'drag_euler_polar3':
            solution, time = drag_euler_polar3(P, tf, dt)


        r = solution[:, 0]
        phi = normalize_angle(solution[:, 1])  
        r0,vphi0,vr0, alpha = P
        phi0=0
        dphi= vphi0/r0
        dr=vr0
        error_r = None
        if method == 'drag_rk2_polar3' or method == 'drag_euler_polar3':
            h0 =r0**2*dphi
            p = h0**2/k
            e= np.sqrt(h0**4/k**2*((dr/(r0**2*dphi)-k/alpha**2*i(h0/alpha-phi0))**2+(1/r0-k/alpha**2*g(h0/alpha-phi0))**2))
            theta0 = phi0-np.arctan2(dr/(r0**2*dphi)  -  k/alpha**2*i(h0/alpha-phi0),(1/r0 - k/alpha**2*g(h0/alpha-phi0))  )

            solution1, time = drag_rk2_polar3(P, tf, dt)
            phi_an = solution1[:, 1]
            z_values = h0/alpha - phi_an

            g_vec = np.vectorize(g)
            g_results = g_vec(z_values)

            if alpha== 0.0:
                r_an = p/(1+ e*np.cos(phi_an))
            else:
                r_an = p / (e * np.cos(phi_an-theta0) + (h0 / alpha)**2* g_results)


            error_r = np.abs(r-r_an)
        else: 
            aksdlja=3

        total_absolute_error = np.mean(np.sqrt(
            (error_r if error_r is not None else 0)**2 
        ))
        errors.append(total_absolute_error)
    return errors



def totalerrors_cartesian3(method, dt_values, P, tf):
    errors = []

    for dt in dt_values:
        if method == 'drag_rk2_cartesian3':
            solution, time = drag_rk2_cartesian3(P, tf, dt)
        elif method == 'drag_euler_cartesian3':
            solution, time = drag_euler_cartesian3(P, tf, dt)

        x = solution[:, 0]
        y = solution[:, 1]  
        r0,vphi0,vr0, alpha = P
        dphi= vphi0/r0
        dr=vr0
        h0=r0**2*dphi
        phi0=0
        theta0 = -np.arctan2(dr/(r0**2*dphi)  -  k/alpha**2*i(h0/alpha),(1/r0 - k/alpha**2*g(h0/alpha))  )
        vx=vr0
        vy=vphi0
        error_x = None
 

        if method == 'drag_rk2_cartesian3' or method == 'drag_euler_cartesian3':
            p = h0**2/k
            e= np.sqrt(h0**4/k**2*((dr/(r0**2*dphi)-k/alpha**2*i(h0/alpha-phi0))**2+(1/r0-k/alpha**2*g(h0/alpha-phi0))**2))


            solution1, time = drag_rk2_polar3(P, tf, dt)
            phi_an = solution1[:, 1]  
            z_values = h0/alpha - phi_an  

            g_vec = np.vectorize(g)
            g_results = g_vec(z_values) 
            if alpha== 0.0:
                r_an = p/(1+ e*np.cos(phi_an))
            else:
                r_an = p / (e * np.cos(phi_an-theta0) + (h0 / alpha)**2 * g_results)

            x_analytical = r_an * np.cos(phi_an)
            y_analytical = r_an * np.sin(phi_an)

            error_x = np.abs(x - x_analytical)
            error_y = np.abs(y - y_analytical)
        else:
            laskdjk= 2

        total_absolute_error = np.mean(np.sqrt(
            (error_x if error_x is not None else 0)**2 + 
            (error_y if error_y is not None else 0)**2 
        ))
        errors.append(total_absolute_error)

    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--r0', type=float, required=True)
    parser.add_argument('--vphi0', type=float, required=True)
    parser.add_argument('--vr0', type=float, required=True)
    parser.add_argument('--alpha', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=False)
    parser.add_argument('--method1', choices=['euler_cartesian2', 'euler_polar2', 'rk2_cartesian2', 'rk2_polar2', 'drag_euler_cartesian3', 'drag_euler_polar3', 'drag_rk2_cartesian3', 'drag_rk2_polar3'])
    parser.add_argument('--method2', choices=['euler_cartesian2', 'euler_polar2', 'rk2_cartesian2', 'rk2_polar2',  'drag_euler_cartesian3', 'drag_euler_polar3', 'drag_rk2_cartesian3', 'drag_rk2_polar3'])
    args = parser.parse_args()

    dt_values = np.logspace(-4, -1, 5)
    P = [args.r0, args.vphi0,args.vr0,args.alpha]

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

#python3 error3.py --e 0 --a 1 --alpha 0.02 --tf 10.0 --method1 drag_rk2_polar3 --method2 drag_rk2_cartesian3