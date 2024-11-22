import argparse
import matplotlib.pyplot as plt
import numpy as np

G, M,m = 1, 1, 1  


def only_drag_rk2_polar(IC, tf, dt, alpha):
    r, phi, vr, vphi = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr, vphi]

    for i in range(1, len(time)):
        if r<0.04:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with rk2_polar.")
            return solution[:i], time[:i]
        else:
            dphi= vphi/r
            dr= vr
            d2r = r*dphi**2 -alpha*dr
            d2phi= -2*dr * dphi/r -alpha*dphi
            k1 = np.array([dr, dphi, d2r, d2phi])

            r_pred = r + dt * k1[0]
            phi_pred = phi + dt * k1[1]
            dr_pred = dr + dt * k1[2]
            dphi_pred = dphi + dt * k1[3]
            d2r_pred = r_pred*dphi_pred**2 -alpha*dr_pred
            d2phi_pred = -2*dr_pred * dphi_pred/r_pred-alpha*dphi_pred  
            k2 = np.array([dr_pred, dphi_pred, d2r_pred, d2phi_pred])

            r += 0.5 * dt *(k1[0] + k2[0])
            phi += 0.5 * dt * (k1[1] + k2[1])
            dr += 0.5 * dt * (k1[2] + k2[2])
            dphi+= 0.5 * dt* (k1[3] + k2[3]) 

        vr=dr
        vphi=r*dphi
        solution[i] = [r, phi, vr, vphi]

    return solution, time

def only_drag_rk2_cartesian(IC, tf, dt, alpha):
    x, y , vx, vy = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [x, y, vx, vy]
    
    for i in range(1, len(time)):
        r = np.sqrt(x**2 + y**2)
        if r<0.1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with rk2_cartesian.")
            return solution[:i], time[:i]
        else:
            ax =  -alpha*vx
            ay = -alpha*vy
            k1 = np.array([vx, vy, ax, ay])


            x_pred = x+dt*k1[0]
            y_pred = y+dt*k1[1]
            vx_pred = vx+dt*k1[2]
            vy_pred = vy+dt*k1[3]

            r_pred = np.sqrt(x_pred**2+y_pred**2)
            ax_pred=  -alpha*vx_pred
            ay_pred =  -alpha*vy_pred
            k2 = np.array([vx_pred,vy_pred, ax_pred, ay_pred])   

            x += 0.5 * dt * (k1[0] + k2[0])
            y += 0.5 * dt * (k1[1] + k2[1])
            vx += 0.5*dt*(k1[2] + k2[2])
            vy += 0.5*dt*(k1[3] + k2[3])

        solution[i] = [x,y,vx,vy]

    return solution, time

def only_drag_euler_cartesian(IC, tf, dt,alpha):
    x, y , vx, vy = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [x, y, vx, vy]
    
    for i in range(1, len(time)):
        r = np.sqrt(x**2 + y**2)
        if r<1e-1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with euler_cartesian.")

            return solution[:i], time[:i]
        else:
            ax =  -alpha*vx
            ay =  -alpha*vy
            
            x += dt*vx
            y += dt*vy

            vx += dt*ax
            vy += dt*ay


        solution[i] = [x,y,vx,vy]

    return solution, time



def only_drag_euler_polar(IC, tf, dt,alpha):
    r, phi, vr, vphi = IC

    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr, vphi]

    for i in range(1, len(time)):
        if r<1e-1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with euler_polar.")

            return solution[:i], time[:i]
        else:
            dphi= vphi/r
            dr= vr

            d2r = r*dphi**2 -alpha*dr
            d2phi= -2*dr * dphi/r -alpha*dphi


            r += dt * vr
            phi += dt * dphi
            vr +=  dt * d2r
            dphi+= dt* d2phi


        vphi=r*dphi
        solution[i] = [r, phi, vr, vphi]

    return solution, time

def totalerrors_polar(method, dt_values, IC, tf,alpha):
    errors = []
    for dt in dt_values:
        try:
            if method == 'only_drag_rk2_polar':
                solution, time = only_drag_rk2_polar(IC, tf, dt, alpha)
            elif method == 'only_drag_euler_polar':
                solution, time = only_drag_euler_polar(IC, tf, dt, alpha)

            r = solution[:, 0]
            phi = solution[:, 1]
            r0 = IC[0]
            phi0 = IC[1]
            vr0 = IC[2]
            vphi0 = IC[3]
            r_an = -vr0 * np.exp(-alpha * time) / alpha + r0 + vr0 / alpha
            #r_an = r0*np.exp(-alpha*time)/alpha
            phi_an = phi0
            error_phi = np.abs(phi - phi_an)
            error_r = np.abs(r - r_an)
            total_absolute_error = np.mean(np.sqrt(error_phi**2 + error_r**2))
            errors.append(total_absolute_error)
        except Exception as e:
            print(f"Error during simulation with dt={dt}: {e}")
            errors.append(np.nan)
    return errors






def totalerrors_cartesian(method, dt_values, IC, tf,alpha):
    errors = []

    for dt in dt_values:
        try:
            if method == 'only_drag_rk2_cartesian':
                solution, time = only_drag_rk2_cartesian(IC, tf, dt,alpha)
            if method == 'only_drag_euler_cartesian':
                solution, time = only_drag_euler_cartesian(IC, tf, dt,alpha)


            x = solution[:, 0]
            y = solution[:, 1]  
            x0 = IC[0]
            y0 = IC[1]
            vx0 = IC[2]
            vy0= IC[3]
            error_x = None
            error_y = None
            x_an = -vx0*np.exp(-alpha*time)/alpha + x0 + vx0/alpha
            y_an = -vy0*np.exp(-alpha*time)/alpha + y0 + vy0/alpha

            error_x = np.abs(x - x_an)
            error_y = np.abs(y - y_an)
            total_absolute_error = np.mean(np.sqrt(
                (error_x if error_x is not None else 0)**2 + 
                (error_y if error_y is not None else 0)**2
            ))
            errors.append(total_absolute_error)
        except Exception as e:
            print(f"Error during simulation with dt={dt}: {e}")
            errors.append(np.nan)
    return errors

def totalerrors_polar2(method, alpha_values, IC, tf,dt):
    errors = []
    for alpha in alpha_values:
        try:
            if method == 'only_drag_rk2_polar':
                solution, time = only_drag_rk2_polar(IC, tf, dt, alpha)
            elif method == 'only_drag_euler_polar':
                solution, time = only_drag_euler_polar(IC, tf, dt, alpha)

            r = solution[:, 0]
            phi = solution[:, 1]
            r0 = IC[0]
            phi0 = IC[1]
            vr0 = IC[2]
            vphi0 = IC[3]
            r_an = -vr0 * np.exp(-alpha * time) / alpha + r0 + vr0 / alpha
            #r_an = r0*np.exp(-alpha*time)/alpha
            phi_an = phi0
            error_phi = np.abs(phi - phi_an)
            error_r = np.abs(r - r_an)
            total_absolute_error = np.mean(np.sqrt(error_phi**2 + error_r**2))
            errors.append(total_absolute_error)
        except Exception as e:
            print(f"Error during simulation with dt={dt}: {e}")
            errors.append(np.nan)
    return errors

def totalerrors_cartesian2(method, alpha_values, IC, tf,dt):
    errors = []
    for alpha in alpha_values:
        try:
            if method == 'only_drag_rk2_cartesian':
                solution, time = only_drag_rk2_cartesian(IC, tf, dt,alpha)
            if method == 'only_drag_euler_cartesian':
                solution, time = only_drag_euler_cartesian(IC, tf, dt,alpha)


            x = solution[:, 0]
            y = solution[:, 1]  
            x0 = IC[0]
            y0 = IC[1]
            vx0 = IC[2]
            vy0= IC[3]
            error_x = None
            error_y = None
            x_an = -vx0*np.exp(-alpha*time)/alpha + x0 + vx0/alpha
            y_an = -vy0*np.exp(-alpha*time)/alpha + y0 + vy0/alpha

            error_x = np.abs(x - x_an)
            error_y = np.abs(y - y_an)
            total_absolute_error = np.mean(np.sqrt(
                (error_x if error_x is not None else 0)**2 + 
                (error_y if error_y is not None else 0)**2
            ))
            errors.append(total_absolute_error)
        except Exception as e:
            print(f"Error during simulation with dt={dt}: {e}")
            errors.append(np.nan)
    return errors






def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=float, required=True)
    parser.add_argument('--c20', type=float, required=True)
    parser.add_argument('--vc10', type=float, required=True)
    parser.add_argument('--vc20', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=True)

    parser.add_argument('--alpha', type=float, required=True)
    parser.add_argument('--method1', choices=['only_drag_euler_cartesian', 'only_drag_euler_polar', 'only_drag_rk2_cartesian', 'only_drag_rk2_polar'])
    parser.add_argument('--method2', choices=['only_drag_euler_cartesian', 'only_drag_euler_polar', 'only_drag_rk2_cartesian', 'only_drag_rk2_polar'])
    args = parser.parse_args()

    alpha_values = np.logspace(-5, 1, 7)
    dt_values = np.logspace(-4, -1, 10)
    IC = [args.c10, args.c20, args.vc10, args.vc20]
    alpha = args.alpha
    plt.figure(figsize=(12, 6))
    if 'polar' in args.method1:
        errors1 = totalerrors_polar(args.method1, dt_values, IC, args.tf,alpha)
    else:
        errors1 = totalerrors_cartesian(args.method1,dt_values, IC, args.tf,alpha )
    plt.loglog(dt_values, errors1, marker='o', linestyle='--',color='blue', label=f'{args.method1}')
    errors2 = None
    if args.method2:
        if 'polar' in args.method2:
            errors2 = totalerrors_polar(args.method2, dt_values, IC, args.tf,alpha)
        else:
            errors2 = totalerrors_cartesian(args.method2,dt_values, IC, args.tf,alpha )    
        plt.loglog(dt_values, errors2, marker='o', linestyle='--', color='red',label=f'{args.method2}')



    plt.xlabel('dt')
    plt.ylabel('Absolute error')
    plt.title('Relation between error and dt')
    plt.grid(True)
    plt.legend()

    plt.show()
    tf, dt = args.tf,args.dt
    plt.figure(figsize=(12, 6))
    if 'polar' in args.method1:
        errors1 = totalerrors_polar2(args.method1,alpha_values,  IC, tf,dt)
    else:
        errors1 = totalerrors_cartesian2(args.method1, alpha_values,  IC, tf,dt)
    plt.loglog(alpha_values, errors1, marker='o', linestyle='--', color='blue', label=f'{args.method1}')

    if args.method2:
        if 'polar' in args.method2:
            errors2 = totalerrors_polar2(args.method2, alpha_values,  IC, tf,dt)
        else:
            errors2 = totalerrors_cartesian2(args.method2, alpha_values,  IC, tf,dt)
        plt.loglog(alpha_values, errors2, marker='o', linestyle='--', color='red', label=f'{args.method2}')
    plt.xlabel('Alpha')
    plt.ylabel('Absolute error')
    plt.title('Relation between error and alpha')
    plt.grid(True)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()

#python3 only_drag_error.py --c10 30 --c20 0.0 --vc10 -1.0 --vc20 0.0 --tf 3 --alpha 0.1 --dt 0.1 --method1 only_drag_euler_polar --method2 only_drag_euler_cartesian