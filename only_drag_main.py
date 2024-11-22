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



def polar_to_cartesian(r, phi):
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return x, y

def plot_orbit(time1, solution1, method_name1, solution2=None, method_name2=None):
    plt.figure(figsize=(12, 6))

    if 'polar' in method_name1 and 'analytical' not in method_name1:
        r1 = solution1[:, 0]
        phi1 = solution1[:, 1]
        x1, y1 = polar_to_cartesian(r1, phi1)
    else:
        x1 = solution1[:, 0]
        y1 = solution1[:, 1]
    plt.plot(x1, y1, color='purple', label=f'Orbit with {method_name1}')

    if solution2 is not None:
        if 'polar' in method_name2 and 'analytical' not in method_name2:
            r2 = solution2[:, 0]
            phi2 = solution2[:, 1]
            x2, y2 = polar_to_cartesian(r2, phi2)
        else:
            x2 = solution2[:, 0]
            y2 = solution2[:, 1]
        plt.plot(x2, y2, color='red', linestyle='--',label=f'Orbit with {method_name2}')
    plt.scatter(0, 0, color='red', s=150, label='Star')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Comparison of Orbits')
    plt.grid(True)
    plt.legend()
    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--c10', type=float, required=True)
    parser.add_argument('--c20', type=float, required=True) 
    parser.add_argument('--vc10', type=float, required=True)    
    parser.add_argument('--vc20', type=float, required=True)
    parser.add_argument('--tf', type=float, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--alpha', type=float, required=True)
    parser.add_argument('--method1', choices=['only_drag_analytical_polar', 'only_draganalytical_cartesian', 'only_drag_euler_cartesian', 'only_drag_euler_polar', 'only_drag_rk2_cartesian', 'only_drag_rk2_polar'])
    parser.add_argument('--method2', choices=['only_drag_analytical_polar', 'only_drag_analytical_cartesian', 'only_drag_euler_cartesian', 'only_drag_euler_polar', 'only_drag_rk2_cartesian', 'only_drag_rk2_polar'])


    args = parser.parse_args()

    IC = [args.c10, args.c20, args.vc10, args.vc20]

    if 'polar' in args.method1:
        solution1, time1 = globals()[args.method1](IC, args.tf, args.dt,args.alpha)
        method_type1 = 'polar'
    else:
        solution1, time1 = globals()[args.method1](IC, args.tf, args.dt,args.alpha)
        method_type1 = 'cartesian'

    solution2, method_type2 = None, None
    if args.method2:
        if 'polar' in args.method2:
            solution2, time2 = globals()[args.method2](IC, args.tf, args.dt,args.alpha)
            method_type2 = 'polar'
        else:
            solution2, time2 = globals()[args.method2](IC, args.tf, args.dt,args.alpha)
            method_type2 = 'cartesian'

    plot_orbit(time1, solution1, args.method1, solution2, args.method2)
# python3 only_drag_main.py --c10 10 --c20 0.291456794 --vc10 -10.0 --vc20 0.0 --tf 10 --alpha 0.1 --dt 0.1 --method1 only_drag_euler_polar --method2 only_drag_rk2_polar

if __name__ == "__main__":
    main()





