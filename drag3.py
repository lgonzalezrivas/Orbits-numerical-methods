import numpy as np

G, M,m = 1, 1, 1  
k=G*M*m

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






def drag_rk2_polar3(P, tf, dt):
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

    r,vphi,vr, alpha = P
    dphi= vphi/r
    phi=0
    dr=vr
    h0=r**2*dphi
    phi = -np.arctan2(dr/(r**2*dphi)  -  k/alpha**2*i(h0/alpha),(1/r - k/alpha**2*g(h0/alpha))  )
    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr,vphi]
    for i in range(1, len(time)):
        if r==-0.1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with rk2_polar.")
            return solution[:i], time[:i]
        else:
            dphi= vphi/r
            dr= vr
            d2r = r*dphi**2 - G*M/r**2 -alpha*dr/r**2
            d2phi= -2*dr * dphi/r -alpha*dphi/r**2
            k1 = np.array([dr, dphi, d2r, d2phi])

            r_pred = r + dt * k1[0]
            phi_pred = phi + dt * k1[1]
            dr_pred = dr + dt * k1[2]
            dphi_pred = dphi + dt * k1[3]
            d2r_pred = r_pred*dphi_pred**2 - G*M/r_pred**2-alpha*dr_pred/r_pred**2
            d2phi_pred = -2*dr_pred * dphi_pred/r_pred-alpha*dphi_pred/r_pred**2
            k2 = np.array([dr_pred, dphi_pred, d2r_pred, d2phi_pred])

            r += 0.5 * dt *(k1[0] + k2[0])
            phi += 0.5 * dt * (k1[1] + k2[1])
            dr += 0.5 * dt * (k1[2] + k2[2])
            dphi+= 0.5 * dt* (k1[3] + k2[3])

        vr=dr
        vphi=r*dphi
        solution[i] = [r, phi, vr, vphi]

    return solution, time


def drag_rk2_cartesian3(P, tf, dt):
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
    r,vphi,vr, alpha = P
    dphi= vphi/r
    phi=0
    dr=vr
    h0=r**2*dphi
    phi = -np.arctan2(dr/(r**2*dphi)  -  k/alpha**2*i(h0/alpha-phi),(1/r - k/alpha**2*g(h0/alpha-phi))  )
    vx=vr*np.cos(phi)-vphi*np.sin(phi)
    vy=vr*np.sin(phi)+vphi*np.cos(phi)
    print(vx    ,vy)
    x=r*np.cos(phi)
    y=r*np.sin(phi)
    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [x, y, vx, vy]
    
    for i in range(1, len(time)):
        if r==-0.1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with rk2_cartesian.")
            return solution[:i], time[:i]
        else:
            ax = -G*M*x/r**3 -alpha*vx/r**2
            ay = -G*M*y/r**3 -alpha*vy/r**2
            k1 = np.array([vx, vy, ax, ay])


            x_pred = x+dt*k1[0]
            y_pred = y+dt*k1[1]
            vx_pred = vx+dt*k1[2]
            vy_pred = vy+dt*k1[3]

            r_pred = np.sqrt(x_pred**2+y_pred**2)
            ax_pred= -G*M*x_pred/r_pred**3 -alpha*vx_pred/r_pred**2
            ay_pred = -G*M*y_pred/r_pred**3  -alpha*vy_pred/r_pred**2
            k2 = np.array([vx_pred,vy_pred, ax_pred, ay_pred])   

            x += 0.5 * dt * (k1[0] + k2[0])
            y += 0.5 * dt * (k1[1] + k2[1])
            vx += 0.5*dt*(k1[2] + k2[2])
            vy += 0.5*dt*(k1[3] + k2[3])

        solution[i] = [x,y,vx,vy]

    return solution, time

def drag_euler_cartesian3(P, tf, dt):
    r0,v0, alpha = P
    phi0= 0
    h0 =r0*v0+alpha*phi0
    p = h0**2/k

    x=r0
    y=0
    vx=0
    vy=v0
    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [x, y, vx, vy]
    
    for i in range(1, len(time)):
        r = np.sqrt(x**2 + y**2)
        if r<1e-1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with euler_cartesian.")

            return solution[:i], time[:i]
        else:
            ax = -G*M*x/r**3 -alpha*vx/r**2
            ay = -G*M*y/r**3 -alpha*vy/r**2
            
            x += dt*vx
            y += dt*vy

            vx += dt*ax
            vy += dt*ay


        solution[i] = [x,y,vx,vy]

    return solution, time



def drag_euler_polar3(P, tf, dt):
    r,v, alpha = P
    vr= 0
    vphi=v
    dphi= v/r
    phi=0 
    time = np.arange(0, tf, dt)
    solution = np.zeros((len(time), 4))
    solution[0] = [r, phi, vr,dphi]

    for i in range(1, len(time)):
        if r<1e-1:
            print(f"Simulation stopped at t={time[i-1]:.2f} because r has reached a value close to zero with euler_polar.")

            return solution[:i], time[:i]
        else:
            dphi= vphi/r
            dr= vr

            d2r = r*dphi**2 - G*M/r**2 -alpha*dr/r**2
            d2phi= -2*dr * dphi/r -alpha*dphi/r**2


            r += dt * vr
            phi += dt * dphi
            vr +=  dt * d2r
            dphi+= dt* d2phi


        vphi=r*dphi
        solution[i] = [r, phi, vr, vphi]

    return solution, time


