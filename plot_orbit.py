import matplotlib.pyplot as plt
import numpy as np

def polar_to_cartesian(r, phi):
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return x, y

def plot_orbit(time1, solution1, method_name1, solution2=None, method_name2=None):
    plt.figure(figsize=(12, 6))

    if 'polar' in method_name1:
        r1 = solution1[:, 0]
        phi1 = solution1[:, 1]
        x1, y1 = polar_to_cartesian(r1, phi1)
    else:
        x1 = solution1[:, 0]
        y1 = solution1[:, 1]
    plt.plot(x1, y1, color='blue', label=f'Orbit with {method_name1}')

    if solution2 is not None:
        if 'polar' in method_name2:
            r2 = solution2[:, 0]
            phi2 = solution2[:, 1]
            x2, y2 = polar_to_cartesian(r2, phi2)
        else:
            x2 = solution2[:, 0]
            y2 = solution2[:, 1]
        plt.plot(x2, y2, color='red', linestyle='--',label=f'Orbit with {method_name2}')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Comparison of Orbits')
    plt.grid(True)
    plt.legend()
    plt.show()
