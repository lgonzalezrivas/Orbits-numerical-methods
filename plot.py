import numpy as np
import matplotlib.pyplot as plt

def plot(time, solution):
    x_pos = solution[:, 0]
    y_pos = solution[:, 1]
    plt.figure(figsize=(12, 6))
    plt.plot(x_pos, y_pos, color='purple')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Orbit with RK2')
    plt.grid(True)
    plt.show()
