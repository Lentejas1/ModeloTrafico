import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

import pandas as pd
import scienceplots

plt.style.use("science")

N = 5000

def plots(caso):
    path_x = f"C:/Users/Francisco Rodríguez/Desktop/ordenar/Universidad/Tercer curso/MN2/C++/TRAFICO/Posiciones{caso}.csv"
    path_v = f"C:/Users/Francisco Rodríguez/Desktop/ordenar/Universidad/Tercer curso/MN2/C++/TRAFICO/Velocidades{caso}.csv"
    path_f = f"C:/Users/Francisco Rodríguez/Desktop/ordenar/Universidad/Tercer curso/MN2/C++/TRAFICO/Fluxe{caso}.csv"


    plt.figure(figsize=(6, 4))
    positions = pd.read_csv(path_x)
    speeds = pd.read_csv(path_v)
    flux = pd.read_csv(path_f)

    time = positions["t (s)"]
    x_0 = positions[" x_0 (m)"]
    x_1 = positions[" x_1 (m)"]
    x_2 = positions[" x_2 (m)"]
    x_3 = positions[" x_3 (m)"]
    x_4 = positions[" x_4 (m)"]

    v_0 = speeds[" v_0 (m/s)"]
    v_1 = speeds[" v_1 (m/s)"]
    v_2 = speeds[" v_2 (m/s)"]
    v_3 = speeds[" v_3 (m/s)"]
    v_4 = speeds[" v_4 (m/s)"]
    
    D = flux["D (1/m)"]
    J = flux[" J (1/s)"]

    plt.scatter(time[0], x_0[0], color="red", label="Car 0")
    plt.scatter(time[0], x_1[0], color="orange", label="Car 1")
    plt.scatter(time[0], x_2[0], color="gold", label="Car 2")
    plt.scatter(time[0], x_3[0], color="green", label="Car 3")
    plt.scatter(time[0], x_4[0], color="skyblue", label="Car 4")

    plt.scatter(time[N - 1], x_0[N - 1], color="red")
    plt.scatter(time[N - 1], x_1[N - 1], color="orange")
    plt.scatter(time[N - 1], x_2[N - 1], color="gold")
    plt.scatter(time[N - 1], x_3[N - 1], color="green")
    plt.scatter(time[N - 1], x_4[N - 1], color="skyblue")

    plt.plot(time, x_0, color="red")
    plt.plot(time, x_1, color="orange")
    plt.plot(time, x_2, color="gold")
    plt.plot(time, x_3, color="green")
    plt.plot(time, x_4, color="skyblue")

    plt.xlabel("$t$ (s)")
    plt.ylabel("$x$ (m)")
    plt.legend()
    plt.show()

    plt.figure(figsize=(6, 4))

    plt.scatter(time[0], v_0[0], color="red", label="Car 0")
    plt.scatter(time[0], v_1[0], color="orange", label="Car 1")
    plt.scatter(time[0], v_2[0], color="gold", label="Car 2")
    plt.scatter(time[0], v_3[0], color="green", label="Car 3")
    plt.scatter(time[0], v_4[0], color="skyblue", label="Car 4")

    plt.scatter(time[N - 1], v_0[N - 1], color="red")
    plt.scatter(time[N - 1], v_1[N - 1], color="orange")
    plt.scatter(time[N - 1], v_2[N - 1], color="gold")
    plt.scatter(time[N - 1], v_3[N - 1], color="green")
    plt.scatter(time[N - 1], v_4[N - 1], color="skyblue")

    plt.plot(time, v_0, color="red")
    plt.plot(time, v_1, color="orange")
    plt.plot(time, v_2, color="gold")
    plt.plot(time, v_3, color="green")
    plt.plot(time, v_4, color="skyblue")

    plt.xlabel("$t$ (s)")
    plt.ylabel("$v$ (m/s)")
    plt.legend()
    plt.show()
    
    plt.figure(figsize=(6, 4))

    plt.scatter(D[0], J[0], color="red", label="flux")
    

    plt.scatter(D[N - 1], J[N - 1], color="blue")

    plt.plot(D, J, color="blue")

    plt.xlabel("$D$ (1/m)")
    plt.ylabel("$J$ (1/s)")
    plt.legend()
    plt.show()
    
    

plots(1)

plots(2)

plots(3)