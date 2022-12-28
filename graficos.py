import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import pandas as pd
plt.style.use("science")

N = 2000
path = "C:/Users/minih/Documents/GitHub/ModeloTrafico/cmake-build-debug/Instantaneo.csv"

plt.figure(figsize=(6, 4))
df = pd.read_csv(path)
time = df["t (s)"]
x_0 = df[" x_0 (m)"]
x_1 = df[" x_1 (m)"]
x_2 = df[" x_2 (m)"]
x_3 = df[" x_3 (m)"]
x_4 = df[" x_4 (m)"]

plt.scatter(time[0], x_0[0], color="red")
plt.scatter(time[0], x_1[0], color="orange")
plt.scatter(time[0], x_2[0], color="gold")
plt.scatter(time[0], x_3[0], color="green")
plt.scatter(time[0], x_4[0], color="skyblue")

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

plt.show()

"""fig, axes = plt.subplots()
graficar0, = plt.scatter([], [])
graficar1, = plt.scatter([], [])
graficar2, = plt.scatter([], [])
graficar3, = plt.scatter([], [])
graficar4, = plt.scatter([], [])

Re = 0
Im = 1

def funcion():
    graficar0.set_data(time, x_0)
    graficar1.set_data(time, x_1)
    graficar2.set_data(time, x_2)
    graficar3.set_data(time, x_3)
    graficar4.set_data(time, x_4)

    return graficar0, graficar1, graficar3, graficar2, graficar4


animation.FuncAnimation(fig, funcion, blit=True)"""
"""# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)


# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    line.set_data([], [])
    line.set_data([], [])
    line.set_data([], [])
    line.set_data([], [])

    return line,


# animation function.  This is called sequentially
def animate(i):
    line.set_data(time[i], x_0[i])
    line1.set_data(time[i], x_1[i])
    line2.set_data(time[i], x_2[i])
    line.set_data(time[i], x_3[i])
    line.set_data(time[i], x_4[i])

    return line,


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('basic_animation.mp4', fps=4, extra_args=['-vcodec', 'libx264'])

plt.show()
"""
