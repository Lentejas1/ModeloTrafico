import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

path = "C:/Users/minih/Documents/GitHub/ModeloTrafico/cmake-build-debug/Instantaneo.csv"

df = pd.read_csv(path)

time = df["t (s)"]
x_0 = df["x_0 (m)"]

plt.scatter(time, x_0)
print(df)

