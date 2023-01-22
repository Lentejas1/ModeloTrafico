# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 21:15:40 2023

@author: Francisco Rodríguez
"""

import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use("science")


d=26 #distancia de seguimiento en la cual los coches empiezan a aminorar la velocidad
L=4
vmax=120/3.6
rhomax=1/L  #densidad máxima
rhocri=1/(L+d)  #densidad en la cual empiezan a bajar la velocidad

def v(rho):
    return vmax*np.log(rhomax/rho)*np.log(rhomax/rhocri)**(-1)  #fórmula de los apuntes

h=2000
rho=np.linspace(rhocri,rhomax,h)

rhoop=rho[0]  #velocidad que hace el flujo máximo
vop=v(rhoop)  #densidad que hace el flujo máximo

for i in range(h):
    if rho[i]*v(rho[i])>=rhoop*vop:
        rhoop=rho[i]
        vop=v(rhoop)
        
print("La velocitat óptima és", vop, rhoop, vop*rhoop)

plt.figure(figsize=(6, 4))
plt.plot(rho/rhocri, v(rho)/vmax,"black", label= r"$J(\rho)$") #esta es la de velocidad en función de densidad
plt.plot([0,1],[1,1],"black")
plt.plot(rho/rhocri, (v(rho)*rho)/(vop*rhoop), "blue", label= r"$v(\rho)$")
plt.plot([0,1],[0,0.73],"blue")
plt.plot([rhoop/rhocri, rhoop/rhocri], [0,3],"--", label= r"$\rho_{opt}$")
plt.ylim(0,1.1)
plt.plot(rho/rhocri,rho*v(rho)/rhocri*v(rhocri), "--", label= r"$\rho_{max}$")
plt.plot([1,1],[0,3],"--", label= r"$\rho_{crit}$")
plt.xlabel(r"$\rho/\rho_{crit}$")
plt.ylabel(r"$v(\rho)/v_{max}$" " " "y" " " r"$J(\rho)/J_{max}$")
plt.legend(loc='upper right')
plt.show()