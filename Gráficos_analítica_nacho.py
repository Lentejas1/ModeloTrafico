import numpy as np
from matplotlib import pyplot as plt

d=16 #distancia de seguimiento en la cual los coches empiezan a aminorar la velocidad
L=4
v_max=120/3.6
rho_max=1/L  #densidad máxima
rho_cri=1/(L+d)  #densidad en la cual empiezan a bajar la velocidad

def v(rho):
    return v_max*np.log(rho_max/rho)*np.log(rho_max/rho_cri)**(-1)  #fórmula de los apuntes

h=2000
rho=np.linspace(rho_cri,rho_max,h)

rho_op=rho[0]  #velocidad que hace el flujo máximo
v_op=v(rho_op)  #densidad que hace el flujo máximo

for i in range(h):
    if rho[i]*v(rho[i])>=rho_op*v_op:
        rho_op=rho[i]
        v_op=v(rho_op)
        
print("La velocitat óptima és", v_op)


plt.plot(rho/rho_cri, v(rho)/v_max,"black") #esta es la de velocidad en función de densidad
plt.plot([0,1],[1,1],"black")
plt.plot(rho/rho_cri, (v(rho)*rho)/(v_op*rho_op), "blue")
plt.plot([0,1],[0,0.88],"blue")
plt.plot([rho_op/rho_cri, rho_op/rho_cri], [0,3],"--")
plt.ylim(0,1.1)
plt.plot(rho/rho_cri,rho*v(rho)/rho_cri*v(rho_cri), "--")
plt.plot([1,1],[0,3],"--")
plt.title("Velocitat y flux en funció de la densitat per la solució analítica")
plt.xlabel(r"$\rho$"" /"r"$\rho$_cri")
plt.ylabel("v("r"$\rho$)""/""v_max""  y  J(v)/J_max")


