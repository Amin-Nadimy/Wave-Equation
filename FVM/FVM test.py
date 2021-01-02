import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sin, cos
#
C= .1
nt = 301
nx = 301
tmax = 0.5
L = 2
dx = L/(nx-1)
c = .1
dt = C*dx/c
range(nx)
u = np.zeros(nx)
u[25:100]=1
selected_data=np.ones((3,nx))
u[0] = np.sin(nt)
u[nx-1]= 0


#Loop
for n in range(nt):
    un=u.copy()
    for i in range(1,nx):
        u[i] = un[i]-c*(dt/dx)*(un[i]-un[i-1])
    if n==1:
        selected_data[0,:] = u.copy()
    if n==150:
        selected_data[1,:] = u.copy()
    if n==300:
        selected_data[2,:] = u.copy()
    # plt.plot(np.linspace(0,L,nx),u)
plt.figure(1)
plt.xlabel('Distrance')
plt.ylabel('u')
plt.plot(np.linspace(0,L,nx), selected_data[0,:])
plt.plot(np.linspace(0,L,nx), selected_data[1,:])
plt.plot(np.linspace(0,L,nx), selected_data[2,:])

#x=np.linspace(0,2,100)
#y=np.sin(x)
#plt.plot(x,y)
