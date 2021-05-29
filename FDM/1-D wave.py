import numpy as np
import matplotlib.pyplot as plt

# C  = 1
# c= -.1
# L=2
# nx = 500                    # distance which the wave travels
# nt = 40                     # time length
# dx = L/(nx-1)
# dt = .04#C*dx/c


# u = np.zeros(nx)
# u[5:30]=1                 # It's just an example IC which I chose
# x  = np.arange(0, nx)*dx#x=np.linspace(0,20,nx)
# plt.plot(x,u)

# #Initial Conditions:

# u[nx-1]= 0

# for n in range(nt):
#     u[0] = np.cos(n*dt)
#     plt.clf()
#     un = u.copy()
#     for i in range(1, nx-1):
#           # Forward euler's method for time and backward for space
#         u[i] = un[i] - c*(dt/dx)*(un[i]-un[i-1])

#     plt.axis([0,1, 0,1.2])
#     plt.plot(x,u)
#     plt.pause(.0001)

#x=np.linspace(0,20,100)
#y=np.cos(x)
#plt.plot(x,y)

###----------------------------------------
# nt,   nx, tmax, xmax, c
# 151, 1001, 0.5, 2.0, 0.5
#C= .1
#nt = 301
#nx = 301
#tmax = 0.5
#L = 2
#dx = L/(nx-1)
#c = .5
#dt = C*dx/c
#
#u = np.zeros((nx,nt))
#
# #boundary conditions
#u[0,:]=u[nx-1,:]=0
#
# #initial conditions
#for i in range(1,nx-1):
#    if i > (nx-1)/4 and i<(nx-1)/2:
#        u[i,0]=1
#    else:
#        u[i,0]=0
#
# #Loop
#for n in range(0,nt-1):
#  
#    for i in range(nx-1):
#        u[i,n+1] = u[i,n]-c*(dt/dx)*(u[i,n]-u[i-1,n])
#plt.figure(1)
#plt.xlabel('Distrance')
#plt.ylabel('u')
#plt.plot(np.linspace(0,L,nx),u)

#-------------------------------------------
# This is the Python code for Finite Difference Mehtod in 1D
C= 0.05                         # Courant number
nt = 250                       # total number of time-steps 
nx = 60                        # total number of displacement 
L = 0.5                           # length of thedomain in 1D
dx = L/(nx-1)                   # dispacement legth each time-step                   
c = 0.1                          # velocity of the wave
dt = C*dx/c                     # length of each time-step

U = np.zeros(nx)                # initialising U matrix
U[nx//5:nx//2]=1                     # defining square wave
selected_elements=np.ones((3,nx)) # dummy vbl used for saving solutions in particular time-step

for n in range(nt):             # marchig time
    Un=U.copy()                 # dummy vbl to save current values of U
    for i in range(nx):         # marching in space
        U[i] = Un[i]-c*(dt/dx)*(Un[i]-Un[i-1])
    if n==1:                    # savin the solution at the begining of the simulation
        selected_elements[0,:] = U.copy()
    if n==nt//2:                    # savin the solution at the middle of the simulation
        selected_elements[1,:] = U.copy()
    if n==nt-1:                    # savin the solution at the end of the simulation
        selected_elements[2,:] = U.copy()

# ----------------------- Plotting the results ------------------------------------
plt.figure(1)
plt.plot(np.linspace(0,L,nx), selected_elements[0,:], label='first timestep')
plt.plot(np.linspace(0,L,nx), selected_elements[1,:], label='mid-timestep')
plt.plot(np.linspace(0,L,nx), selected_elements[2,:], label='last timestep')
plt.xlabel('Distrance')
plt.ylabel('u')
plt.legend()
#----------------------------------------------------------
#
#C= .1
#nt = 301
#nx = 301
#tmax = 0.5
#L = 2
#dx = L/(nx-1)
#c = .5
#dt = C*dx/c
#range(nx)
#u = np.zeros(nx)
#u[25:100]=1
#selected_data=np.ones((3,nx))
#u[0] = 0
#u[nx-1]= 0
#
#
##Loop
#for n in range(nt):
#    un=u.copy()
#    for i in range(1,nx-1):
#        u[i] = un[i]-c*(dt/dx)*(un[i]-un[i-1])
#    if n==1:
#        selected_data[0,:] = u.copy()
#    if n==150:
#        selected_data[1,:] = u.copy()
#    if n==300:
#        selected_data[2,:] = u.copy()
#    # plt.plot(np.linspace(0,L,nx),u)
#plt.figure(1)
#plt.xlabel('Distrance')
#plt.ylabel('u')
#plt.plot(np.linspace(0,L,nx), selected_data[0,:])
#plt.plot(np.linspace(0,L,nx), selected_data[1,:])
#plt.plot(np.linspace(0,L,nx), selected_data[2,:])



        
        

        
        



