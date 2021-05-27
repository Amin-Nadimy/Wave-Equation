import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import seaborn as sb

#------------------------------- 2-D FDM------------------------------------------
C= .05                           # CLF number
nt = 200                        # Number of time steps
nx = 100                         # Number of steps in x-direction
ny = 100                        # Number of steps in y-direction
L = .5                            # Lengths of the domain in x and y directions
dx = L/(nx-1)                    # dispacement legth at each time-step in x-direction
dy = L/(ny-1)                    # dispacement legth at each time-step in y-direction
d_displacement = (dx**2+dy**2)**0.5 # overall displacement of the wave in both directions
x=np.linspace(0, L, nx)
y=np.linspace(0, L, ny)
c_x = 0.1                        # Wave velocity in x-dir
c_y = 0.1                        # Wave velocity in y-di
c = (c_x**2+c_y**2)**0.5            # overall velocity of the wave
dt = C*d_displacement/c          # length of time step. 
U = np.zeros((ny,nx))            # Wave matrix
Un = np.zeros((ny,nx))           # Dummy variable to save current components of U
U1=U2=U3 = U                     # Dummy matrices to plot 3 time steps

U[int(ny*.1):int(ny*.3), int(nx*0.1):int(nx*0.3)]=1     # Defining wave components
X, Y =np.meshgrid(x,y)           # Creating a mesh grid

#--------------------------------- Initial Conds ------------------------------
U[:, 0]= U[:, -1] = 0
U[0,:] = U[-1,:] = 0

#------------------------------------------------------------------------------
for n in range(nt+1):             # marching in time steps
    Un=U.copy()                   # dummy vbl to save current solutions
    U[1:,1:] = Un[1:,1:] - c_y*(dt/dy)*(Un[1:,1:] - Un[:-1,1:]) - c_x*(dt/dx)*(Un[1:,1:] - Un[1:, :-1])

    if n==0:                    # saving the first time step to plot
        U1 = Un
    elif n==int(nt/2):          # saving the half time step to plot
        U2 = Un
    elif n==int(nt*0.99):       # saving the final time step to plot
        U3 = Un

#------------------------------ Plotting the results ----------------------------        
#fetch = open("2D_FDM.txt", 'w')
#fetch.write(2)
#fetch.close()

plt.figure(1)        
ax = plt.gca(projection='3d')
ax.plot_surface(X, Y, U1 , label='first timestep')
ax.plot_surface(X, Y, U2, label='second timestep')
ax.plot_surface(X, Y, U3, label='final timestep')
plt.legend()
ax.set_ylabel('$y$')
ax.set_xlabel('$x$')
ax.set_zlabel('$U$')

plt.show()

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(X, Y, U[:], cmap=plt.cm.viridis)
##ax.plot(x, y, U, label='Square Wave')
#plt.show()
# import pygmsh
# def test():
#     with pygmsh.geo.Geometry() as geom:
#         geom.add_rectangle(0.0, 1.0, 0.0, 1.0, 0.0, 0.1)
#         mesh = geom.generate_mesh()


# z= np.zeros((10,15))
# z[5:8,3:10]=1
# z[9,12]=5
# plt.figure()
# ax = sb.heatmap(z)
# ax.invert_yaxis()
# plt.show()
######################################## 1-D  #################################
#C = .6
#dt = 1
#dx = 1
#nx = 100
#nt = 200
#
#u = np.zeros(nx)
#u[int(nx*0.05):int(nx*0.2)]=1
#u[0] = u[nx-1] = 0
#
#
#for n in range(nt):
#    plt.clf()
#    un = u.copy()
#    for i in range(nx):
#        u[i] = un[i] - C*dt/dx*(un[i]-un[i-1])
#    plt.axis([0,1, 0,1])
#    plt.plot(np.linspace(0,2,nx),u)
#    plt.pause(0.0001)

#================================ Sin (x) =====================================       

# C  = 1
# c= .1
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
#          # Forward euler's method for time and backward for space
#         u[i] = un[i] - c*(dt/dx)*(un[i]-un[i-1])

#     plt.axis([0,1, 0,1.2])
#     plt.plot(x,u)
#     plt.pause(.0001)




































       