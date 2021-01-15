#--------------------------------2D-DG FEM Implicit Single U ---------------------
import numpy as np
import sympy as sy
import seaborn as sb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

#------------------------------- Parameters -----------------------------------
C= .05                                      # CLF number
#p_i_x = 2                                   # degree of polynomial function +1 in x-dir
#p_i_y = 2                                   # degree of polynomial function +1 in y-dir
Np = 4                                      # number rof points at each element
nt = 1 
N_e_r = 4  
N_e_c = 3
nx =  N_e_c * N_e_r * 2                     # total number of x steps               
ny =  N_e_c * N_e_r * 2                     # total number of y steps                                                                 
N = N_e_c * N_e_r                           # number of elements
L = 0.5                                       # x and y lengths
dx = L/(N_e_r)
dy = L/(N_e_c)


c_x = 0.1                                   # Wave velocity in x-dir
c_y = 0.1                                   # Wave velocity in y-dir
#c = (c_x**2+c_y**2)**0.5                    # speed of the wave
dt = C*dx/c_x


U = np.zeros(N*Np)                          # Wave matrix
Un = np.zeros(N*Np)                         # Dummy variable to save current components of U
U1=U2 = U                                # Dummy matrices to plot 3 time steps
#U_plot = np.ones((3,N*Np))                 # A matrix to save 3 time steps used for plotting the results

U[int(N*Np*.3):int(N*Np*.8)]=1              # Defining wave components

#--------------------------------- Initial Conds ------------------------------
#U[:, 0]= U[:, -1] = 0
#U[0,:] = U[-1,:] = 0

###############################################################################
######################## Natural coordinate system ############################
###############################################################################

#--------------------- Interpolation functions phi(i) -------------------------
# defining variable for symbolic integration of interpolation functions
x, y, a, b, x_m, y_m =sy.symbols('x y a b x_m y_m')

# interpolation functions for the rectangular four noded elements
phi_1 = 1/4*(1-x)*(1-y)
phi_2 = 1/4*(1+x)*(1-y)
phi_3 = 1/4*(1+x)*(1+y)
phi_4 = 1/4*(1-x)*(1+y)


X = x_m + dx/2*x                 # X=global x-coordinate, x_m=mid-point, x=natural x-coordinate
Y = y_m + dy/2*y                 # Y=global y-coordinate, y_m=mid-point, y=natural y-coordinate

J_x = sy.diff(X , x)             # Jacobian function (dX/dx)
J_y = sy.diff(Y, y)              # Jacobian function for (dY/dy) 

interpolation_func = np.array(([phi_1], [phi_2], [phi_3], [phi_4]))  # Matrix form of shape functions
dinterpolation_func_dx = sy.diff(interpolation_func,x)               # differentialsof phi_i functions  with respect to x used in K matrix
dinterpolation_func_dy = sy.diff(interpolation_func,y)               # differentials of phi_i functions with respect to y used in K matrix

#--------- constructing M based in Natural coordinate system -------------------

sub_M = np.zeros((Np,Np))                      # local mass matrix with the size of number of degree of freedom in each element

for i in range(Np):
    for j in range(Np):
        C = interpolation_func[j,0]*interpolation_func[i,0]*J_x*J_y # int(phi_j * phi_i * J_x * J_y * dxdy)
        #print(C,i,j)    
        sub_M[i,j] = sy.integrate(C, (x, -1, 1), (y,-1,1))          # constructing upper diagonal matrix 
        sub_M[j,i] = sub_M[i,j]                                     # filling the lower diagonal entities.

M = np.kron(np.eye(N), sub_M)                                       # Creating global mass matrix
# plt.spy(M)                                    # Useful command to crreat eye matrix with another sum_M repeated on the diagonal

M_inv = np.linalg.inv(M)

# #--------------------------Stifness Matrix 'K' in Equation 44 -----------------
### Constructing sub_K matrix using natural coordinate system

## manually calculated sub_K matrix using local coordinate system 
sub_K=np.zeros((Np,Np))
sub_K=(dt/6)*np.array(([-dy*c_x-dx*c_y , -dy*c_x-dx*c_y/2 , -dy*c_x/2-dx*c_y/2 , -dy*c_x/2-dx*c_y],
                        [dy*c_x-dx*c_y/2 , dy*c_x-dx*c_y , dy*c_x/2-dx*c_y , dy*c_x-dx*c_y/2], 
                        [dy*c_x/2+dx*c_y/2 , dy*c_x/2+dx*c_y , dy*c_x+dx*c_y , dy*c_x+dx*c_y/2],
                        [-dy*c_x/2+dx*c_y , -dy*c_x/2+dx*c_y/2 , -dy*c_x+dx*c_y/2 , -dy*c_x+dx*c_y]))


### constructing global K matrix using kron function
K = np.kron(np.eye(N), sub_K)             # Creating global stiffness matrix
# plt.spy(K)                              

#-------------------------------Flux in Equation 44----------------------------

sub_flux = np.zeros((Np , Np*N_e_r+2))          # initialisation of local flux matrix for one element
sub_flux[0,1] = sub_flux[1,0] = -c_y*dt                         # these 4 lines define fluxes at the boundaries of one element
sub_flux[2, Np*N_e_r] = dt*(c_x+c_y)
sub_flux[0,-7]= sub_flux[3,-6]=-c_x*dt
sub_flux[3, Np*N_e_r+1] = sub_flux[Np-3, Np*N_e_r-1] = c_y*dt

F = np.zeros((Np*N_e_r*N_e_c , Np*N_e_r*N_e_c+Np*N_e_r-2))          # initialisation of global flux matrix
i=0
j=0
for n in range(N_e_c*N_e_r):
    F[i:i+sub_flux.shape[0], j:j+sub_flux.shape[1]]+=sub_flux       # copying local element into the global matrix
    i+= Np
    j+= Np
F=F[:,N_e_r*Np-2:]                              # applying boundary condition in y-dir
i=Np*N_e_r+3

for n in range(N_e_c-1):                        # applying bc in x-dir
    F[i,i-6]=0                                  # in upwind flux left boundary of the first element of each row = 0
    i+= N_e_r*Np                                # movig to the next row

    
    

#plt.spy(F)  

# sub_flux = np.zeros((N_e_r*Np , N_e_r*Np + N_e_r))          # flux matrix for only the one row including boundary points
#                                                         # that is why it is not nXn. later on, it will become nXn.
# i = 0
# # considering only y-dir in one row
# for n in range(N_e_r):
#     sub_flux[i,i]= -c_y*dt
#     sub_flux[i+3*N_e_r,i+4*N_e_r] = c_y*dt
#     i+=1
    
# # considering only x-dir in one row, which is exactly as 1D
# fx = np.zeros((N_e_r*2,N_e_r*2+1))
# i=0
# j=1
# for n in range(N_e_r):
#     fx[i,i]= -c_x*dt
#     fx[j,j+1] = c_x*dt
#     i+=2
#     j+=2
# fx=fx[:,1:]                                             # deleting the boundary points 

# sub_flux[N_e_r:N_e_r+fx.shape[0], N_e_r*2:N_e_r*2+fx.shape[1]]+=fx   # Constructing overal flux matrix for one row only

# F = np.zeros((Np*N_e_r*N_e_c,Np*N_e_r*N_e_c+N_e_r))                  # expanding sub_F to create flux for all rows and columns
# i=0
# j=0
# for n in range(N_e_c):
#     F[i:i+sub_flux.shape[0], j:j+sub_flux.shape[1]]+=sub_flux
#     i+=N_e_r*Np
#     j+=N_e_r*Np
# F = F[:,N_e_r:]    

#--------------------------------RHS Constant in Equation 44-------------------
RHS_cst = (M + K - F)

##-Matrix method----------------------------------------------------------------
##Mrching forward in time
Un=np.zeros(nx)                     # dummy vbl to save current values of U (U^t) 
for n in range(nt):                 # Marching in time
    Un = U.copy()
    RHS = RHS_cst.dot(Un)           # saving U^t to be used at the next timestep calculation
    U=np.linalg.solve(M,RHS)        # solving for U(t+1)
    
    if n==1:                        # saving U at timestep 1 to plot
        U1=U
    elif n==nt//2:                  # saving U at half timestep to plot
        U2=U

# separating U in x and y directions
j=0
U_x = U_y = np.zeros(N_e_r*N_e_c*2)                            # empty array to save for U_y
                            # empty array to save for U_x
for i in U:
    print(i)
    U_y[j]=U[j]
    j+=1
    U_x[j]=U[j]
    j+=1

U_y = np.reshape(U_y, (N_e_c*2, N_e_r))         # reshaping U_y into a N_e_c*2 by N_e_r matrix
U_x = np.reshape(U_x, (N_e_c, N_e_r*2))         # reshaping U_y into a N_e_c by N_e_r*2 matrix
#U = np.zeros((N_e_c*2, N_e_r*2))               # creating the master U matrix in 2D

#------------------------------plot initiation --------------------------------
x = np.linspace(0, L, N_e_c*N_e_r*4)                  # initialising the discretisation of the domain in x-dir
#print(x)
# x = [ele for ele in x for i in range(2)]        # doubling the points for DG
# x=x[1:-1]                                       # excluding boundary points
# #print(x)
#plt.plot(x, U_x[0:7], label='Timestep 1')

y = np.linspace(0, L, N_e_c+1)                  # initialising the discretisation of the domain in y-dir
# #print(y)
# y = [ele for ele in y for i in range(2)]        # doubling the points for DG
# y=y[1:-1]                                       # excluding boundary points
# #print(y)
plt.plot(x,U)
X, Y =np.meshgrid(x,y)                          # Creating a mesh grid
#plt.plot(X,Y,U)


# plt.figure()        
# ax = plt.gca(projection='3d')
# ax.plot_surface(X, Y, U1 , label='t=0')
# ax.plot_surface(X, Y, U2, label='t=nt/2')
# ax.plot_surface(X, Y, U3, label='t=final')
# ax.set_ylabel('$y$')
# ax.set_xlabel('$x$')
# ax.set_zlabel('$U$')
# plt.legend()
# plt.show()

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(X, Y, U[:], cmap=plt.cm.viridis)
# #ax.plot(x, y, U, label='Square Wave')
# plt.show()













###############################################################################
########################## Local coordinate system ############################
###############################################################################

#----------------------- Interpolation functions phi(i) ----------------------
# x=sy.Symbol('x')
# y =sy.Symbol('y')
# phi_1 = (1-x/dx)*(1-y/dy)
# phi_2 = (x/dx)*(1-y/dy)
# phi_3 = (x*y)/(dx*dy)
# phi_4 = (y/dy)*(1-x/dx)

# interpolation_func = np.array(([phi_1], [phi_2], [phi_3], [phi_4]))
# dinterpolation_func_dx = sy.diff(interpolation_func,x)
# dinterpolation_func_dy = sy.diff(interpolation_func,y)

# #----------Mass Matrix 'M' in Equation 44 in consistent form ------------------
# t1 = time.time()                            # starting for timing the M_diag_inv calculation
# sub_M = np.zeros((Np,Np))                 # local mass matrix

# for i in range(Np):
#     for j in range(i, Np):
#         C = interpolation_func[i,0]*interpolation_func[j,0]
#         #print(sy.integrate(C, (x, 0, a), (y,0,b)))
#         sub_M[i,j] = sy.integrate(C, (x, 0, dx), (y,0,dy)) 
#         sub_M[j,i] = sub_M[i,j]
# M=np.zeros((N*Np,N*Np))                         # generating global mass matrix
# i=0
# while i<N*Np:
#     M[i:i+Np, i:i+Np]= sub_M[0:Np,0:Np]
#     i+=Np
# t2 = time.time()                            # end point of M_diag_inv generation
# print(str(t2-t1))
        
# #--------------------------Stifness Matrix 'K' in Equation 44 -----------------

# sub_K=np.zeros((Np,Np))                   # local stifness matrix
# K = np.zeros((N*Np,N*Np))                 # generating global stifness matix
# for i in range(Np):
#     for j in range(Np):
#         C_x = dt*c_x*interpolation_func[i,0]*(dinterpolation_func_dx[j,0])
#         C_y = dt*c_y*interpolation_func[i,0]*(dinterpolation_func_dy[j,0])
#         sub_K[i,j] = sy.integrate(C_x, (x, 0, dx), (y,0,dy)) + sy.integrate(C_y, (x, 0, dx), (y,0,dy))

# i=0
# while i<N*Np:
#     K[i:i+Np, i:i+Np]= sub_K[0:Np,0:Np]
#     i+=Np 
    
#------------------Exact area and define integral with python (YouTube)--------

# # x=np.linspace(-1,3,1000)
# # def f(x): return x**2
# # plt.plot(x,f(x))
# # plt.axhline(color= 'black')
# # plt.fill_between(x, f(x), where=[x>0 and x<2 for x in x], color='green', alpha = 0.3)
# # plt.show()
# # x=sy.Symbol('x')
# # # print(sy.integrate(f(x), (x)))
# M=np.zeros((4,4))
 
# # c_x = 1
# # c_y = 1

# a =1
# b =1
# x =sy.Symbol('x')
# y =sy.Symbol('y')

# phi_1 = (1-x/a)*(1-y/b)
# phi_2 = (x/a)*(1-y/b)
# phi_3 = (x*y)/(a*b)
# phi_4 = (y/b)*(1-x/a)

# A = np.array(([phi_1], [phi_2], [phi_3], [phi_4]))
# # dA_dx = sy.diff(A,x)
# # dA_dy = sy.diff(A,y)

# #A_difff=sy.integrate(A_by_dA_x, (x, 0, a), (y,0,b))
# for i in range(4):
#     for j in range(4):
#         C = A[i,0]*A[j,0]
#         #print(sy.integrate(C, (x, 0, a), (y,0,b)))
#         M[i,j] = sy.integrate(C, (x, 0, a), (y,0,b))

# x, y, dx, dy =sy.symbols('x y dx dy')

# phi_1 = 1/4*(1-x)*(1-y)
# phi_2 = 1/4*(1+x)*(1-y)
# phi_3 = 1/4*(1+x)*(1+y)
# phi_4 = 1/4*(1-x)*(1+y)

# interpolation_func = np.array(([phi_1], [phi_2], [phi_3], [phi_4]))
# dinterpolation_func_dx = sy.diff(interpolation_func,x)
# dinterpolation_func_dy = sy.diff(interpolation_func,y)

# sub_K2=np.zeros((Np,Np))                   # local stifness matrix
                       
# for i in range(Np):
#     for j in range(Np):
#         C_x = dt*c_x*interpolation_func[i,0]*(dinterpolation_func_dx[j,0])
#         C_y = dt*c_y*interpolation_func[i,0]*(dinterpolation_func_dy[j,0])
#         sub_K2[i,j] = sy.integrate(C_x, (x, 0, dx), (y,0,dy)) + sy.integrate(C_y, (x, 0, dx), (y,0,dy))

# K = np.zeros((N*Np,N*Np))   # generating global stifness matix
# i=0
# while i<N*Np:
#     K[i:i+Np, i:i+Np]= sub_K[0:Np,0:Np]
#     i+=Np