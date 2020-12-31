#--------------------------------2D-DG FEM Implicit Single U ---------------------
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

#------------------------------- Parameters -----------------------------------
C= .005                                      # CLF number
p_i_x = 2                                   # degree of polynomial function +1 in x-dir
p_i_y = 2                                   # degree of polynomial function +1 in y-dir
Np = 4                                      # numbe rof points in each element
nt = 1000                                   # Number of time steps
nx = 21                                     # Number of x steps
ny = 6                                      # Number of y steps
N = (nx-1)*(ny-1)                           # number of elements
L = 2                                       # x and y lengths
dx = L/(nx-1)
dy = L/(ny-1)


c_x = 0.1                                   # Wave velocity in x-dir
c_y = 0.1                                   # Wave velocity in y-dir
c = (c_x**2+c_y**2)**0.5                    # speed of the wave
dt = C*dx/c
U = np.zeros(N*2*Np)                   # Wave matrix
Un = np.zeros(N*2*Np)                  # Dummy variable to save current components of U
U1=U2=U3 = U                                # Dummy matrices to plot 3 time steps
#U_plot = np.ones((3,N*Np))    # A matrix to save 3 time steps used for plotting the results

U[int(N*2*Np*.3):int(N*2*Np*.8)]=1              # Defining wave components

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

X = x_m + dx/2*x                # X=global x-coordinate, x_m=mid-point, x=natural x-coordinate
Y = y_m + dy/2*y                # Y=global y-coordinate, y_m=mid-point, y=natural y-coordinate

dX = sy.diff(X , x)             # Jacobian function (dX/dx)
dY = sy.diff(Y, y)              # Jacobian function for (dY/dy) 

interpolation_func = np.array(([phi_1], [phi_2], [phi_3], [phi_4]))  # Matrix form of shape functions
dinterpolation_func_dx = sy.diff(interpolation_func,x)               # differentialsof phi_i functions  with respect to x used in K matrix
dinterpolation_func_dy = sy.diff(interpolation_func,y)               # differentials of phi_i functions with respect to y used in K matrix

#--------- constructing M based in Natural coordinate system -------------------
t1 = time.time() 
sub_M = np.zeros((Np,Np))                      # local mass matrix with the size of number of degree of freedom in each element

for i in range(Np):
    for j in range(i, Np):
        C = interpolation_func[j,0]*interpolation_func[i,0]* 1/dX * 1/dY # int(phi_j * phi_i * J_x * J_y * dxdy)
        #print(sy.integrate(C, (x, -1, 1), (y,-1,1)))
        sub_M[i,j] = sy.integrate(C, (x, -1, 1), (y,-1,1))          # constructing upper diagonal matrix 
        sub_M[j,i] = sub_M[i,j]                                     # filling the lower diagonal entities.


sub_M2 =  np.zeros((2*Np,2*Np))   
for i in range(Np):
    for j in range(Np):
        sub_M2[i*2, j*2] = sub_M[i,j]

for i in range(Np):
    for j in range(Np):
        sub_M2[i*2+1, j*2+1] = sub_M[i,j]

M = np.kron(np.eye(N), sub_M2)             # Creating global mass matrix
# plt.spy(M)                                    # Useful command to crreat eye matrix with another sum_M repeated on the diagonal


t2 = time.time()                            # end point of M_diag_inv generation
print(str(t2-t1))

# #--------------------------Stifness Matrix 'K' in Equation 44 -----------------

# sub_K=np.zeros((Np,Np))                        # local stifness matrix
                     
# for i in range(Np):                            # Creating general local stiffness matrix
#     for j in range(Np):
#         C_x = dt*c_x*interpolation_func[i,0]*(dinterpolation_func_dx[j,0])*1/dX*dX*dY
#         C_y = dt*c_y*interpolation_func[i,0]*(dinterpolation_func_dy[j,0])*1/dY*dX*dY
#         sub_K[i,j] = sy.integrate(C_x, (x, -1, 1), (y,-1,1)) + sy.integrate(C_y, (x, -1, 1), (y,-1,1))

# sub_K2 =  np.zeros((2*Np,2*Np))                # Expanding general local stiffness matrix
# for i in range(Np):
#     for j in range(Np):
#         sub_K2[i*2, j*2] = sub_K[i,j]

# for i in range(Np):
#     for j in range(Np):
#         sub_K2[i*2+1, j*2+1] = sub_K[i,j]

# K = np.kron(np.eye(N*2*Np), sub_K2)             # Creating global stiffness matrix
# plt.spy(K)                                    # Useful command to crreat eye matrix with another sum_M repeated on the diagonal

sub_K=(dt/6)*np.array(([-dy*c_x,0,-dy*c_x,0,-dy*c_x/2,0,-dy*c_x/2,0],[0, -dx*c_y,0,-dx*c_y/2,0,-dx*c_y/2,0,-dx*c_y],
                       [dy*c_x,0,dy*c_x,0,dy*c_x/2,0,dy*c_x,0], [0, -dx*c_y/2,0,-dx*c_y,0,-dx*c_y,0,-dx*c_y/2],
                       [dy*c_x/2,0,dy*c_x/2,0,dy*c_x,0,dy*c_x,0],[0, dx*c_y/2,0,dx*c_y,0,dx*c_y,0,dx*c_y/2],
                       [-dy*c_x/2,0,-dy*c_x/2,0,-dy*c_x,0,-dy*c_x,0],[0, dx*c_y,0,dx*c_y/2,0,dx*c_y/2,0,dx*c_y]))                        # local stifness matrix

K = np.kron(np.eye(N), sub_K)             # Creating global stiffness matrix
# plt.spy(K)                                    # Useful command to crreat eye matrix with another sum_M repeated on the diagonal
    
###############################################################################
###############################################################################
###############################################################################

#-------------------------------Flux in Equation 44----------------------------
sub_F = dt*np.array([[-c_x,0,0,0,0,0,0,0], 
                     [0,-c_y,0,0,0,0,0,0], 
                     [0,0,c_x,0,0,0,0,0 ], 
                     [0,0,0,-c_y,0,0,0,0],
                     [0,0,0,0, c_x,0,0,0],
                     [0,0,0,0,0, c_y,0,0],
                     [0,0,0,0,0,0,-c_x,0],
                     [0,0,0,0,0,0,0, c_y]])

F = np.kron(np.eye(N), sub_F)             # Creating global flux matrix
# plt.spy(K)                                    # Useful command to crreat eye matrix with another sum_M repeated on the diagonal

#--------------------------------RHS Constant in Equation 44-------------------
RHS_cst = (M + K - F)

##-Matrix method----------------------------------------------------------------
##Mrching forward in time
Un=np.zeros(nx)                     # dummy vbl to save current values of U (U^t) 
t3 = time.time()
for n in range(nt):                 # Marching in time
    Un = U.copy()
    RHS = RHS_cst.dot(Un)           # saving U^t to be used in the next time step calculation
    U=np.linalg.solve(M,RHS)    
    
    if n==1:
        U1 = U.copy()                     # saving U(t=1)
    if n==int(nt/2):
        U2 = U.copy()                    # saving U(t=nt/2)
    if n==int(nt*0.99):
        U3 = U.copy()                    # saving U(t= almost the end to time steps)
t4 = time.time()
#------------------------------plot initiation --------------------------------
x=np.linspace(0,2,N*Np)
plt.figure(1)
plt.axis([0,L, -1,2])
plt.plot(x, U1, label='Timestep 1')
plt.plot(x, U2, label='Timestep 0.5*nx')
plt.plot(x, U3, label='Timestep 0.9*nx')
plt.xlabel('Distrance')
plt.ylabel('U')
plt.legend()
plt.title(f'Simulation Duration: {round((t4-t3)/60, 2)} minutes')















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