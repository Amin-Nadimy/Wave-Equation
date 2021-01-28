import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sympy as sy
from scipy.special.orthogonal import p_roots


points, weights = np.polynomial.legendre.leggauss(3)


def gauss(f, a, b, points, weights):
    x = np.zeros(len(points))
    for i in range(len(points)):
        x[i] = (b+a)/2 + (b-a)/2 *points[i]
    answer = 0
    for i in range(len(points)):
        answer = answer + (b-a)/2 * (weights[i]*f(x[i]))

    return answer#(b-a)/2 * (weights[0]*f(x[0]) + weights[1]*f(x[1]) + weights[2]*f(x[2]))   

C = 0.05
c = 0.1
L = 0.3
N = 400
nx = 2*N
dx = L/(nx-1)
dt = C*dx/c
nt = 2000
nt2 = 1000

h_e = L/N           # length of element

################### Mass matrix M1 = int ( phi_i phi_j ) ######################
M1_00 = lambda x: 1/2*(1-x) *1/2 *(1-x) *h_e/2
M1_01 = M1_10 = lambda x: 1/2 * (1-x) * 1/2 * (1+x)* h_e/2
M1_11 = lambda x: 1/4 * (1+x)**2* h_e/2

sub_M1 = np.zeros((2,2))
sub_M1[0,0] = gauss(M1_00, -1, 1, points, weights)
sub_M1[0,1] = sub_M1[1,0] = gauss(M1_01, -1, 1, points, weights)
sub_M1[1,1] = gauss(M1_11, -1, 1, points, weights)

#################### K1 mweightstrix = int(c phi_i d_phi_j/dx  ) ####################
K_00 = lambda x: dt*c*1/2*(1-x) *-1/2 
K_01 = lambda x: dt*c*1/2*(1-x) *-1/2 
K_10 = lambda x: dt*c*1/2*(1+x) *1/2
K_11 = lambda x: dt*c*1/2*(1+x) *1/2 

sub_K1 = np.zeros((2,2))
sub_K1[0,0] = gauss(K_00, -1, 1, points, weights)
sub_K1[0,1] = gauss(K_01, -1, 1, points, weights)
sub_K1[1,0] = gauss(K_10, -1, 1, points, weights)
sub_K1[1,1] = gauss(K_11, -1, 1, points, weights)

############################# K2 = int ( nu d_phi_j/dx ) ######################
U = np.zeros((2))
U[:] = 0
#U[1,0] = 0
Un = np.zeros((2))
Un[0] = Un[1] = 1
#Un[1:3,0] = 1

a = 1
#for n in range(3):
for i in range(3):
    j_1 = lambda x: dt* a* -1/2* c * h_e/(8*c) * ( U[0]/dt*(1-x)/2 -Un[0]/dt  * ( (1-x)/2 + c/h_e) )**2 /  (  (U[0]-Un[0]/dt  * (1-x)/2  )**2  +  (-Un[0]/h_e)**2 )
    j_2 = lambda x: dt* a*  1/2* c * h_e/(8*c) * ( U[1]/dt*(1+x)/2 -Un[1]/dt  * ( (1+x)/2 - c/h_e) )**2 /  (  (U[1]-Un[1]/dt  * (1+x)/2  )**2  +  (Un[1]/h_e)**2 )

    sub_K2 = np.zeros((2,2))
    sub_K2[0,0] = sub_K2[1,0] = gauss(j_1, -1, 1, points, weights) 
    sub_K2[0,1] = sub_K2[1,1] = gauss(j_2, -1, 1, points, weights)
    
    sub_K = sub_K2 - sub_K1


## if Un = 0 gives nan value for U due to denominator becomes zero. IF we want to calculate for all U.
# nx = 4
# U = np.zeros((nx))
# U[:] = 0
# #U[1,0] = 0
# Un = np.zeros((nx))
# Un[1] = Un[2] = 1
# #Un[1:3,0] = 1

# sub_K2 = np.zeros((nx,nx))
# #for n in range(3):
# j=0
# for i in range(3):
#     while j < nx:
#         j_1 = lambda x: -1/2* dx/(8*c) * ( U[j]/dx*(1-x)/2 -Un[j]/dt  * ( (1-x)/2 - c/dx) )**2 /  (  (U[j]-Un[j]/dt  * (1-x)/2  )**2  +  (-Un[j]/dx)**2 )
#         j_2 = lambda x:  1/2* dx/(8*c) * ( U[j+1]/dx*(1+x)/2 -Un[j+1]/dt  * ( (1+x)/2 + c/dx) )**2 /  (  (U[j+1]-Un[j+1]/dt  * (1+x)/2  )**2  +  (Un[j+1]/dx)**2 )
# #        i+=2    
#         sub_K2[j,j] = sub_K2[j+1,j] = gauss(j_1, -1, 1, E, A) *dt
#         sub_K2[j,j+1] = sub_K2[j+1,j+1] = gauss(j_2, -1, 1, E, A) *dt
#         j+=2
        
#         print(j)
    
    
    
    ############################### M2 = int ( nu phi_j  ) ########################
    
    j_1 = lambda x: a* h_e/2 * (1-x)/2 * h_e/(8*c) * ( U[0]/dt*(1-x)/2 -Un[0]/dt  * ( (1-x)/2 - c/h_e) )**2 /  (  (U[0]-Un[0]/dt  * (1-x)/2  )**2  +  (-Un[0]/h_e)**2 )
    j_2 = lambda x: a* h_e/2 * (1+x)/2 * h_e/(8*c) * ( U[1]/dt*(1-x)/2 -Un[1]/dt  * ( (1-x)/2 + c/h_e) )**2 /  (  (U[1]-Un[1]/dt  * (1-x)/2  )**2  +  (-Un[1]/h_e)**2 )
    
    sub_M2 = np.zeros((2,2))
    sub_M2[0,0] = sub_M2[1,0]=gauss(j_1, -1, 1, points, weights)
    sub_M2[0,1] = sub_M2[1,1]=gauss(j_2, -1, 1, points, weights)
    
    sub_M = sub_M1 + sub_M2
    
    ##############################################################################
    sub_F = np.zeros((2,2))
    sub_F[1,1] = c*dt
    
    ##############################################################################
    rhs_cst = (sub_M - sub_K -sub_F)
    RHS = rhs_cst.dot(Un)
    
    U = np.linalg.solve(sub_M,RHS)


K = np.kron(np.eye(N), sub_K)
M = np.kron(np.eye(N), sub_M)
F = np.zeros((nx,nx+1))
i=0
j=1
while i<=nx-1:
    F[i,i]= c*dt
    F[j,j+1] = -c*dt
    i+=2
    j+=2
F=F[:,1:]                            # excludig left boundary to get nx by nx matrix

RHS_cst = (M - K + F)

U = np.zeros(nx)            # U is a square wave between 0 <U< 1
U[0] = U[nx-1] = 0          # Dirichlet BC
U[int(L*nx*0.3):int(L*nx*1.7)]=1
Un=np.zeros(nx)                     # dummy vbl to save current values of U (U^t)
U_plot = np.ones((3,nx))    # A matrix to save 3 time steps used for plotting the results 

for n in range(nt):                 # Marching in time
    Un = U.copy()
    RHS = RHS_cst.dot(Un)           # saving U^t to be used in the next time step calculation
    U=np.linalg.solve(M,RHS) 
    
    if n==1:
        U_plot[0,:] = U.copy()      # saving U(t=1)
    if n==int(nt/2):
        U_plot[1,:] = U.copy()      # saving U(t=nt/2)
    if n==int(nt*0.99):
        U_plot[2,:] = U.copy()      # saving U(t= almost the end to time steps)








U3 = np.zeros(nx)            # U is a square wave between 0 <U< 1
U_plot3 = np.ones((3,nx))    # A matrix to save 3 time steps used for plotting the results

#------------------------------------------------------------------------------
# Boundary Conditions
U3[0] = U3[nx-1] = 0          # Dirichlet BC

#------------------------------------------------------------------------------
# Initial conditions
U3[int(L*nx*0.3):int(L*nx*1.7)]=1

#------------------------------- Interpolation functions ----------------------
x_bar=sy.Symbol('x_bar')                     # defining local x symbol for creating interpolation functions and their integrations
phi_1= 1-x_bar/dx                            # interpolation function 1
phi_2 = x_bar/dx                             # interpolation function 2

A = np.array(([phi_1], [phi_2]))             # Creaing matrix A(1,2) for the basis of constructing Mass and stiffness matrices 
A=np.transpose(A)
B = np.array(([phi_1], [phi_2]))             # Creaing matrix B(2,1) for the basis of constructing Mass and stiffness matrices
A_by_B=A*B                                   # putting both interpolation functions into a matrix (2 by 2)

#--------------------------------Mass Matrix 'M' in Equation 29 ---------------

# sub_M=np.zeros((N_i,N_i))                    # starting for timing the M_diag_inv calculation
# for i in range(N_i):
#     for j in range(N_i):
#         sub_M[i,j]= sy.integrate(A_by_B[i,j], (x_bar,0,dx))
sub_M3 = np.array([[dx/3,dx/6],[dx/6,dx/3]]) # local mass matrix
M3=np.zeros((nx,nx))                          # generating global mass matrix
i=0
while i<nx:
    M3[i:i+2, i:i+2]= sub_M3[0:2,0:2]
    i+=2

#--------------------------Stifness Matrix 'K' in Equation 29 -----------------
sub_K3=np.array([[-c*dt/2,-c*dt/2],[c*dt/2,c*dt/2]]) # local stifness matrix
K3 = np.zeros((nx,nx))                               # generating global stifness matix
i=0
while i<nx:
    K3[i:i+2, i:i+2]= sub_K3[0:2,0:2]
    i+=2

#-------------------------------Flux in Equation 29----------------------------


#--------------------------------RHS Constant in Equation 29-------------------
RHS_cst3 = (M3 + K3 + F)

##-Matrix method----------------------------------------------------------------
##Mrching forward in time
Un3=np.zeros(nx)                     # dummy vbl to save current values of U (U^t) 
for n in range(nt2):                 # Marching in time
    Un3 = U3.copy()
    RHS3 = RHS_cst3.dot(Un3)           # saving U^t to be used in the next time step calculation
    U3=np.linalg.solve(M3,RHS3)    
    
    if n==1:
        U_plot3[0,:] = U3.copy()      # saving U(t=1)
    if n==int(nt2/2):
        U_plot3[1,:] = U3.copy()      # saving U(t=nt/2)
    if n==int(nt2*0.99):
        U_plot3[2,:] = U3.copy()      # saving U(t= almost the end to time steps)
#------------------------------plot initiation --------------------------------
x = np.linspace(0, L, nx//2)
x = [ele for ele in x for i in range(2)]

plt.figure(1)
plt.axis([0,L, -1,2])
#plt.plot(x, U_plot3[0,:], label='DG-FEM Timestep 1')
#plt.plot(x, U_plot3[1,:], label='Timestep 0.5 nt')
plt.plot(x, U_plot3[2,:], label='DG-FEM Final Timestep')

#.plot(x, U_plot[0,:], label='DPG-FEM Timestep 1')
#plt.plot(x, U_plot[1,:], label='Timestep 0.5 nt')
plt.plot(x, U_plot[2,:], label='DPG-FEM Final Timestep')
plt.xlabel('Distrance')
plt.ylabel('U')
plt.legend()



# # exact solution with tolerant =0
# f = lambda x: 1/4 *(1-x)**2
# sp.integrate.quadrature(f, -1, 1.0)

# ###### generates gaussian points by itself
# def gausss(f,n,a,b):
#     [x,w] = p_roots(n+1)
#     G=0.5*(b-a)*sum(w*f(x))
#     return G

# def func(x):
#     return x**2+3*x**6
# print(gausss(func, 4,-1,1))

# # general Gaussian Quadrature method
# E = np.array([-0.774597, 0.000000, 0.774597])
# A = np.array([0.555556, 0.888889, 0.555556])

# def gauss(f, a, b, E, A):
#     x = np.zeros(len(E))
#     for i in range(len(E)):
#         x[i] = (b+a)/2 + (b-a)/2 *E[i]

#     return (b-a)/2 * (A[0]*f(x[0]) + A[1]*f(x[1]) + A[2]*f(x[2]))

# f = lambda x: 1/4 * (1-x)**2
# a = -1; b = 1

# areaGau = gauss(f, a, b, E, A)
# print("Gaussian integral: ", areaGau)

########## Gaussian quadrature int by specifying order ######
# from scipy.integrate import fixed_quad as G_Q
# def func(x,y):
#     return x**2+3*x**6+y
# G_Q(func, -1, 1, n=4)

########## hard coding points weight for Gaussian Quadrature int and calculating int ######
# points, weights = np.polynomial.legendre.leggauss(3)
# def gausss(f,n,a,b):
#     [x,w] = p_roots(n+1)
#     G=0.5*(b-a)*sum(w*f(x))
#     return G

import sympy as sy

def Jacobian(v_str, f_list):
    vars = sy.symbols(v_str)
    f = sy.sympify(f_list)
    J = sy.zeros(len(f),len(vars))
    for i, fi in enumerate(f):
        for j, s in enumerate(vars):
            J[j,i] = sy.diff(fi, s)
    return J

x, y, a, b, x_m, y_m =sy.symbols('x y a b x_m y_m')

# interpolation functions for the rectangular four noded elements
phi_1 = 1/4*(1-x)*(1-y)
phi_2 = 1/4*(1+x)*(1-y)
phi_3 = 1/4*(1+x)*(1+y)
phi_4 = 1/4*(1-x)*(1+y)

dx = 0.125
dy = 0.16667

X = x_m + dx/2*x                 # X=global x-coordinate, x_m=mid-point, x=natural x-coordinate
Y = y_m + dy/2*y                 # Y=global y-coordinate, y_m=mid-point, y=natural y-coordinate

jacobian = Jacobian('x y', [x_m + dx/2*x,y_m + dy/2*y ])
det = float(sy.det(jacobian))



























