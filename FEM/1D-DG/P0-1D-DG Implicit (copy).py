#--------------------------------1D-DG FEM Implicit----------------------------
import numpy as np
import matplotlib.pyplot as plt
import sympy as sy
import time

nx = 100                   # total number of nodes(degree of freedom)
nt = 20                   # total number of time steps
N_i = 2                     # number of interpolation functions
L =  0.5                    # Totla length
C = .05                     # Courant number
c = .1                      # Wave velocity
#dx = L/(nx//2)               # Distace stepping size
dx = L/(nx//2)               # Distace stepping size
dt = C*dx/c                 # Time stepping size
#x  = np.arange(0, nx)*dx    # or x=np.linspace(0,L,nx)
#print(x)
x = np.linspace(0, L, (nx+2)//2+1)
x = [ele for ele in x for i in range(2)]
#print(x)
x=x[1:-1]
U = np.zeros(nx+2)            # U is a square wave between 0 <U< 1
U_plot = np.ones((3,nx+2))    # A matrix to save 3 time steps used for plotting the results

#------------------------------------------------------------------------------
# Boundary Conditions
U[0] = U[nx-1] = 0          # Dirichlet BC

#------------------------------------------------------------------------------
# Initial conditions
U[int(L*nx*0.37):int(L*nx*1.27)]=1


xi = [-1/np.sqrt(3), 1/np.sqrt(3)] # quad points
w = [1,1] # quad points weights

# shape funcs
phi = {0: 1/4,#lambda xi: 1/2*(1-xi),
       1: 1/4}#lambda xi: 1/2*(1+xi)}

# derivatives of shape funcs
d_dx = {0: 0,
        1: 0}

M = np.zeros((nx,nx))
K = np.zeros((nx,nx))
surf = np.zeros((nx,nx+1))

# global node numbers of each elements
def glob(e):
    node_1=(e+1)*2-2
    node_2= (e+1)*2-1
    return [node_1, node_2]


# coordinates of each point in an element
def coord(e):
        x1 = x[glob(e)[0]]
        x2 = x[glob(e)[1]]
        return [x1, x2]


n_hat ={0:-1,
        1:1}

for e in range(nx//2):
    # det_jacobian
    jac=0  
    for i in range(2):
        jac = jac+ d_dx[i]*coord(e)[i]
    for i_n in range(2):
        global_i = glob(e)[i_n]
        for j_n in range(2):
            global_j = glob(e)[j_n]
            M[global_i , global_j] = 0
            K[global_i , global_j] = 0
            for g in range(len(xi)):
                M[global_i, global_j] = M[global_i , global_j] + 0.25 * dx/2

M_inv = np.linalg.inv(M)

L_xi = [-1,1]
for e in range(nx//2):
    #for s in range(2):
    for s_i in range(2):
        global_si = glob(e)[s_i]
        for s_j in range(2):
            global_sj = glob(e)[s_j]
            
            surf[global_si, glob(e-1)[1]+1] = (n_hat[0] *c*dt* phi[s_i] *
                                                               phi[s_j]) # 99,98
            
            surf[global_si, global_sj+1] = (n_hat[s_i]*c*dt* phi[s_i] *
                                                             phi[s_j])#99,100
surf = -surf[:,1:]

RHS_cst = (M  + surf)

##-Matrix method----------------------------------------------------------------
#Mrching forward in time
#Un=np.zeros(nx)                     # dummy vbl to save current values of U (U^t)
t3 = time.time()
for n in range(nt):                 # Marching in time
    Un = U.copy()
    RHS = RHS_cst.dot(Un[1:nx+1])           # saving U^t to be used in the next time step calculation
    U[1:nx+1]=np.linalg.solve(M,RHS)

    if n==1:
        U_plot[0,:] = U.copy()      # saving U(t=1)
    if n==int(nt/2):
        U_plot[1,:] = U.copy()      # saving U(t=nt/2)
    if n==int(nt*0.99):
        U_plot[2,:] = U.copy()      # saving U(t= almost the end to time steps)
t4 = time.time()
##------------------------------plot initiation --------------------------------
plt.figure(1)
plt.axis([0,L, -1,2])
#plt.plot(x, U_plot[0,:], label='Timestep 1')
#plt.plot(x, U_plot[1,:], label='Timestep 0.5 nt')
plt.plot(x, U_plot[2,:], label='Final Timestep')
plt.xlabel('x-Distrance')
plt.ylabel('U')
plt.legend()
#plt.title(f'Simulation Duration: {round((t4-t3)/60, 2)} minutes')                
        
        
        
        
##------------------------------- Interpolation functions ----------------------
#x_bar=sy.Symbol('x_bar')                     # defining local x symbol for creating interpolation functions and their integrations
#phi_1= 1-x_bar/dx                            # interpolation function 1
#phi_2 = x_bar/dx                             # interpolation function 2
#
#A = np.array(([phi_1], [phi_2]))             # Creaing matrix A(1,2) for the basis of constructing Mass and stiffness matrices
#A=np.transpose(A)
#B = np.array(([phi_1], [phi_2]))             # Creaing matrix B(2,1) for the basis of constructing Mass and stiffness matrices
#A_by_B=A*B                                   # putting both interpolation functions into a matrix (2 by 2)
#
##--------------------------------Mass Matrix 'M' in Equation 29 ---------------
#t1 = time.time()
## sub_M=np.zeros((N_i,N_i))                    # starting for timing the M_diag_inv calculation
## for i in range(N_i):
##     for j in range(N_i):
##         sub_M[i,j]= sy.integrate(A_by_B[i,j], (x_bar,0,dx))
#sub_M = np.array([[dx/3,dx/6],[dx/6,dx/3]]) # local mass matrix
#M=np.zeros((nx,nx))                          # generating global mass matrix
#i=0
#while i<nx:
#    M[i:i+2, i:i+2]= sub_M[0:2,0:2]
#    i+=2
#t2 = time.time()                            # end point of M_diag_inv generation
#print(str(t2-t1))
#
##--------------------------Stifness Matrix 'K' in Equation 29 -----------------
#sub_K=np.array([[-c*dt/2,-c*dt/2],[c*dt/2,c*dt/2]]) # local stifness matrix
#K = np.zeros((nx,nx))                               # generating global stifness matix
#i=0
#while i<nx:
#    K[i:i+2, i:i+2]= sub_K[0:2,0:2]
#    i+=2
#
##-------------------------------Flux in Equation 29----------------------------
#F = np.zeros((nx,nx+1))
#i=0
#j=1
#while i<=nx-1:
#    F[i,i]= c*dt
#    F[j,j+1] = -c*dt
#    i+=2
#    j+=2
#F=F[:,1:]                            # excludig left boundary to get nx by nx matrix
#
##--------------------------------RHS Constant in Equation 29-------------------
#RHS_cst = (M + K + F)

##-Matrix method----------------------------------------------------------------
##Mrching forward in time
#Un=np.zeros(nx)                     # dummy vbl to save current values of U (U^t)
#t3 = time.time()
#for n in range(nt):                 # Marching in time
#    Un = U.copy()
#    RHS = RHS_cst.dot(Un)           # saving U^t to be used in the next time step calculation
#    U=np.linalg.solve(M,RHS)
#
#    if n==1:
#        U_plot[0,:] = U.copy()      # saving U(t=1)
#    if n==int(nt/2):
#        U_plot[1,:] = U.copy()      # saving U(t=nt/2)
#    if n==int(nt*0.99):
#        U_plot[2,:] = U.copy()      # saving U(t= almost the end to time steps)
#t4 = time.time()
#------------------------------plot initiation --------------------------------
#plt.figure(1)
#plt.axis([0,L, -1,2])
#plt.plot(x, U_plot[0,:], label='Timestep 1')
#plt.plot(x, U_plot[1,:], label='Timestep 0.5 nt')
#plt.plot(x, U_plot[2,:], label='Timestep 0.9 nt')
#plt.xlabel('Distrance')
#plt.ylabel('U')
#plt.legend()
#plt.title(f'Simulation Duration: {round((t4-t3)/60, 2)} minutes')        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        