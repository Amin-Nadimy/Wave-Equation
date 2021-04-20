#-----------------------------------------1D-DG FEM----------------------------------
# This is the Python code for Discontinuous-Galerkin FEM. This code uses the Explicit scheme with the lumped mass matrix for simulation.
# It refers to mathematical formulations in the note in Section 8.3.

import numpy as np
import matplotlib.pyplot as plt
import time

nx = 100                    # total number of nodes(degree of freedom)
nt = 20                   # total number of time steps
L =  0.5                    # Totla length
C = .05                     # Courant number
c = .1                      # Wave velocity
dx = L/(nx//2)               # Distace stepping size
dt = C*dx/c                 # Time stepping size
x  = np.arange(0, nx)*dx    # or x=np.linspace(0,2,nx)
U = np.zeros(nx)            # U is a square wave between 0 <U< 1
U_plot = np.ones((3,nx))    # A matrix to save 3 time steps used for plotting the results

#------------------------------------------------------------------------------------
# Boundary Conditions
U[0] = U[nx-1] = 0                   # Dirichlet BC

#------------------------------------------------------------------------------------
# Initial conditions
U[int(L*nx*0.2):int(L*nx*0.5)]=1     # defining square wave shape

#-----------------------Explicit using lumped mass matrix---------------------------
# Un=np.zeros(nx)                      # dummy vbl to save current values of U (U^t) 
# t3 = time.time() 
# for n in range(nt):                  # Marching in time
#     Un = U.copy()                    # saving U^t to be used in the next time step calculation
#     U[1] = Un[1] - c* dt/dx *(Un[1]+Un[2])
#     i=2                              
#     while i<nx-1:                    
#         U[i] = Un[i] + (-1)**i*c*dt/dx*(Un[i] + Un[i-1]) + (-1)**(i+1)*c*Un[i]
#         i +=1
#         if i==nx-1:             
#             StopIteration
#         else:
#             U[i] = Un[i] + (-1)**i*c*dt/dx*(Un[i] + Un[i+1]) + (-1)**(i+1)*c*Un[i-1]
#         i +=1
    
#     if n==1:
#         U_plot[0,:] = U.copy()      # saving U(t=1)
#     if n==int(nt/2):
#         U_plot[1,:] = U.copy()      # saving U(t=nt/2)
#     if n==int(nt*0.99):
#         U_plot[2,:] = U.copy()      # saving U(t= almost the end to time steps)
# t4 = time.time() 
#---------------------------Mass Matrix 'M' in Equation 55 --------------------------
t1 = time.time()                    # starting for timing the M_diag_inv calculation
sub_M = np.array([[dx/3,dx/6],[dx/6,dx/3]]) # local mass matrix
M=np.zeros((nx-2,nx-2))                         # generating global mass matrix
i=0
while i<nx-2:
    M[i:i+2, i:i+2]= sub_M[0:2,0:2]
    i+=2
M_inv=np.linalg.inv(M)
t2 = time.time()                            # end point of M_diag_inv generation
print(str(t2-t1))
#--------------------------Stifness Matrix 'K' in Equation 55------------------------
sub_K=np.array([[-c*dt/2,-c*dt/2],[c*dt/2,c*dt/2]]) # local stifness matrix
K = np.zeros((nx-2,nx-2))                               # generating global stifness matix
i=0
while i<nx-2:
    K[i:i+2, i:i+2]= sub_K[0:2,0:2]
    i+=2

# #-------------------------------Flux matrix in Equation 55--------------------------
F = np.zeros((nx-2,nx-1))
i=0
j=1
while i<=nx-3:
    F[i,i]= c*dt
    F[j,j+1] = -c*dt
    i+=2
    j+=2
F=F[:,1:]                            # excludig left boundary to get nx by nx matrix

#--------------------------------RHS in equation 29----------------------------------
RHS_cst = M_inv.dot((M + K + F))

#----------------------------Explicit using consistance mass matrix------------------
Un=np.zeros(nx)                      # dummy vbl to save current values of U (U^t) 
t3 = time.time() 
for n in range(nt):                  # Marching in time
    Un = U.copy()                    # saving U^t to be used in the next time step calculation
    U[1] = RHS_cst[0,0]*Un[1] + RHS_cst[0,1]*Un[2]
    U[2] = RHS_cst[1,0]*Un[1] + RHS_cst[1,1]*Un[2]
    i=3
    j=1
    while i<nx-2:
        U[i] = RHS_cst[i-1,j]*Un[i-1] + RHS_cst[i-1,j+1]*Un[i] + RHS_cst[i-1,j+2]*Un[i+1]
        i+=1
        U[i] = RHS_cst[i-1,j]*Un[i-2] + RHS_cst[i-1,j+1]*Un[i-1] + RHS_cst[i-1,j+2]*Un[i]
        i+=1
        j+=2
    if n==1:
        U_plot[0,:] = U.copy()       # saving U(t=1)
    if n==int(nt/2):
        U_plot[1,:] = U.copy()       # saving U(t=nt/2)
    if n==int(nt*0.99):
        U_plot[2,:] = U.copy()       # saving U(t= almost the end to time steps)
t4 = time.time()              
#------------------------------plot initiation --------------------------------------
plt.figure(1)
plt.axis([0,L, -1,2])
plt.plot(x, U_plot[0,:], label='Timestep 1')
plt.plot(x, U_plot[1,:], label='Timestep 0.5 nt')
plt.plot(x, U_plot[2,:], label='Timestep 0.9 nt')
plt.xlabel('Distrance')
plt.ylabel('U')
plt.legend()
plt.title(f'Simulation Duration: {round((t4-t3)/60, 2)} minutes')


# an_array = np.array([[1,0,0],[0,0,-1]])
# repetitions = 3
# B=np.eye(4,4)
# repeats_array = np.tile(an_array, (repetitions, 1))

# print(repeats_array)

# # Your initial 100 x 100 matrix 
# a = np.zeros((100, 100))


# for i in range(10):
#   b = np.ones((10, 10)) * (i + 1)
#   a[i*10:(i+1)*10,i*10:(i+1)*10] = b


# a = np.ones((100,100)); b = np.ones((10,10))*2;

# np.diagonal(a)[:] = np.ravel(b)

# i=0
# j=0
# a=np.zeros((6,6))
# b=np.array([[1,2],[3,4]])

# print(b.shape)
# while i<5:
#     a[i:i+2, j:j+2]= b
#     i+=2
#     j+=2


# nx=4
# b=np.array([[1,0,0],[0,0,-1]])
# a = np.zeros((nx,nx+1))
# i=0
# j=0
# while i<=int(nx+1)/2:
#     a[i:i+2, j:j+3]= b
#     i+=2
#     j+=2
# F=a[1:,1:]

# v = np.arange(36).reshape(6,6)
# ptotal = np.arange(6)
# sum(v[ptotal,:])
# v.sum(axis=1)

# M_eigenvalues, M_eigenvectors = np.linalg.eig(M)    # Calculatingeigenvalues and vectors of M
# M_eig_inv = np.linalg.inv(M_eigenvectors)           # calculating inverse of M eigenvector
# M_diag = M_eig_inv.dot(M).dot(M_eigenvectors)       # Diagonalising M
# M_diag_inv = np.linalg.inv(M_diag)                  # inverse of M_diag used in Equ 30


