# #--------------------------------1D-DG FEM Implicit----------------------------
# import numpy as np
# import matplotlib.pyplot as plt
# import sympy as sy
# import time

# nx = 500                   # total number of nodes(degree of freedom)
# nt = 1000                   # total number of time steps
# N_i = 2                     # number of interpolation functions
# L =  0.5                    # Totla length
# C = .05                     # Courant number
# c = .1                      # Wave velocity
# dx = L/(nx-1)               # Distace stepping size
# dt = C*dx/c                 # Time stepping size
# x  = np.arange(0, nx)*dx    # or x=np.linspace(0,2,nx)
# U = np.zeros(nx)            # U is a square wave between 0 <U< 1
# U_plot = np.ones((3,nx))    # A matrix to save 3 time steps used for plotting the results

# #------------------------------------------------------------------------------
# # Boundary Conditions
# U[0] = U[nx-1] = 0          # Dirichlet BC

# #------------------------------------------------------------------------------
# # Initial conditions
# U[int(L*nx*0.2):int(L*nx*0.5)]=1

# #------------------------------- Interpolation functions ----------------------
# x_bar=sy.Symbol('x_bar')                     # defining local x symbol for creating interpolation functions and their integrations
# phi_1= 1-x_bar/dx                            # interpolation function 1
# phi_2 = x_bar/dx                             # interpolation function 2

# A = np.array(([phi_1], [phi_2]))             # Creaing matrix A(1,2) for the basis of constructing Mass and stiffness matrices 
# A=np.transpose(A)
# B = np.array(([phi_1], [phi_2]))             # Creaing matrix B(2,1) for the basis of constructing Mass and stiffness matrices
# A_by_B=A*B                                   # putting both interpolation functions into a matrix (2 by 2)

# #--------------------------------Mass Matrix 'M' in Equation 29 ---------------
# t1 = time.time()
# sub_M=np.zeros((N_i,N_i))                    # starting for timing the M_diag_inv calculation
# for i in range(N_i):
#     for j in range(N_i):
#         sub_M[i,j]= sy.integrate(A_by_B[i,j], (x_bar,0,dx))
# #sub_M = np.array([[dx/3,dx/6],[dx/6,dx/3]]) # local mass matrix
# M=np.zeros((nx,nx))                          # generating global mass matrix
# i=0
# while i<nx:
#     M[i:i+2, i:i+2]= sub_M[0:2,0:2]
#     i+=2
# t2 = time.time()                            # end point of M_diag_inv generation
# print(str(t2-t1))

# #--------------------------Stifness Matrix 'K' in Equation 29 -----------------
# sub_K=np.array([[-c*dt/2,-c*dt/2],[c*dt/2,c*dt/2]]) # local stifness matrix
# K = np.zeros((nx,nx))                               # generating global stifness matix
# i=0
# while i<nx:
#     K[i:i+2, i:i+2]= sub_K[0:2,0:2]
#     i+=2

# #-------------------------------Flux in Equation 29----------------------------
# F = np.zeros((nx,nx+1))
# i=0
# j=1
# while i<=nx-1:
#     F[i,i]= c*dt
#     F[j,j+1] = -c*dt
#     i+=2
#     j+=2
# F=F[:,1:]                            # excludig left boundary to get nx by nx matrix

# #--------------------------------RHS Constant in Equation 29-------------------
# RHS_cst = (M + K + F)

# ##-Matrix method----------------------------------------------------------------
# ##Mrching forward in time
# Un=np.zeros(nx)                     # dummy vbl to save current values of U (U^t) 
# t3 = time.time()
# for n in range(nt):                 # Marching in time
#     Un = U.copy()
#     RHS = RHS_cst.dot(Un)           # saving U^t to be used in the next time step calculation
#     U=np.linalg.solve(M,RHS)    
    
#     if n==1:
#         U_plot[0,:] = U.copy()      # saving U(t=1)
#     if n==int(nt/2):
#         U_plot[1,:] = U.copy()      # saving U(t=nt/2)
#     if n==int(nt*0.99):
#         U_plot[2,:] = U.copy()      # saving U(t= almost the end to time steps)
# t4 = time.time()
# #------------------------------plot initiation --------------------------------
# plt.figure(1)
# plt.axis([0,L, -1,2])
# plt.plot(x, U_plot[0,:], label='Timestep 1')
# plt.plot(x, U_plot[1,:], label='Timestep 0.5 nt')
# plt.plot(x, U_plot[2,:], label='Timestep 0.9 nt')
# plt.xlabel('Distrance')
# plt.ylabel('U')
# plt.legend()
# plt.title(f'Simulation Duration: {round((t4-t3)/60, 2)} minutes')



import numpy as np
import matplotlib.pyplot as plt
import sympy as sy
import time

N_e_r = 20  
N_e_c = 30
N_p = 4
nx =  N_e_c * N_e_r * 2                
ny =  N_e_c * N_e_r * 2
nt = 1000                   # total number of time steps
N_i = 2                     # numbe rof interpolation functions
L =  0.5                    # Totla length
C = .05                     # Courant number
c = .1 

# flux = np.zeros((2*nx,2*ny))

# i = N_e_c+1
# j = i+1

# for n in range(N_e_c):
#     flux[i,i]= 1
#     flux[i+3*N_e_c,i+3*N_e_c] = -1
#     i+=1
#     j+=1
# print(i)
# plt.spy(flux)
# flux[i+1,i] = 0

sub_flux = np.zeros((N_e_r*N_p,N_e_r*N_p+N_e_r))

i = 0
# applying y-dir
for n in range(N_e_r):
    sub_flux[i,i]= -1
    sub_flux[i+3*N_e_r,i+4*N_e_r] = 1
    i+=1
#applying x-dir
# j = N_e_r+1
# for n in range(N_e_c):
#     sub_flux[j,j+N_e_r] = 1
#     sub_flux[j+1,j+N_e_r]=-1
#     j+=2
# sub_flux[3*N_e_r,4*N_e_r-1]=0
fx = np.zeros((N_e_r*2,N_e_r*2+1))
i=0
j=1
for n in range(N_e_r):
    fx[i,i]= 1
    fx[j,j+1] = -1
    i+=2
    j+=2
fx=fx[:,1:]

sub_flux[N_e_r:N_e_r+fx.shape[0], N_e_r*2:N_e_r*2+fx.shape[1]]+=fx



flux = np.zeros((N_p*N_e_r*N_e_c,N_p*N_e_r*N_e_c+N_e_r))
i=0
j=0
for n in range(N_e_c):
    flux[i:i+sub_flux.shape[0], j:j+sub_flux.shape[1]]+=sub_flux
    i+=N_e_r*N_p
    j+=N_e_r*N_p
flux = flux[:,N_e_r:]
# plt.spy(sub_flux)
plt.spy(flux)    







