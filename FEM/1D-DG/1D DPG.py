#--------------------------------1D-DPG FEM Implicit----------------------------
import numpy as np
import matplotlib.pyplot as plt
import sympy as sy
import time

nx = 300                   # total number of nodes(degree of freedom)
nt = 500                   # total number of time steps
L =  0.5                    # Totla length
C = .05                     # Courant number
c = .1                      # Wave velocity
dx = L/(nx-1)               # Distace stepping size
dt = C*dx/c                 # Time stepping size
x  = np.arange(0, 2*nx)*dx    # or x=np.linspace(0,2,nx)
U = np.zeros(2*(nx-1))            # U is a square wave between 0 <U< 1
U_plot = np.zeros((3,2*nx))    # A matrix to save 3 time steps used for plotting the results

#------------------------------------------------------------------------------
# Boundary Conditions
U[0] = U[nx-1] = 0          # Dirichlet BC

#------------------------------------------------------------------------------
# Initial conditions
U[int(L*nx*0.1):int(L*nx*0.7)]=1

#------------------- Mass Matrix 'M' for Petrov-Galerkin method ---------------
t1 = time.time()                    # starting for timing the M_diag_inv calculation
# sub_M = np.array([[dx/3-1/2, dx/6-1/2],
#                   [dx/6+1/2, dx/3+1/2]]) # local mass matrix
a=1
M_prime=-dx/16
sub_M = np.array([[dx/3-a*M_prime*1/2  , dx/6-a*M_prime*1/2],
                  [dx/6+a*M_prime*1/2  , dx/3+a*M_prime*1/2]]) # local mass matrix
# sub_M = np.array([[dx/3, dx/6],
#                   [dx/6, dx/3]]) # local mass matrix

M = np.kron(np.eye((nx-1)), sub_M)

t2 = time.time()                            # end point of M_diag_inv generation
print(str(t2-t1))

#-------------- Stifness Matrix 'K' for Petrov-Galerking method ---------------
b=1
K_prime = c*dt*dx/4
sub_K=np.array([[-c*dt/2+b*K_prime/2, -c*dt/2-b*K_prime/2],
                [ c*dt/2-b*K_prime/2,  c*dt/2+b*K_prime/2]]) # local stifness matrix
# sub_K=np.array([[-c/2-1/8, -c/2+1/8],
#                 [c/2+1/8, c/2-1/8]])*dt # local stifness matrix
# sub_K=np.array([[-c/2, -c/2],
#                 [c/2, c/2]])*dt # local stifness matrix
K = np.kron(np.eye((nx-1)), sub_K)

# #-------------------------------Flux in Equation 29----------------------------
F = np.zeros((2*(nx-1),2*(nx-1)+1))
i=0
j=1
while i<=nx-1:
    F[i,i]= c*dt
    F[j,j+1] = -c*dt
    i+=2
    j+=2
F=F[:,1:]                            # excludig left boundary to get nx by nx matrix

#--------------------------------RHS Constant in Equation 29-------------------
RHS_cst = (M + K + F)

##-Matrix method----------------------------------------------------------------
##Mrching forward in time
Un=np.zeros(nx)                     # dummy vbl to save current values of U (U^t) 
t3 = time.time()
for n in range(nt):                 # Marching in time
    Un = U.copy()
    RHS = RHS_cst.dot(Un)           # saving U^t to be used in the next time step calculation
    U=np.linalg.solve(M,RHS)    
    
    if n==1:
        U_plot[0,1:-1] = U.copy()      # saving U(t=1)
    if n==int(nt/2):
        U_plot[1,1:-1] = U.copy()      # saving U(t=nt/2)
    if n==int(nt*0.99):
        U_plot[2,1:-1] = U.copy()      # saving U(t= almost the end to time steps)
t4 = time.time()

# plt.figure(1)
# plt.axis([0,L, -1,2])
# plt.plot(x, U_plot[0,:], label='DPG Timestep 1')
# plt.plot(x, U_plot[1,:], label='DPG Timestep 0.5*nx')
# plt.plot(x, U_plot[2,:], label='DPG Timestep 0.9*nx')

# plt.xlabel('Distrance')
# plt.ylabel('U')
# plt.legend()
# plt.title(f'Simulation Duration: {round((t4-t3)/60, 2)} minutes')

####################### adding standard DG to compare with DPG ################
U2 = np.zeros(2*(nx-1))            # U is a square wave between 0 <U< 1
U_plot2 = np.ones((3,2*nx))    # A matrix to save 3 time steps used for plotting the results

#------------------------------------------------------------------------------
# Boundary Conditions
U2[0] = U2[nx-1] = 0          # Dirichlet BC

#------------------------------------------------------------------------------
# Initial conditions
U2[int(L*nx*0.1):int(L*nx*0.7)]=1

 #--------------------------------Mass Matrix 'M' in Equation 29 ---------------
t1 = time.time()

sub_M2 = np.array([[dx/3, dx/6],
                   [dx/6, dx/3]]) # local mass matrix

M2 = np.kron(np.eye((nx-1)), sub_M2)
t2 = time.time()                            # end point of M_diag_inv generation
print(str(t2-t1))

 #--------------------------Stifness Matrix 'K' in Equation 29 -----------------
sub_K2=np.array([[-c/2, -c/2],
                 [c/2, c/2]])*dt # local stifness matrix
K2 = np.kron(np.eye((nx-1)), sub_K2)

 #--------------------------------RHS Constant in Equation 29-------------------
RHS_cst2 = (M2 + K2 + F)

 ##-Matrix method----------------------------------------------------------------
 ##Mrching forward in time
Un=np.zeros(nx)                     # dummy vbl to save current values of U (U^t) 
t3 = time.time()
for n in range(nt):                 # Marching in time
    Un = U2.copy()
    RHS = RHS_cst2.dot(Un)           # saving U^t to be used in the next time step calculation
    U2=np.linalg.solve(M2,RHS)    
    
    if n==1:
        U_plot2[0,1:-1] = U2.copy()      # saving U(t=1)
    if n==int(nt/2):
        U_plot2[1,1:-1] = U2.copy()      # saving U(t=nt/2)
    if n==int(nt*0.99):
        U_plot2[2,1:-1] = U2.copy()      # saving U(t= almost the end to time steps)
t4 = time.time()
# #------------------------------plot initiation --------------------------------
plt.figure(1)
plt.axis([0,L, -.3,1.3])

#plt.plot(x, U_plot2[0,:], label='DG Timestep 1')
# plt.plot(x, U_plot2[1,:], label='DG Timestep 0.5 nt')
plt.plot(x, U_plot2[2,:], label='DG Timestep 450')

#plt.plot(x, U_plot[0,:], label='DPG Timestep 1')
#plt.plot(x, U_plot[1,:], label='DPG Timestep 0.5 nt')
plt.plot(x, U_plot[2,:], label='DPG Timestep 450')


plt.xlabel('Distrance')
plt.ylabel('U')
plt.legend()
plt.title(f'Simulation Duration: {round((t4-t3)/60, 2)} minutes')

