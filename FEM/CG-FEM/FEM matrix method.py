import numpy as np
import matplotlib.pyplot as plt

nx = 100                    # total number of nodes(degree of freedom)
nt = 1000                   # total number of time steps
L =  0.5                     # Totla length
C = .001                     # Courant number
c = .1                      # Wave velocity
dx = L/(nx-1)               # Distace stepping size
dt = C*dx/c                 # Time stepping size
x  = np.arange(0, nx)*dx    # or x=np.linspace(0,2,nx)
U = np.zeros(nx)            # U is a square wave between 0 <U< 1
U_plot = np.ones((3,nx))    # A matrix to save 3 time steps used for plotting the results

#------------------------------------------------------------------------------
#Boundary Conditions
#U[0] = U[nx-1] = 0          # Dirichlet BC

#------------------------------------------------------------------------------
#Initial conditions
U[int(L*nx*0.2):int(L*nx*0.5)]=1

#------------------------------------------------------------------------------
# Matrix A it is the LHS (mass matrix) shown in Equation 14
# Matrix A it is the LHS (mass matrix) shown in Equation 14
A = np.zeros((nx-2,nx-2))       
for i in range(nx-3):
    for j in range (nx-3):
        if j==i:
            A[i,j] = dx/3*2
        elif j==i+1:
            A[i,j] = dx/6
        elif j==i-1:
            A[i,j] = dx/6
        else:
            A[i,j] = 0
            
#corner cells
A[0,0] = A[nx-3, nx-3] = dx/3
A[1,0] = A[0,1] = A[nx-3, nx-4] = A[nx-4, nx-3] = dx/6

#A_inv= np.linalg.inv(A)         # Inverse of matrix A used in Equation 15

#------------------------------------------------------------------------------
B = np.zeros((nx-2, nx-2))          # Matrix B in Equation 15
for i in range(nx-2):
    for j in range(nx-2):
        if  i==j+1:
            B[i,j] = 0.5
        elif i==j-1:
            B[i,j] = -0.5
B[0,0]=-0.5
B[nx-3,nx-3]=0.5

#----------------------------A.U^{t+1}=(A+c.dt.B).U^{t}---------equ.27-------
# E = Matrix in formula 27. E= A+c.dt.B
E = A+c*dt*B

#-Matrix method----------------------------------------------------------------
#Marching forward in time
#rhs=np.zeros(nx)
Un=np.zeros(nx-2)
for n in range(1, nt):                  # Marching in time
    Un = U[1:nx-1].copy()                       # un is a dummy vbl for the current known quantity
    #for j in range(1, nx):              # Marching in elements (fundamental difference)
    rhs = np.matmul(E,Un)               # making Eu in the RHS in equ 27
    U[1:nx-1] = np.linalg.solve(A, rhs)     # linear solver to solve equ 27
    if n==1:
        U_plot[0,:] = U.copy()
    if n==int(nt/2):
        U_plot[1,:] = U.copy()
    if n==int(nt*0.99):
        U_plot[2,:] = U.copy()

plt.figure(1)
plt.axis([0,L, -1,1.5])

# Plotting FEM-----------------------------------------------------------------
plt.plot(x, U_plot[0,:], label='Timestep 1')
plt.plot(x, U_plot[1,:], label='Timestep 0.5*nx')
plt.plot(x, U_plot[2,:], label='Timestep 0.9*nx')
plt.xlabel('Distrance')
plt.ylabel('U')
plt.legend()

          


            


#------------------------------------------------------------------------------
