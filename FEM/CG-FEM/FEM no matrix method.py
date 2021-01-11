import numpy as np
import matplotlib.pyplot as plt

nx = 500                        # distance which the wave travels
nt = 1000                       # total number of time steps
L =  0.5                        # Totla length
C = 0.01                        # Courant number
c = 0.1                         # Wave velocity
dx = L/(nx-1)                   # Distace stepping size
dt = C*dx/c                     # Time stepping size
x  = np.arange(0, nx)*dx        # or x=np.linspace(0,2,nx)
U = np.zeros(nx)                # U is a square wave between 0 <U< 1
U_FD = np.ones(nx)*1.5          # U is a square wave between 0 <U< 1
U_plot = np.ones((3,nx))        # A matrix to save 3 time steps used for plotting the results of FEM
U_FD_plot = np.ones((3,nx))     # An matrix to save 3 time steps used for plotting the results of FDM

#------------------------------------------------------------------------------
#Boundary Conditions
U[0] = U[nx-1] = 0              # Dirichlet BC for 
U_FD[0] = U_FD[nx-1] = 1.5      # Dirichlet BC for FDM

#------------------------------------------------------------------------------
#Initial conditions
U[int(.1*L*nx):int(.5*L*nx)]=1
U_FD[int(.1*L*nx):int(.5*L*nx)]=2.5

#------------------------------------------------------------------------------
# Matrix A it is the LHS (mass matrix) shown in Equation 14
A = np.zeros((nx,nx))       
for i in range(nx-1):
    for j in range (nx-1):
        # defining mass matrix manually
        if j==i:
            A[i,j] = dx/3*2
        elif j==i+1:
            A[i,j] = dx/6
        elif j==i-1:
            A[i,j] = dx/6
        else:
            A[i,j] = 0
            
# defining corner cells
A[0,0] = A[nx-1, nx-1] = dx/3
A[1,0] = A[0,1] = A[nx-1, nx-2] = A[nx-2, nx-1] = dx/6

A_inv= np.linalg.inv(A)          # Inverse of matrix A used in Equation 15

#------------------------------------------------------------------------------
B = np.zeros((nx, nx))           # Matrix B in Equation 15
for i in range(nx):
    for j in range(nx):
        # defining stiffness matrix manually
        if  i==j+1:
            B[i,j] = 0.5
        elif i==j-1:
            B[i,j] = -0.5

# defining corner cells
B[0,0]=-0.5
B[nx-1,nx-1]=0.5

#------------------------- Explicit scheme ----------------------------------------
# u(t+1)= u(t)+cdt A^(-1) B u(t)
dummy=A_inv.dot(B)               # A dummy vbl to make the RHS of Equation 15 clearer
                                 # Dummy is a nx by nx matrix
for n in range(1,nt):
    Un=U.copy()
    Un_FD=U_FD.copy()
    dummy2 = dummy.dot(Un)       # Dummy 2 creats a nx by 1 matrix
    for j in range(1,nx):
        U[j]=Un[j] + c*dt*dummy2[j]     # Overall equation
        U_FD[j] = U_FD[j]-c*(dt/dx)*(Un_FD[j]-Un_FD[j-1])
    
        # saving the solutions in particular time steps
    if n==1:
        U_plot[0,:] = U.copy()
        U_FD_plot[0,:] = U_FD.copy()
    if n==int(nt/2):
        U_plot[1,:] = U.copy()
        U_FD_plot[1,:] = U_FD.copy()
    if n==int(nt*0.95):
        U_plot[2,:] = U.copy()
        U_FD_plot[2,:] = U_FD.copy()

# ------------------------- Plotting the solutions ------------------------------    
plt.figure(1)
plt.axis([0,L, -1,3])
plt.xlabel('Distrance')
plt.ylabel('U')

# ----------------------- Plotting FDM and FEM for comparision ------------------
# Plotting FEM-----------------------------------------------------------------
plt.plot(x, U_plot[0,:], label='Timestep 1 FE')
plt.plot(x, U_plot[1,:], label='Timestep 0.5 nt FE')
plt.plot(x, U_plot[2,:], label='Timestep 0.9 nt FE')

# Plotting FDM-----------------------------------------------------------------
plt.plot(x, U_FD_plot[0,:], label='Timestep 1 FD')
plt.plot(x, U_FD_plot[1,:], label='Timestep 0.5 nt FD')
plt.plot(x, U_FD_plot[2,:], label='Timestep 0.9 nt FD')

plt.legend()





























