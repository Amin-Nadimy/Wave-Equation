import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sympy as sy

# double integral with Gaussian Quadrature
# generates gauusian points and weights based on the degree of polynomila 2n-1. e.g. 3 generates 3 points guadrature.
points, weights = np.polynomial.legendre.leggauss(3)        

def gauss(f, a, b, points, weights):
    y=x = np.zeros(len(points))
    for i in range(len(points)):
        x[i] = (b+a)/2 + (b-a)/2 *points[i]
        y[i] = (b+a)/2 + (b-a)/2 *points[i]
    answer =0
    for i in range(len(points)):
        for j in range(len(points)):
            answer =  answer+(b-a)/2 * (weights[i]*weights[j]*f(x[i], y[j]))
    return answer

C= .05                                      # CLF number
Np = 4                                      # number rof points at each element
nt = 20 
N_e_r = 4  
N_e_c = 3
nx =  N_e_c * N_e_r * 2                     # total number of x steps               
ny =  N_e_c * N_e_r * 2                     # total number of y steps                                                                 
N = N_e_c * N_e_r                           # number of elements
L = 0.5                                       # x and y lengths
dx = L/(N_e_r)
dy = L/(N_e_c)
c_x = 0.1                                   # Wave velocity in x-dir
c_y = 0.1                                 # Wave velocity in y-dir
dt = C*dx*dy/(c_x*dy+c_y*dx)

U = np.zeros(N*Np)                          # Wave matrix
Un = np.zeros(N*Np)                         # Dummy variable to save current components of U
U1=U2 = U                                # Dummy matrices to plot 3 time steps
#U_plot = np.ones((3,N*Np))                 # A matrix to save 3 time steps used for plotting the results

U[int(N*Np*.3):int(N*Np*.8)]=1              # Defining wave components

#--------------------------------- Initial Conds ------------------------------
U[0]= U[-1] = 0
#--------------------- Interpolation functions phi(i) -------------------------
x, y, a, b, x_m, y_m =sy.symbols('x y a b x_m y_m')
# phi_1 = 1/4*(1-x)*(1-y)
# phi_2 = 1/4*(1+x)*(1-y)
# phi_3 = 1/4*(1-x)*(1+y)
# phi_4 = 1/4*(1+x)*(1+y)

def Jacobian(v_str, f_list):
    vars = sy.symbols(v_str)
    f = sy.sympify(f_list)
    J = sy.zeros(len(f),len(vars))
    for i, fi in enumerate(f):
        for j, s in enumerate(vars):
            J[j,i] = sy.diff(fi, s)
    return J

jacobian = Jacobian('x y', [x_m + dx/2*x,y_m + dy/2*y])
det = float(sy.det(jacobian))
################### Mass matrix M1 = int ( phi_i phi_j ) ######################
sub_M_00 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1-y)*det
sub_M_01 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1-y)*det
sub_M_02 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1+y)*det
sub_M_03 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1+y)*det

sub_M=np.eye(4)*gauss(sub_M_00, -1, 1, points, weights)
sub_M = sub_M + np.diag(np.ones(3)*gauss(sub_M_01, -1, 1, points, weights),1)
sub_M = sub_M + np.diag(np.ones(3)*gauss(sub_M_01, -1, 1, points, weights),-1)
sub_M = sub_M + np.diag(np.ones(2)*gauss(sub_M_02, -1, 1, points, weights),2)
sub_M = sub_M + np.diag(np.ones(2)*gauss(sub_M_02, -1, 1, points, weights),-2)
sub_M[3,0] = sub_M[0,3] = gauss(sub_M_03, -1, 1, points, weights)
M = np.kron(np.eye(N), sub_M)                                       # Creating global mass matrix

############################## stiffness matrix K ############################
sub_K_00 = lambda x,y: c_x*dt * -1/4*(1-y) * 1/4*(1-x)*(1-y) * dy/2 + c_y*dt *-1/4*(1-x) * 1/4*(1-x)*(1-y) * dx/2
sub_K_01 = lambda x,y: c_x*dt * -1/4*(1-y) * 1/4*(1+x)*(1-y) * dy/2 + c_y*dt *-1/4*(1-x) * 1/4*(1+x)*(1-y) * dx/2
sub_K_02 = lambda x,y: c_x*dt * -1/4*(1-y) * 1/4*(1-x)*(1+y) * dy/2 + c_y*dt *-1/4*(1-x) * 1/4*(1-x)*(1+y) * dx/2
sub_K_03 = lambda x,y: c_x*dt * -1/4*(1-y) * 1/4*(1+x)*(1+y) * dy/2 + c_y*dt *-1/4*(1-x) * 1/4*(1+x)*(1+y) * dx/2

sub_K_10 = lambda x,y: c_x*dt * 1/4*(1-y) * 1/4*(1-x)*(1-y) * dy/2 + c_y*dt *-1/4*(1+x) * 1/4*(1-x)*(1-y) * dx/2
sub_K_11 = lambda x,y: c_x*dt * 1/4*(1-y) * 1/4*(1+x)*(1-y) * dy/2 + c_y*dt *-1/4*(1+x) * 1/4*(1+x)*(1-y) * dx/2
sub_K_12 = lambda x,y: c_x*dt * 1/4*(1-y) * 1/4*(1-x)*(1+y) * dy/2 + c_y*dt *-1/4*(1+x) * 1/4*(1-x)*(1+y) * dx/2
sub_K_13 = lambda x,y: c_x*dt * 1/4*(1-y) * 1/4*(1+x)*(1+y) * dy/2 + c_y*dt *-1/4*(1+x) * 1/4*(1+x)*(1+y) * dx/2

sub_K_20 = lambda x,y: c_x*dt * -1/4*(1+y) * 1/4*(1-x)*(1-y) * dy/2 + c_y*dt *1/4*(1-x) * 1/4*(1-x)*(1-y) * dx/2
sub_K_21 = lambda x,y: c_x*dt * -1/4*(1+y) * 1/4*(1+x)*(1-y) * dy/2 + c_y*dt *1/4*(1-x) * 1/4*(1+x)*(1-y) * dx/2
sub_K_22 = lambda x,y: c_x*dt * -1/4*(1+y) * 1/4*(1-x)*(1+y) * dy/2 + c_y*dt *1/4*(1-x) * 1/4*(1-x)*(1+y) * dx/2
sub_K_23 = lambda x,y: c_x*dt * -1/4*(1+y) * 1/4*(1+x)*(1+y) * dy/2 + c_y*dt *1/4*(1-x) * 1/4*(1+x)*(1+y) * dx/2

sub_K_30 = lambda x,y: c_x*dt * 1/4*(1+y) * 1/4*(1-x)*(1-y) * dy/2 + c_y*dt *1/4*(1+x) * 1/4*(1-x)*(1-y) * dx/2
sub_K_31 = lambda x,y: c_x*dt * 1/4*(1+y) * 1/4*(1+x)*(1-y) * dy/2 + c_y*dt *1/4*(1+x) * 1/4*(1+x)*(1-y) * dx/2
sub_K_32 = lambda x,y: c_x*dt * 1/4*(1+y) * 1/4*(1-x)*(1+y) * dy/2 + c_y*dt *1/4*(1+x) * 1/4*(1-x)*(1+y) * dx/2
sub_K_33 = lambda x,y: c_x*dt * 1/4*(1+y) * 1/4*(1+x)*(1+y) * dy/2 + c_y*dt *1/4*(1+x) * 1/4*(1+x)*(1+y) * dx/2

sub_K=np.zeros((4,4))
sub_K[0,0]=gauss(sub_K_00, -1, 1, points, weights)
sub_K[0,1]=gauss(sub_K_01, -1, 1, points, weights)
sub_K[0,2]=gauss(sub_K_02, -1, 1, points, weights)
sub_K[0,3]=gauss(sub_K_03, -1, 1, points, weights)

sub_K[1,0]=gauss(sub_K_10, -1, 1, points, weights)
sub_K[1,1]=gauss(sub_K_11, -1, 1, points, weights)
sub_K[1,2]=gauss(sub_K_12, -1, 1, points, weights)
sub_K[1,3]=gauss(sub_K_13, -1, 1, points, weights)

sub_K[2,0]=gauss(sub_K_20, -1, 1, points, weights)
sub_K[2,1]=gauss(sub_K_21, -1, 1, points, weights)
sub_K[2,2]=gauss(sub_K_22, -1, 1, points, weights)
sub_K[2,3]=gauss(sub_K_23, -1, 1, points, weights)

sub_K[3,0]=gauss(sub_K_30, -1, 1, points, weights)
sub_K[3,1]=gauss(sub_K_31, -1, 1, points, weights)
sub_K[3,2]=gauss(sub_K_32, -1, 1, points, weights)
sub_K[3,3]=gauss(sub_K_33, -1, 1, points, weights)

K = np.kron(np.eye(N), sub_K)            # Creating global stiffness matrix

################################ flux wiht #########################################
# x-dir line connecting points 1 and 3 with positive upwind flux n_x = -1
def gauss_1_3(f, a, b, points, weights):
    y = np.zeros(3)
    for i in range(3):
        y[i] = (b+a)/2 + (b-a)/2 *points[i]
        x=-1

    return (b-a)/2 * (weights[0]*f(x,y[0]) + weights[1]*f(x,y[1]) + weights[2]*f(x,y[2]))

sub_F_x_00 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1-y)*dy/2*c_x*dt*(-1)
sub_F_x_01 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1-y)*dy/2*c_x*dt*(-1)
sub_F_x_02 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1+y)*dy/2*c_x*dt*(-1)
sub_F_x_03 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1+y)*dy/2*c_x*dt*(-1)

sub_F_x_20 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1-y)*dy/2*c_x*dt*(-1)
sub_F_x_21 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1-y)*dy/2*c_x*dt*(-1)
sub_F_x_22 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1+y)*dy/2*c_x*dt*(-1)
sub_F_x_23 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1+y)*dy/2*c_x*dt*(-1)

sub_F_x=np.zeros((Np,Np))

sub_F_x[0,0]=gauss_1_3(sub_F_x_00, -1, 1, points, weights)
sub_F_x[0,1]=gauss_1_3(sub_F_x_01, -1, 1, points, weights)
sub_F_x[0,2]=gauss_1_3(sub_F_x_02, -1, 1, points, weights)
sub_F_x[0,3]=gauss_1_3(sub_F_x_03, -1, 1, points, weights)

sub_F_x[2,0]=gauss_1_3(sub_F_x_20, -1, 1, points, weights)
sub_F_x[2,1]=gauss_1_3(sub_F_x_21, -1, 1, points, weights)
sub_F_x[2,2]=gauss_1_3(sub_F_x_22, -1, 1, points, weights)
sub_F_x[2,3]=gauss_1_3(sub_F_x_23, -1, 1, points, weights)

# x-dir line connecting points 2 and 4 with positive upwind flux n_x = +1
def gauss_2_4(f, a, b, points, weights):
    y = np.zeros(3)
    for i in range(3):
        y[i] = (b+a)/2 + (b-a)/2 *points[i]
        x=1

    return (b-a)/2 * (weights[0]*f(x,y[0]) + weights[1]*f(x,y[1]) + weights[2]*f(x,y[2]))

sub_F_x_10 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1-y)*dy/2*-c_x*dt*(+1)
sub_F_x_11 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1-y)*dy/2*-c_x*dt*(+1)
sub_F_x_12 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1+y)*dy/2*c_x*dt*(+1)
sub_F_x_13 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1+y)*dy/2*c_x*dt*(+1)

sub_F_x_30 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1-y)*dy/2*-c_x*dt*(+1)
sub_F_x_31 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1-y)*dy/2*-c_x*dt*(+1)
sub_F_x_32 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1+y)*dy/2*c_x*dt*(+1)
sub_F_x_33 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1+y)*dy/2*c_x*dt*(+1)

sub_F_x[1,0]=gauss_2_4(sub_F_x_10, -1, 1, points, weights)
sub_F_x[1,1]=gauss_2_4(sub_F_x_11, -1, 1, points, weights)
sub_F_x[1,2]=gauss_2_4(sub_F_x_12, -1, 1, points, weights)
sub_F_x[1,3]=gauss_2_4(sub_F_x_13, -1, 1, points, weights)

sub_F_x[3,0]=gauss_2_4(sub_F_x_30, -1, 1, points, weights)
sub_F_x[3,1]=gauss_2_4(sub_F_x_31, -1, 1, points, weights)
sub_F_x[3,2]=gauss_2_4(sub_F_x_32, -1, 1, points, weights)
sub_F_x[3,3]=gauss_2_4(sub_F_x_33, -1, 1, points, weights)

################################## y -dir #####################################
# y-dir line connecting points 1 weightsnd 2 with positive upwind flux n_y = -1
def gauss_1_2(f, a, b, points, weights):
    x = np.zeros(3)
    for i in range(3):
        x[i] = (b+a)/2 + (b-a)/2 *points[i]
        y=-1

    return (b-a)/2 * (weights[0]*f(x[0],y) + weights[1]*f(x[1],y) + weights[2]*f(x[2],y))

sub_F_y_00 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1-y)*det*-c_y*dt*(-1)
sub_F_y_01 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1-y)*det*-c_y*dt*(-1)
sub_F_y_02 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1+y)*det*c_y*dt*(-1)
sub_F_y_03 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1+y)*det*c_y*dt*(-1)

sub_F_y_10 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1-y)*det*-c_y*dt*(-1)
sub_F_y_11 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1-y)*det*-c_y*dt*(-1)
sub_F_y_12 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1+y)*det*c_y*dt*(-1)
sub_F_y_13 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1+y)*det*c_y*dt*(-1)

sub_F_y=np.zeros((Np,Np))
sub_F_y[0,0]=gauss_1_2(sub_F_y_00, -1, 1, points, weights)
sub_F_y[0,1]=gauss_1_2(sub_F_y_01, -1, 1, points, weights)
sub_F_y[0,2]=gauss_1_2(sub_F_y_02, -1, 1, points, weights)
sub_F_y[0,3]=gauss_1_2(sub_F_y_03, -1, 1, points, weights)

sub_F_y[1,0]=gauss_1_2(sub_F_y_10, -1, 1, points, weights)
sub_F_y[1,1]=gauss_1_2(sub_F_y_11, -1, 1, points, weights)
sub_F_y[1,2]=gauss_1_2(sub_F_y_12, -1, 1, points, weights)
sub_F_y[1,3]=gauss_1_2(sub_F_y_13, -1, 1, points, weights)


# y-dir line connecting points 3 and 4 with positive upwind flux n_y = +1
def gauss_3_4(f, a, b, points, weights):
    x = np.zeros(3)
    for i in range(3):
        x[i] = (b+a)/2 + (b-a)/2 *points[i]
        y=+1

    return (b-a)/2 * (weights[0]*f(x[0],y) + weights[1]*f(x[1],y) + weights[2]*f(x[2],y))

sub_F_y_20 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1-y)*dx/2*(-c_y)*dt
sub_F_y_21 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1-y)*dx/2*-c_y*dt
sub_F_y_22 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1+y)*dx/2*c_y*dt
sub_F_y_23 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1+y)*dx/2*c_y*dt

sub_F_y_30 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1-y)*dx/2*-c_y*dt
sub_F_y_31 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1-y)*dx/2*-c_y*dt
sub_F_y_32 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1+y)*dx/2*c_y*dt
sub_F_y_33 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1+y)*dx/2*c_y*dt

sub_F_y[2,0]=gauss_3_4(sub_F_y_20, -1, 1, points, weights)
sub_F_y[2,1]=gauss_3_4(sub_F_y_21, -1, 1, points, weights)
sub_F_y[2,2]=gauss_3_4(sub_F_y_22, -1, 1, points, weights)
sub_F_y[2,3]=gauss_3_4(sub_F_y_23, -1, 1, points, weights)

sub_F_y[3,0]=gauss_3_4(sub_F_y_30, -1, 1, points, weights)
sub_F_y[3,1]=gauss_3_4(sub_F_y_31, -1, 1, points, weights)
sub_F_y[3,2]=gauss_3_4(sub_F_y_32, -1, 1, points, weights)
sub_F_y[3,3]=gauss_3_4(sub_F_y_33, -1, 1, points, weights)

sub_F = sub_F_x + sub_F_y

sub_F = np.zeros((2*N_e_r+2 , 2*N_e_r+2))          # initialisation of global flux matrix
sub_F[0,0] = sub_F_x[0,0] + sub_F_y[0,0]
sub_F[0,1] = sub_F_x[0,1] + sub_F_y[0,1]
sub_F[1,0] = sub_F_x[1,0] + sub_F_y[1,0]
sub_F[1,1] = sub_F_x[1,1] + sub_F_y[1,1]

sub_F[0,-1] = sub_F_x[0,3] + sub_F_y[0,3]
sub_F[0,-2] = sub_F_x[0,2] + sub_F_y[0,2]
sub_F[1,N_e_r*2+1] = sub_F_x[1,3] + sub_F_y[1,3]
sub_F[1,N_e_r*2] = sub_F_x[1,2] + sub_F_y[1,2]


sub_F[-2,0] = sub_F_x[2,0] + sub_F_y[2,0]
sub_F[-1,0] = sub_F_x[3,0] + sub_F_y[3,0]
sub_F[-2,1] = sub_F_x[2,1] + sub_F_y[2,1]
sub_F[-1,1] = sub_F_x[3,1] + sub_F_y[3,1]

sub_F[-2,-2] = sub_F_x[2,2] + sub_F_y[2,2]
sub_F[-1,-2] = sub_F_x[3,2] + sub_F_y[3,2]
sub_F[-2,-1] = sub_F_x[2,3] + sub_F_y[2,3]
sub_F[-1,-1] = sub_F_x[3,3] + sub_F_y[3,3]



##############################################################################

# RHS_cst = (M + K - F)
# for n in range(nt):                 # Marching in time
#     Un = U.copy()
#     RHS = RHS_cst.dot(Un)           # saving U^t to be used at the next timestep calculation
#     U=np.linalg.solve(M,RHS)        # solving for U(t+1)
    
#     if n==1:                        # saving U at timestep 1 to plot
#         U1=U
#     elif n==nt//2:                  # saving U at half timestep to plot
#         U2=U


# single integral with Gaussia Quadrature
# import numpy as np

# E = np.array([-0.774597, 0.000000, 0.774597])
# A = np.array([0.555556, 0.888889, 0.555556])

# def gauss(f, a, b, E, A):
#     x = np.zeros(3)
#     for i in range(3):
#         x[i] = (b+a)/2 + (b-a)/2 *E[i]

#     return (b-a)/2 * (A[0]*f(x[0]) + A[1]*f(x[1]) + A[2]*f(x[2]))


# f = lambda x: 2 + np.sin(x)
# a = 0.0; b = np.pi/2

# areaGau = gauss(f, a, b, E, A)
# print("Gaussian integral: ", areaGau)

# sub_F_00 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1-y)*det*-c_x*dt + 1/4*(1-x)*(1-y)*1/4*(1-x)*(1-y)*det*-c_y*dt
# sub_F_01 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1-y)*det*-c_x*dt + 1/4*(1-x)*(1-y)*1/4*(1+x)*(1-y)*det*-c_y*dt
# sub_F_02 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1+y)*det*c_x*dt + 1/4*(1-x)*(1-y)*1/4*(1-x)*(1+y)*det*c_y*dt
# sub_F_03 = lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1+y)*det*c_x*dt + 1/4*(1-x)*(1-y)*1/4*(1+x)*(1+y)*det*c_y*dt

# sub_F_10 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1-y)*det*-c_x*dt + 1/4*(1+x)*(1-y)*1/4*(1-x)*(1-y)*det*-c_y*dt
# sub_F_11 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1-y)*det*-c_x*dt + 1/4*(1+x)*(1-y)*1/4*(1+x)*(1-y)*det*-c_y*dt
# sub_F_12 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1+y)*det*c_x*dt + 1/4*(1+x)*(1-y)*1/4*(1-x)*(1+y)*det*c_y*dt
# sub_F_13 = lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1+y)*det*c_x*dt + 1/4*(1+x)*(1-y)*1/4*(1+x)*(1+y)*det*c_y*dt

# sub_F_20 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1-y)*det*-c_x*dt + 1/4*(1-x)*(1+y)*1/4*(1-x)*(1-y)*det*-c_y*dt
# sub_F_21 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1-y)*det*-c_x*dt + 1/4*(1-x)*(1+y)*1/4*(1+x)*(1-y)*det*-c_y*dt
# sub_F_22 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1+y)*det*c_x*dt + 1/4*(1-x)*(1+y)*1/4*(1-x)*(1+y)*det*c_y*dt
# sub_F_23 = lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1+y)*det*c_x*dt + 1/4*(1-x)*(1+y)*1/4*(1+x)*(1+y)*det*c_y*dt

# sub_F_30 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1-y)*det*-c_x*dt + 1/4*(1+x)*(1+y)*1/4*(1-x)*(1-y)*det*-c_y*dt
# sub_F_31 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1-y)*det*-c_x*dt + 1/4*(1+x)*(1+y)*1/4*(1+x)*(1-y)*det*-c_y*dt
# sub_F_32 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1+y)*det*c_x*dt + 1/4*(1+x)*(1+y)*1/4*(1-x)*(1+y)*det*c_y*dt
# sub_F_33 = lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1+y)*det*c_x*dt + 1/4*(1+x)*(1+y)*1/4*(1+x)*(1+y)*det*c_y*dt

# sub_F=np.zeros((Np,Np))
# sub_F[0,0]=gauss(sub_F_00, -1, 1, E, A)
# sub_F[0,1]=gauss(sub_F_01, -1, 1, E, A)
# sub_F[0,2]=gauss(sub_F_02, -1, 1, E, A)
# sub_F[0,3]=gauss(sub_F_03, -1, 1, E, A)

# sub_F[1,0]=gauss(sub_F_10, -1, 1, E, A)
# sub_F[1,1]=gauss(sub_F_11, -1, 1, E, A)
# sub_F[1,2]=gauss(sub_F_12, -1, 1, E, A)
# sub_F[1,3]=gauss(sub_F_13, -1, 1, E, A)

# sub_F[2,0]=gauss(sub_F_20, -1, 1, E, A)
# sub_F[2,1]=gauss(sub_F_21, -1, 1, E, A)
# sub_F[2,2]=gauss(sub_F_22, -1, 1, E, A)
# sub_F[2,3]=gauss(sub_F_23, -1, 1, E, A)

# sub_F[3,0]=gauss(sub_F_30, -1, 1, E, A)
# sub_F[3,1]=gauss(sub_F_31, -1, 1, E, A)
# sub_F[3,2]=gauss(sub_F_32, -1, 1, E, A)
# sub_F[3,3]=gauss(sub_F_33, -1, 1, E, A)
