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
nt = 30 
N_e_r = 4  
N_e_c = 3
#nx =  N_e_c * N_e_r * 2                     # total number of x steps               
#ny =  N_e_c * N_e_r * 2                     # total number of y steps                                                                 
N = N_e_c * N_e_r                           # number of elements
L = 0.5                                       # x and y lengths
dx = L/(N_e_r)
dy = L/(N_e_c)
c_x = 0                                   # Wave velocity in x-dir
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
x, y, x_m, y_m =sy.symbols('x y x_m y_m')
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
sub_M=np.eye(4)*gauss(lambda x,y:1/4*(1-x)*(1-y)*1/4*(1-x)*(1-y)*det, -1, 1, points, weights)
sub_M[0,1] = gauss(lambda x,y: 1/4*(1-x)*(1-y) * 1/4*(1+x)*(1-y)*det, -1, 1, points, weights)
sub_M[0,2] = gauss(lambda x,y: 1/4*(1-x)*(1-y) * 1/4*(1-x)*(1+y)*det, -1, 1, points, weights)
sub_M[0,3] = gauss(lambda x,y: 1/4*(1-x)*(1-y) * 1/4*(1+x)*(1+y)*det, -1, 1, points, weights)

sub_M[1,0] = gauss(lambda x,y: 1/4*(1+x)*(1-y) * 1/4*(1-x)*(1-y)*det, -1, 1, points, weights)
sub_M[1,2] = gauss(lambda x,y: 1/4*(1+x)*(1-y) * 1/4*(1-x)*(1+y)*det, -1, 1, points, weights)
sub_M[1,3] = gauss(lambda x,y: 1/4*(1+x)*(1-y) * 1/4*(1+x)*(1+y)*det, -1, 1, points, weights)

sub_M[2,0] = gauss(lambda x,y: 1/4*(1-x)*(1+y) * 1/4*(1-x)*(1-y)*det, -1, 1, points, weights)
sub_M[2,1] = gauss(lambda x,y: 1/4*(1-x)*(1+y) * 1/4*(1+x)*(1-y)*det, -1, 1, points, weights)
sub_M[2,3] = gauss(lambda x,y: 1/4*(1-x)*(1+y) * 1/4*(1+x)*(1+y)*det, -1, 1, points, weights)

sub_M[3,0] = gauss(lambda x,y: 1/4*(1+x)*(1+y) * 1/4*(1-x)*(1-y)*det, -1, 1, points, weights)
sub_M[3,1] = gauss(lambda x,y: 1/4*(1+x)*(1+y) * 1/4*(1+x)*(1-y)*det, -1, 1, points, weights)
sub_M[3,2] = gauss(lambda x,y: 1/4*(1+x)*(1+y) * 1/4*(1-x)*(1+y)*det, -1, 1, points, weights)

M = np.kron(np.eye(N), sub_M)                                       # Creating global mass matrix

############################## stiffness matrix K ############################
sub_K=np.zeros((4,4))
sub_K[0,0]=gauss(lambda x,y: c_x*dt * -1/4*(1-y) * 1/4*(1-x)*(1-y) * dy/2 + c_y*dt *-1/4*(1-x) * 1/4*(1-x)*(1-y) * dx/2, -1, 1, points, weights)
sub_K[0,1]=gauss(lambda x,y: c_x*dt * -1/4*(1-y) * 1/4*(1+x)*(1-y) * dy/2 + c_y*dt *-1/4*(1-x) * 1/4*(1+x)*(1-y) * dx/2, -1, 1, points, weights)
sub_K[0,2]=gauss(lambda x,y: c_x*dt * -1/4*(1-y) * 1/4*(1-x)*(1+y) * dy/2 + c_y*dt *-1/4*(1-x) * 1/4*(1-x)*(1+y) * dx/2, -1, 1, points, weights)
sub_K[0,3]=gauss(lambda x,y: c_x*dt * -1/4*(1-y) * 1/4*(1+x)*(1+y) * dy/2 + c_y*dt *-1/4*(1-x) * 1/4*(1+x)*(1+y) * dx/2, -1, 1, points, weights)

sub_K[1,0]=gauss(lambda x,y: c_x*dt * 1/4*(1-y) * 1/4*(1-x)*(1-y) * dy/2 + c_y*dt *-1/4*(1+x) * 1/4*(1-x)*(1-y) * dx/2, -1, 1, points, weights)
sub_K[1,1]=gauss(lambda x,y: c_x*dt * 1/4*(1-y) * 1/4*(1+x)*(1-y) * dy/2 + c_y*dt *-1/4*(1+x) * 1/4*(1+x)*(1-y) * dx/2, -1, 1, points, weights)
sub_K[1,2]=gauss(lambda x,y: c_x*dt * 1/4*(1-y) * 1/4*(1-x)*(1+y) * dy/2 + c_y*dt *-1/4*(1+x) * 1/4*(1-x)*(1+y) * dx/2, -1, 1, points, weights)
sub_K[1,3]=gauss(lambda x,y: c_x*dt * 1/4*(1-y) * 1/4*(1+x)*(1+y) * dy/2 + c_y*dt *-1/4*(1+x) * 1/4*(1+x)*(1+y) * dx/2, -1, 1, points, weights)

sub_K[2,0]=gauss(lambda x,y: c_x*dt * -1/4*(1+y) * 1/4*(1-x)*(1-y) * dy/2 + c_y*dt *1/4*(1-x) * 1/4*(1-x)*(1-y) * dx/2, -1, 1, points, weights)
sub_K[2,1]=gauss(lambda x,y: c_x*dt * -1/4*(1+y) * 1/4*(1+x)*(1-y) * dy/2 + c_y*dt *1/4*(1-x) * 1/4*(1+x)*(1-y) * dx/2, -1, 1, points, weights)
sub_K[2,2]=gauss(lambda x,y: c_x*dt * -1/4*(1+y) * 1/4*(1-x)*(1+y) * dy/2 + c_y*dt *1/4*(1-x) * 1/4*(1-x)*(1+y) * dx/2, -1, 1, points, weights)
sub_K[2,3]=gauss(lambda x,y: c_x*dt * -1/4*(1+y) * 1/4*(1+x)*(1+y) * dy/2 + c_y*dt *1/4*(1-x) * 1/4*(1+x)*(1+y) * dx/2, -1, 1, points, weights)

sub_K[3,0]=gauss(lambda x,y: c_x*dt * 1/4*(1+y) * 1/4*(1-x)*(1-y) * dy/2 + c_y*dt *1/4*(1+x) * 1/4*(1-x)*(1-y) * dx/2, -1, 1, points, weights)
sub_K[3,1]=gauss(lambda x,y: c_x*dt * 1/4*(1+y) * 1/4*(1+x)*(1-y) * dy/2 + c_y*dt *1/4*(1+x) * 1/4*(1+x)*(1-y) * dx/2, -1, 1, points, weights)
sub_K[3,2]=gauss(lambda x,y: c_x*dt * 1/4*(1+y) * 1/4*(1-x)*(1+y) * dy/2 + c_y*dt *1/4*(1+x) * 1/4*(1-x)*(1+y) * dx/2, -1, 1, points, weights)
sub_K[3,3]=gauss(lambda x,y: c_x*dt * 1/4*(1+y) * 1/4*(1+x)*(1+y) * dy/2 + c_y*dt *1/4*(1+x) * 1/4*(1+x)*(1+y) * dx/2, -1, 1, points, weights)

K = np.kron(np.eye(N), sub_K)            # Creating global stiffness matrix

################################ flux #########################################
# x-dir line connecting points 1 and 3 with positive upwind flux n_x = -1
def gauss_1_3(f, a, b, points, weights):
    y = np.zeros(3)
    for i in range(3):
        y[i] = (b+a)/2 + (b-a)/2 *points[i]
        x=-1

    return (b-a)/2 * (weights[0]*f(x,y[0]) + weights[1]*f(x,y[1]) + weights[2]*f(x,y[2]))


sub_F_x=np.zeros((Np,Np))

sub_F_x[0,0]=gauss_1_3(lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1-y)*dy/2*c_x*dt*(-1), -1, 1, points, weights)
# sub_F_x[0,1]=gauss_1_3(lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1-y)*dy/2*c_x*dt*(-1), -1, 1, points, weights)
sub_F_x[0,2]=gauss_1_3(lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1+y)*dy/2*c_x*dt*(-1), -1, 1, points, weights)
# sub_F_x[0,3]=gauss_1_3(lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1+y)*dy/2*c_x*dt*(-1), -1, 1, points, weights)

sub_F_x[2,0]=gauss_1_3(lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1-y)*dy/2*c_x*dt*(-1), -1, 1, points, weights)
sub_F_x[2,1]=gauss_1_3(lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1-y)*dy/2*c_x*dt*(-1), -1, 1, points, weights)
sub_F_x[2,2]=gauss_1_3(lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1+y)*dy/2*c_x*dt*(-1), -1, 1, points, weights)
sub_F_x[2,3]=gauss_1_3(lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1+y)*dy/2*c_x*dt*(-1), -1, 1, points, weights)

# x-dir line connecting points 2 and 4 with positive upwind flux n_x = +1
def gauss_2_4(f, a, b, points, weights):
    y = np.zeros(3)
    for i in range(3):
        y[i] = (b+a)/2 + (b-a)/2 *points[i]
        x=1

    return (b-a)/2 * (weights[0]*f(x,y[0]) + weights[1]*f(x,y[1]) + weights[2]*f(x,y[2]))


# sub_F_x[1,0]=gauss_2_4(lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1-y)*dy/2*c_x*dt*(+1), -1, 1, points, weights)
sub_F_x[1,1]=gauss_2_4(lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1-y)*dy/2*c_x*dt*(+1), -1, 1, points, weights)
# sub_F_x[1,2]=gauss_2_4(lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1+y)*dy/2*c_x*dt*(+1), -1, 1, points, weights)
sub_F_x[1,3]=gauss_2_4(lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1+y)*dy/2*c_x*dt*(+1), -1, 1, points, weights)

sub_F_x[3,0]=gauss_2_4(lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1-y)*dy/2*c_x*dt*(+1), -1, 1, points, weights)
sub_F_x[3,1]=gauss_2_4(lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1-y)*dy/2*c_x*dt*(+1), -1, 1, points, weights)
sub_F_x[3,2]=gauss_2_4(lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1+y)*dy/2*c_x*dt*(+1), -1, 1, points, weights)
sub_F_x[3,3]=gauss_2_4(lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1+y)*dy/2*c_x*dt*(+1), -1, 1, points, weights)

################################## y -dir #####################################
# y-dir line connecting points 1 weightsnd 2 with positive upwind flux n_y = -1
def gauss_1_2(f, a, b, points, weights):
    x = np.zeros(3)
    for i in range(3):
        x[i] = (b+a)/2 + (b-a)/2 *points[i]
        y=-1

    return (b-a)/2 * (weights[0]*f(x[0],y) + weights[1]*f(x[1],y) + weights[2]*f(x[2],y))


sub_F_y=np.zeros((Np,Np))
sub_F_y[0,0]=gauss_1_2(lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1-y)*det*c_y*dt*(-1), -1, 1, points, weights)
sub_F_y[0,1]=gauss_1_2(lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1-y)*det*c_y*dt*(-1), -1, 1, points, weights)
sub_F_y[0,2]=gauss_1_2(lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1-x)*(1+y)*det*c_y*dt*(-1), -1, 1, points, weights)
sub_F_y[0,3]=gauss_1_2(lambda x,y: 1/4*(1-x)*(1-y)*1/4*(1+x)*(1+y)*det*c_y*dt*(-1), -1, 1, points, weights)

sub_F_y[1,0]=gauss_1_2(lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1-y)*det*c_y*dt*(-1), -1, 1, points, weights)
sub_F_y[1,1]=gauss_1_2(lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1-y)*det*c_y*dt*(-1), -1, 1, points, weights)
sub_F_y[1,2]=gauss_1_2(lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1-x)*(1+y)*det*c_y*dt*(-1), -1, 1, points, weights)
sub_F_y[1,3]=gauss_1_2(lambda x,y: 1/4*(1+x)*(1-y)*1/4*(1+x)*(1+y)*det*c_y*dt*(-1), -1, 1, points, weights)


# y-dir line connecting points 3 and 4 with positive upwind flux n_y = +1
def gauss_3_4(f, a, b, points, weights):
    x = np.zeros(3)
    for i in range(3):
        x[i] = (b+a)/2 + (b-a)/2 *points[i]
        y=+1

    return (b-a)/2 * (weights[0]*f(x[0],y) + weights[1]*f(x[1],y) + weights[2]*f(x[2],y))


sub_F_y[2,0]=gauss_3_4(lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1-y)*dx/2*c_y*dt*(+1), -1, 1, points, weights)
sub_F_y[2,1]=gauss_3_4(lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1-y)*dx/2*c_y*dt*(+1), -1, 1, points, weights)
sub_F_y[2,2]=gauss_3_4(lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1-x)*(1+y)*dx/2*c_y*dt*(+1), -1, 1, points, weights)
sub_F_y[2,3]=gauss_3_4(lambda x,y: 1/4*(1-x)*(1+y)*1/4*(1+x)*(1+y)*dx/2*c_y*dt*(+1), -1, 1, points, weights)

sub_F_y[3,0]=gauss_3_4(lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1-y)*dx/2*c_y*dt*(+1), -1, 1, points, weights)
sub_F_y[3,1]=gauss_3_4(lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1-y)*dx/2*c_y*dt*(+1), -1, 1, points, weights)
sub_F_y[3,2]=gauss_3_4(lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1-x)*(1+y)*dx/2*c_y*dt*(+1), -1, 1, points, weights)
sub_F_y[3,3]=gauss_3_4(lambda x,y: 1/4*(1+x)*(1+y)*1/4*(1+x)*(1+y)*dx/2*c_y*dt*(+1), -1, 1, points, weights)

sub_F = sub_F_x + sub_F_y

sub_F2_x = np.zeros((2*N_e_r+2 , Np*N_e_r+2))          # initialisation of global flux 
sub_F2_x[-2,2*N_e_r] = sub_F2_x[0,-2] = sub_F_x[0,2]
sub_F2_x[-1,2*N_e_r+1] = sub_F2_x[1,-1] = sub_F_x[1,3]
sub_F2_x[-2,Np*N_e_r-1] = sub_F2_x[0,2*N_e_r-1] = sub_F_x[0,0]
sub_F2_x[2*N_e_r+1 , Np*N_e_r+1] = sub_F2_x[1,2*N_e_r+1] = sub_F_x[1,1] 

sub_F2_y = np.zeros((2*N_e_r+2 , Np*N_e_r+2))          # initialisation of global flux 
sub_F2_y[0,0] = sub_F2_y[1,1] = sub_F_y[0,0]
sub_F2_y[0,1] = sub_F2_y[1,0] = sub_F_y[0,1]    
sub_F2_y[0,2*N_e_r+1] = sub_F2_y[1,2*N_e_r] =  sub_F_y[2,3]
sub_F2_y[2*N_e_r+1 , Np*N_e_r+1] = sub_F2_y[-2 , -2] = sub_F_y[3,3]

sub_F2 = sub_F2_x + sub_F2_y

F = np.zeros((Np*N_e_r*N_e_c , Np*N_e_r*N_e_c+2*N_e_r))          # initialisation of global flux matrix
i=0
j=0

for n in range(N_e_c):
    for n in range(N_e_r):
        F[i:i+sub_F2.shape[0], j:j+sub_F2.shape[1]]+=sub_F2       # copying local element into the global matrix
        i+= 2
        j+= 2
    i+=N_e_r*2
    j+=N_e_r*2

i=0
j=N_e_r*2-1
while i<Np*N_e_r*N_e_c:
    F[i,j]=0
    i+=N_e_r*2
    j+=N_e_r*2
F= F[:,2*N_e_r:]

#plt.spy(F)
##############################################################################

RHS_cst = (M + K - F)
for n in range(nt):                 # Marching in time
    Un = U.copy()
    RHS = RHS_cst.dot(Un)           # saving U^t to be used at the next timestep calculation
    U=np.linalg.solve(M,RHS)        # solving for U(t+1)
    
    if n==1:                        # saving U at timestep 1 to plot
        U1=U
    elif n==nt//2:                  # saving U at half timestep to plot
        U2=U






import matplotlib.pyplot as plt
import matplotlib.collections
import numpy as np


def showMeshPlot(nodes, elements, values):

    y = nodes[:,0]
    z = nodes[:,1]

    def quatplot(y,z, quatrangles, values, ax=None, **kwargs):

        if not ax: ax=plt.gca()
        yz = np.c_[y,z]
        verts= yz[quatrangles]
        pc = matplotlib.collections.PolyCollection(verts, **kwargs)
        pc.set_array(values)
        ax.add_collection(pc)
        ax.autoscale()
        return pc

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    pc = quatplot(y,z, np.asarray(elements), values, ax=ax, 
             edgecolor="crimson", cmap="rainbow")
    fig.colorbar(pc, ax=ax)        
    ax.plot(y,z, marker="o", ls="", color="crimson")

    ax.set(title='This is the plot for: quad', xlabel='Y Axis', ylabel='Z Axis')

    plt.show()

nodes = np.array([[0,0],    [0,dx],    [0,dx],    [0,2*dx],    [0,2*dx],    [0,3*dx],    [0,3*dx],    [0,4*dx],
               [dy,0],   [dy,dx],   [dy,dx],   [dy,2*dx],   [dy,2*dx],   [dy,3*dx],   [dy,3*dx],   [dy,4*dx],
               [dy,0], [dy,dx], [dy,dx], [dy,2*dx], [dy,2*dx], [dy,3*dx], [dy,3*dx], [dy,4*dx],
               [2*dy,0], [2*dy,dx], [2*dy,dx], [2*dy,2*dx], [2*dy,2*dx], [2*dy,3*dx], [2*dy,3*dx], [2*dy,4*dx],
               [2*dy,0], [2*dy,dx], [2*dy,dx], [2*dy,2*dx], [2*dy,2*dx], [2*dy,3*dx], [2*dy,3*dx], [2*dy,4*dx],
               [3*dy,0], [3*dy,dx], [3*dy,dx], [3*dy,2*dx], [3*dy,2*dx], [3*dy,3*dx], [3*dy,3*dx], [3*dy,4*dx]])

elements = np.array([[0,1,8,9],[2,3,10,11], [4,5,12,13],[6,7,14,15], [16,17,24,25], [18,19,26,27],
                     [20,21,28,29], [22,23,30,31], [32,33,40,41], [34,35,42,43], [36,37,44,45], [38,39,46,47]])

showMeshPlot(nodes, elements, U)






# import matplotlib.pyplot as plt
# import numpy as np
# import matplotlib.collections
# dx=0.125
# dy=0.1667
# xy = np.array([[0,0],    [0,dx],    [0,dx],    [0,2*dx],    [0,2*dx],    [0,3*dx],    [0,3*dx],    [0,4*dx],
#                [dy,0],   [dy,dx],   [dy,dx],   [dy,2*dx],   [dy,2*dx],   [dy,3*dx],   [dy,3*dx],   [dy,4*dx],
#                [dy,0], [dy,dx], [dy,dx], [dy,2*dx], [dy,2*dx], [dy,3*dx], [dy,3*dx], [dy,4*dx],
#                [2*dy,0], [2*dy,dx], [2*dy,dx], [2*dy,2*dx], [2*dy,2*dx], [2*dy,3*dx], [2*dy,3*dx], [2*dy,4*dx],
#                [2*dy,0], [2*dy,dx], [2*dy,dx], [2*dy,2*dx], [2*dy,2*dx], [2*dy,3*dx], [2*dy,3*dx], [2*dy,4*dx],
#                [3*dy,0], [3*dy,dx], [3*dy,dx], [3*dy,2*dx], [3*dy,2*dx], [3*dy,3*dx], [3*dy,3*dx], [3*dy,4*dx]])

# elements = np.array([[0,1,8,9],[2,3,10,11], [4,5,12,13],[6,7,14,15], [16,17,24,25], [18,19,26,27],
#                      [20,21,28,29], [22,23,30,31], [32,33,40,41], [34,35,42,43], [36,37,44,45], [38,39,46,47]])

# x = np.degrees(xy[:, 0])
# y = np.degrees(xy[:, 1])

# def quatplot(x,y, quatrangles, ax=None, **kwargs):
#     if not ax: ax=plt.gca()
#     xy = np.c_[x,y]
#     verts=xy[quatrangles]
#     pc = matplotlib.collections.PolyCollection(verts, **kwargs)
#     ax.add_collection(pc)
#     ax.autoscale()


# plt.figure()
# plt.gca().set_aspect('equal')

# quatplot(x,y, elements, ax=None, color="crimson", facecolor="None")
# plt.plot(x,y, marker="o", ls="", color="crimson")

# plt.xlabel('x')
# plt.ylabel('y')
# for i, (xi,yi) in enumerate(np.degrees(xy)):
#     plt.text(yi,xi,i, size=8)

# plt.show()









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
