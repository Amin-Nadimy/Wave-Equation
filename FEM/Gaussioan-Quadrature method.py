import numpy as np
import scipy as sp

E = np.array([-0.774597, 0.000000, 0.774597])
A = np.array([0.555556, 0.888889, 0.555556])

def gauss(f, a, b, E, A):
    x = np.zeros(len(E))
    for i in range(len(E)):
        x[i] = (b+a)/2 + (b-a)/2 *E[i]

    return (b-a)/2 * (A[0]*f(x[0]) + A[1]*f(x[1]) + A[2]*f(x[2]))
dx = 0.001001

############################ Mass matrix #####################################
M_00 = lambda x: 1/2*(1-x) *1/2 *(1-x) *dx/2
M_01 = M_10 = lambda x: 1/2 * (1-x) * 1/2 * (1+x)* dx/2
M_11 = lambda x: 1/4 * (1+x)**2* dx/2

sub_M = np.zeros((2,2))
sub_M[0,0] = gauss(M_00, -1, 1, E, A)
sub_M[0,1] = sub_M[1,0] = gauss(M_01, -1, 1, E, A)
sub_M[1,1] = gauss(M_11, -1, 1, E, A)
print(sub_M)


############################### K matrix #####################################
K_00 = lambda x: 1/2*(1-x) *-1/2
K_01 = lambda x: 1/2*(1-x) *1/2
K_10 = lambda x: 1/2*(1+x) *-1/2
K_11 = lambda x: 1/2*(1+x) *1/2

sub_K = np.zeros((2,2))
sub_K[0,0] = gauss(K_00, -1, 1, E, A)
sub_K[0,1] = gauss(K_01, -1, 1, E, A)
sub_K[1,0] = gauss(K_10, -1, 1, E, A)
sub_K[1,1] = gauss(K_11, -1, 1, E, A)
print(sub_K)

############################# 3 ##############################################
U = 0
Un = 1
c = 0.1
dt = 0.0005
dx = 0.001001
j_1 = lambda x: -dx/16/c * ( ((x-1)/2/dt - c/dx)**2  /  ( ((x-1)/dt)**2 + (1/dx)**2) )
j_2 = lambda x:  (U/dt * (1/2*(1-x)) - Un/dt*(1/2*(1-x) + c/dx  ) )**2   *  -dx/16/c  /  ( ((U-Un)/dt * 1/2*(1-x))**2 + (Un/dx)**2 )       

a=gauss(j_1, -1, 1, E, A)
b=gauss(j_2, -1, 1, E, A)
print(a)
print(b)





# # exact solution with tolerant =0
# f = lambda x: 1/4 *(1-x)**2
# sp.integrate.quadrature(f, -1, 1.0)


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