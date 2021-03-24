## 2D new method
#points, weights = np.polynomial.legendre.leggauss(3)
#points.shape = [points.shape[0], 1]
import numpy as np
import matplotlib.pyplot as plt
quadrature_points, weights = np.polynomial.legendre.leggauss(3)
C = 0.05
c_x = 0.1
c_y = 0.1
dx = 0.125
dy = 0.1667
dt = C*dx*dy/(c_x*dy+c_y*dx)
N_e_r = 4
N_e_c= 3
local_node_no = 4
total_element=12
total_sloc = 4
total_nodes = total_element * local_node_no
no_of_qp = 9 # degree of polynomial**2
M = np.zeros((total_nodes, total_nodes))
K = np.zeros((total_nodes, total_nodes))
F0 = np.zeros((total_nodes, total_nodes))

#xi =  [-0.7745967,          0,  0.7745967, -0.7745967,   0, 0.7745967, -0.7745967,         0, 0.7745967]
#eta = [-0.7745967, -0.7745967, -0.7745967,          0,   0,         0,  0.7745967, 0.7745967, 0.7745967]
#w_xi = [      5/9,        8/9,        5/9,        5/9, 8/9,       5/9,        5/9,       8/9,       5/9]
#w_eta= [      5/9,        5/9,        5/9,        8/9, 8/9,       8/9,        5/9,       5/9,       5/9]
#-------------------- shape func, ddx and ddy ---------------------------------
shape_func = {0:lambda xi,eta: 1/4*(1-xi)*(1-eta),
              1:lambda xi,eta: 1/4*(1+xi)*(1-eta),
              2:lambda xi,eta: 1/4*(1-xi)*(1+eta),
              3:lambda xi,eta: 1/4*(1+xi)*(1+eta)}

ddx_shape_func = {0:lambda xi,eta: -1/4*(1-eta),
                  1:lambda xi,eta:  1/4*(1-eta),
                  2:lambda xi,eta: -1/4*(1+eta),
                  3:lambda xi,eta:  1/4*(1+eta)}

ddy_shape_func = {0:lambda xi,eta: -1/4*(1-xi),
                  1:lambda xi,eta: -1/4*(1+xi),
                  2:lambda xi,eta:  1/4*(1-xi),
                  3:lambda xi,eta:  1/4*(1+xi)}

p=np.array([1,2,3,4]).reshape(4,1)

#------------------------------ global node numbering -------------------------
def global_no(e,N_e_r):
    col=int(np.ceil(e/N_e_r))
    glob_no = np.array([(col-1)*2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*(e-1)+2, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+2])
    return glob_no


#global_no(7,4)
loc_to_glob_list=[]
for i in range(1,total_element+1):
    loc_to_glob_list.append(global_no(i,N_e_r))

loc_to_glob_list[0][2]
#------------------------------ Main structure of the code --------------------
for element_no in range (total_element):    # elementsnumbers starts from 1 
    def coordinates(element_no):            # global node coordinates
        col =int(np.ceil(element_no/N_e_r))
        row = int(element_no - (col -1) * N_e_r)
        co_ordinates = np.array(([dx*(row-1), dy*(col-1)],
                                    [dx*row  , dy*(col-1)],
                                    [dx*(row-1), dy*(col)],
                                    [dx*row  , dy*(col)]))
        return co_ordinates


    def jacobian(xi, eta, element_no):
        a = 1/4 * np.matrix(([-(1-eta), 1-eta, -(1+eta), 1+eta],
                                   [-(1-xi), -(1+xi), 1-xi, 1+xi])) 
        jacobi = a.dot(coordinates(element_no))
        return jacobi

    
    def det(jacobian):
        return np.linalg.det(jacobian)
    
    
    def inverse(jacobian):
        return np.linalg.inv(jacobian)
    
    
    for iloc in range(local_node_no):     
        global_i = loc_to_glob_list[element_no][iloc]
        for jloc in range(local_node_no):   
            global_j = loc_to_glob_list[element_no][jloc]           
#           for g in range(no_of_qp):
            def gauss(f, quadrature_points, weights):
                xi=eta = np.zeros(len(quadrature_points))
                for i in range(len(quadrature_points)):
                    xi[i] =quadrature_points[i]
                    eta[i] =quadrature_points[i]
                answer = 0
                    
                for i in range(len(quadrature_points)):
                    for j in range(len(quadrature_points)):     
                        answer =  answer + weights[i] * weights[j] * f(xi[i], eta[j]) * det(jacobian(xi[i], eta[j], element_no))
                return answer
                
            def ff(xi, eta):
                return shape_func[iloc](xi,eta)*shape_func[jloc](xi,eta)    
            M[global_i-1,global_j-1] = gauss(lambda xi,eta: shape_func[iloc](xi,eta)*shape_func[jloc](xi,eta), quadrature_points, weights)
            
            def K_gauss(f, quadrature_points, weights):
                xi=eta = np.zeros(len(quadrature_points))
                for i in range(len(quadrature_points)):
                    xi[i] =quadrature_points[i]
                    eta[i] =quadrature_points[i]
                answer = 0
                    
                for i in range(len(quadrature_points)):
                    for j in range(len(quadrature_points)):     
                        answer =  answer + weights[i] * weights[j] * f(xi[i], eta[j])* det(jacobian(xi[i], eta[j], element_no))
                return answer
            K[global_i-1,global_j-1] = K_gauss(lambda xi,eta: c_x*dt*shape_func[jloc](xi,eta)*ddx_shape_func[iloc](xi,eta)* inverse(jacobian(xi,eta,element_no))[0,0]+
                                                              c_y*dt*shape_func[jloc](xi,eta)*ddy_shape_func[iloc](xi,eta)*inverse(jacobian(xi,eta,element_no))[1,1], quadrature_points, weights)
# plt.spy(K)            
    
#------------------------- surface integration ---------------------------------
#  do iface =1,nface
#! for surface iface...
#      do siloc=1,snloc ! local row siloc for surface element
#          sinod=   ! global node number
#          do sjloc=1,snloc ! local coln sjloc for surface element
#             sjnod=   ! global node number
#
#               do sgi=1,sngi ! loop over the quadrature points on surface
# L_quadrature_points, L_weights = np.polynomial.legendre.leggauss(2) 
x_quadrature_points = np.array([-0.5, 0.5,  1  ,  1  , 0.5, -0.5, -1  , -1  ])
y_quadrature_points = np.array([-1  ,-1  , -0.5,  0.5, 1  ,  1  ,  0.5, -0.5])
L_weights=np.array([1,1,1,1,1,1,1,1])

domain_norm = [0,0,1]

dx_dxi = lambda eta: 1/4*((eta-1)*coordinates(e)[0,0] +(1-eta)*coordinates(e)[1,0] -(1+eta)*coordinates(e)[2,0] +(1+eta)*coordinates(e)[3,0])
dy_dxi = lambda eta: 1/4*((eta-1)*coordinates(e)[0,1] +(1-eta)*coordinates(e)[1,1] -(1+eta)*coordinates(e)[2,1] +(1+eta)*coordinates(e)[3,1])
dx_deta = lambda xi: 1/4*((xi-1)*coordinates(e)[0,0] -(xi+1)*coordinates(e)[1,0] +(1-xi)*coordinates(e)[2,0] +(1+xi)*coordinates(e)[3,0])
dy_deta = lambda xi: 1/4*((xi-1)*coordinates(e)[0,1] -(xi+1)*coordinates(e)[1,1] +(1-xi)*coordinates(e)[2,1] +(1+xi)*coordinates(e)[3,1])

for e in range(total_element):      # element numbering starts from 0
    
    # calculating the center of each element used to calculate the sign of normal
    def e_centre(e):
        x_centre =0
        y_centre =0
        for i in range(local_node_no):
            x_centre = x_centre + coordinates(e+1)[i,0]         
            y_centre = y_centre + coordinates(e+1)[i,1]

        z_centre = 0
        e_centre = [x_centre/local_node_no, y_centre/local_node_no, z_centre]
        return e_centre
    #---------------------------------------------------------------------------------
            
    for s in range(total_sloc):                     # number of surfaces = no.elements * no.faces in each e
        for siloc in range(total_sloc):             # = 4
            
            # function giving global node numbers of 2 points of each face------------
            def s_glob_node(e,siloc):                
                global_nod ={0: [loc_to_glob_list[e][0] , loc_to_glob_list[e][1]],
                             1: [loc_to_glob_list[e][1] , loc_to_glob_list[e][3]],
                             2: [loc_to_glob_list[e][3] , loc_to_glob_list[e][2]],
                             3: [loc_to_glob_list[e][2] , loc_to_glob_list[e][0]]}
                return global_nod[siloc]
            #-------------------------------------------------------------------------
            # Jacobian
            jac = {0: np.sqrt(dx_dxi(-1)**2 + dy_dxi(-1)**2),
                   1: np.sqrt(dx_deta(1)**2 + dy_deta(1)**2),
                   2: np.sqrt(dx_dxi(1)**2 + dy_dxi(1)**2),
                   3: np.sqrt(dx_deta(-1)**2 + dy_deta(-1)**2)}
            
            # normal to the boundary lines -------------------------------------------
            n_hat = {0: np.sign(np.cross([dx_dxi(-1),dy_dxi(-1),0] , domain_norm)[1]),                  # n_ds of the line (-1,-1) and (-1,1)
                     1: np.sign(np.cross([dx_deta(1),dy_deta(1),0] , domain_norm)[0]),                  # n_ds of the line (-1,1)  and (1,1)
                     2: np.sign(np.cross(domain_norm,                [dx_dxi(1),dy_dxi(1),0]))[1],      # n_ds of the line (1,1)   and (1,-1)
                     3: np.sign(np.cross(domain_norm,                [dx_deta(-1),dy_deta(-1),0]))[0]}  # n_ds of the line (-1,1)  and (-1,-1)
            
            #-------------------------------------------------------------------------
            
            sinod = s_glob_node(e,siloc)                # gives two nodes numbrs of each i-surface
            
            for sjloc in range(total_sloc):             # = 4
                sjnod = s_glob_node(e,sjloc)            # gives two nodes numbrs of each j-surface
          
                # boundary Gaussian integration --------------------------------------
                def L_gauss(f, x_quadrature_points, y_quadrature_points, L_weights):  
                    xi=eta = np.zeros(len(x_quadrature_points))
                    for i in range(len(x_quadrature_points)):
                        xi[i] =x_quadrature_points[i]
                        eta[i] =y_quadrature_points[i]
                    answer = 0
                    for i in range(len(x_quadrature_points)):
                        answer =  answer + L_weights[i] * f(xi[i], eta[i]) * n_hat[sjloc] *jac[siloc]
                    return answer
                #--------------------------------------------------------------------
                F0[s_glob_node(e,siloc)[0]-1, s_glob_node(e,sjloc)[0]-1] = F0[s_glob_node(e,siloc)[1]-1, s_glob_node(e,sjloc)[1]-1] = L_gauss(lambda xi,eta: dt * c_x * shape_func[siloc](xi,eta)*shape_func[sjloc](xi,eta), x_quadrature_points, y_quadrature_points, L_weights)
                # print(s_glob_node(e,siloc)[0]-1, s_glob_node(e,sjloc)[0]-1, s_glob_node(e,siloc)[1]-1, s_glob_node(e,sjloc)[1]-1)


# plt.spy(F0)













            # vector from centre to one node on a boundary line
            # r = {0: np.subtract([coordinates(e+1)[0,0], coordinates(e+1)[0,1],0] , e_centre(e)),
            #      1: np.subtract([coordinates(e+1)[1,0], coordinates(e+1)[1,1],0] , e_centre(e)),
            #      2: np.subtract([coordinates(e+1)[3,0], coordinates(e+1)[3,1],0] , e_centre(e)),
            #      3: np.subtract([coordinates(e+1)[2,0], coordinates(e+1)[2,1],0] , e_centre(e))}
            
            # sign (norm . r)
            # n_hat = np.sign(np.dot(s_norm[siloc], r[siloc]))
            # print(e, siloc, n_hat[siloc])














# e=1
# x_centre = 0
# y_centre = 0
# for i in range(len(coordinates(e))):
#     x_centre = x_centre + coordinates(e)[i][0] / len(coordinates(e))
#     y_centre = y_centre + coordinates(e)[i][1] / len(coordinates(e)) 
# e_centre = [x_centre, y_centre]
# print(e_centre)             

# r1=[3,0,0]
# r2=[0.5,0,1]
# print(np.dot(r2,r1))
# print(np.linalg.norm(np.cross(r1,r2)))
#print(np.sign(np.cross(r1,r2))[1])
#s_glob_node(1-1,0)
#s_glob_node(1-1,1)
#s_glob_node(1-1,2)
#s_glob_node(1-1,3)
#coordinates(1)
#
#e=0
#
#0.5*abs(np.linalg.norm([coordinates(e)[0]-coordinates(e)[2]]))
#  


#
##------------------------ J, det and coordinates ------------------------------
#class Element:
#    def jacobian(self, xi, eta, element_no):
#        self.xi = xi
#        self.eta = eta
##        self.coordinates = coordinates
#        a = 1/4 * np.matrix(([-(1-eta), 1-eta, -(1+eta), 1+eta],
#                               [-(1-xi), -(1+xi), 1-xi, 1+xi])) 
#        jacobi = a.dot(self.coordinates(element_no))
#        return jacobi
#    
#    def det(self, jacobian):
#        self.jacobian = jacobian
#        return np.linalg.det(jacobian)
#    
#    
#    def coordinates(self, element_no):
#        self.element_no = element_no
#        col =int(np.ceil(element_no/N_e_r))
#        row = int(element_no - (col -1) * N_e_r)
#        co_ordinates = np.matrix(([dx*(row-1), dy*(col-1)],
#                                [dx*row  , dy*(col-1)],
#                                [dx*(row-1), dy*(col)],
#                                [dx*row  , dy*(col)]))
#        return co_ordinates
#
#
#element = Element()
#a= element.coordinates(12)
#b=element.jacobian(-1,1,12)
#element.det(element.jacobian(-1,1,12))
#
#
#
#

#class Element:
##    def __init__(self, element_no):
##        self.element_no = element_no
##        self.col =int(np.ceil(self.element_no/N_e_r))`
#    def coordinates(self):
#        col =int(np.ceil(self.element_no/N_e_r))
#        row = int(self.element_no - (col -1) * N_e_r)
#        co_ordinates = np.matrix(([dx*(row-1), dy*(col-1)],
#                                [dx*row  , dy*(col-1)],
#                                [dx*(row-1), dy*(col)],
#                                [dx*row  , dy*(col)]))
#        return co_ordinates



############# function giving global coordinates of faces #####################
#            def s_glob_node(e,i_j_loc):
#                if i_j_loc==0:
#                    a=[loc_to_glob_list[e][0], loc_to_glob_list[e][1]]
#                elif i_j_loc==1:
#                    a=[loc_to_glob_list[e][1], loc_to_glob_list[e][3]]
#                elif i_j_loc==2:
#                    a=[loc_to_glob_list[e][3], loc_to_glob_list[e][2]]
#                else:
#                    a=[loc_to_glob_list[e][2], loc_to_glob_list[e][0]]
#                return a
###############################################################################




#loc_to_glob_list={}
#for i in range(1,total_element+1):
#    loc_to_glob_list[i]= (global_no(i,4))



#
#xi =  [-sqrt(0.6),          0,  sqrt(0.6), -sqrt(0.6),   0, sqrt(0.6), -sqrt(0.6),         0, sqrt(0.6)]
#eta = [-sqrt(0.6), -sqrt(0.6), -sqrt(0.6),          0,   0,         0,  sqrt(0.6), sqrt(0.6), sqrt(0.6)]
#w_xi = [      5/9,        8/9,        5/9,        5/9, 8/9,       5/9,        5/9,       8/9,       5/9]
#w_eta= [      5/9,        5/9,        5/9,        8/9, 8/9,       8/9,        5/9,       5/9,       5/9]