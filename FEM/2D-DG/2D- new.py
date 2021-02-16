## 2D new method
#points, weights = np.polynomial.legendre.leggauss(3)
#points.shape = [points.shape[0], 1]
import numpy as np
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
F = np.zeros((total_nodes, total_nodes))

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

#loc_to_glob_list[0][2]
#------------------------------ Main structure of the code --------------------
for element_no in range (total_element):    # elementsnumbers starts from 1 
    def coordinates(element_no):            # global node coordinates
        col =int(np.ceil(element_no/N_e_r))
        row = int(element_no - (col -1) * N_e_r)
        co_ordinates = np.matrix(([dx*(row-1), dy*(col-1)],
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
#spy(K)            

#------------------------- surface integration ---------------------------------
#  do iface =1,nface
#! for surface iface...
#      do siloc=1,snloc ! local row siloc for surface element
#          sinod=   ! global node number
#          do sjloc=1,snloc ! local coln sjloc for surface element
#             sjnod=   ! global node number
#
#               do sgi=1,sngi ! loop over the quadrature points on surface
L_quadrature_points, L_weights = np.polynomial.legendre.leggauss(2) 

for e in range(total_element):      # element numbering starts from 0
    for s in range(total_sloc):
        for siloc in range(total_sloc):
            def s_glob_node(e,i_j_loc):
                if i_j_loc==0:
                    a=[loc_to_glob_list[e][0], loc_to_glob_list[e][1]]
                elif i_j_loc==1:
                    a=[loc_to_glob_list[e][1], loc_to_glob_list[e][3]]
                elif i_j_loc==2:
                    a=[loc_to_glob_list[e][3], loc_to_glob_list[e][2]]
                else:
                    a=[loc_to_glob_list[e][2], loc_to_glob_list[e][0]]
                return a
            
            
            sinod = s_glob_node(e,siloc)
            
            for sjloc in range(total_sloc):
                sjnod = s_glob_node(e,jloc)
          
                def L_gauss(f, L_quadrature_points, L_weights):
                    answer = 0
                        
                    for i in range(len(L_quadrature_points)):
                        answer =  answer + L_weights[i] * f(L_quadrature_points[i]) * det(jacobian(L_quadrature_points[i], element_no))
                    return answer
            
            
            if siloc or sjloc == 0:
                eta = -1
                n_dline = -0.5*abs(np.linalg.norm(coordinates(e)[1]-coordinates(e)[0]))
            elif siloc or sjloc == 2:
                eta = 1
                n_dline = 0.5*abs(np.linalg.norm(coordinates(e)[1]-coordinates(e)[0]))
            elif siloc or sjloc == 1:
                xi = 1
                n_dline = 0.5*abs(np.linalg.norm(coordinates(e)[3]-coordinates(e)[1]))
            elif siloc or sjloc == 3:
                xi = -1
                n_dline = -0.5*abs(np.linalg.norm(coordinates(e)[3]-coordinates(e)[1]))
            
            
            F[global_i-1,global_j-1] = L_gauss(lambda xi,eta: shape_func[siloc](xi,eta)*shape_func[sjloc](xi,eta)*n_dline, L_quadrature_points, L_weights)    
                
                #print(sjloc, 0.5*abs(np.linalg.norm([coordinates(e)[siloc]-coordinates(e)[sjloc]])))


#r1=[3,0,0]
#r2=[0,0,1]
#print(np.linalg.norm(np.cross(r2,r1)))
#print(np.linalg.norm(np.cross(r1,r2)))
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
##        self.col =int(np.ceil(self.element_no/N_e_r))
##        self.row = int(self.element_no - (col -1) * N_e_r)
##        co_ordinates = np.matrix(([dx*(row-1), dy*(col-1)],
##                                [dx*row  , dy*(col-1)],
##                                [dx*(row-1), dy*(col)],
##                                [dx*row  , dy*(col)]))
#    def jacobian(self, xi, eta):
#        self.xi = xi
#        self.eta = eta
##        self.coordinates = coordinates
#        a = 1/4 * np.matrix(([-(1-eta), 1-eta, -(1+eta), 1+eta],
#                               [-(1-xi), -(1+xi), 1-xi, 1+xi])) 
#        jacobi = a.dot(self.coordinates)
#        return jacobi
#    
#    def det(self, jacobian):
#        self.jacobian = jacobian
#        return np.linalg.det(jacobian)
#    
#    
#    def coordinates(self):
#        col =int(np.ceil(self.element_no/N_e_r))
#        row = int(self.element_no - (col -1) * N_e_r)
#        co_ordinates = np.matrix(([dx*(row-1), dy*(col-1)],
#                                [dx*row  , dy*(col-1)],
#                                [dx*(row-1), dy*(col)],
#                                [dx*row  , dy*(col)]))
#        return co_ordinates









#loc_to_glob_list={}
#for i in range(1,total_element+1):
#    loc_to_glob_list[i]= (global_no(i,4))



#
#xi =  [-sqrt(0.6),          0,  sqrt(0.6), -sqrt(0.6),   0, sqrt(0.6), -sqrt(0.6),         0, sqrt(0.6)]
#eta = [-sqrt(0.6), -sqrt(0.6), -sqrt(0.6),          0,   0,         0,  sqrt(0.6), sqrt(0.6), sqrt(0.6)]
#w_xi = [      5/9,        8/9,        5/9,        5/9, 8/9,       5/9,        5/9,       8/9,       5/9]
#w_eta= [      5/9,        5/9,        5/9,        8/9, 8/9,       8/9,        5/9,       5/9,       5/9]