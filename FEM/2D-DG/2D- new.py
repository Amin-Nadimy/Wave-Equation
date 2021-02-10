## 2D new method
#points, weights = np.polynomial.legendre.leggauss(3)
#points.shape = [points.shape[0], 1]
import numpy as np
quadrature_points = [-np.sqrt(0.6),0, np.sqrt(0.6)]
weights = [5/9, 8/9, 5/9]
dx = 0.125
dy = 0.1667
N_e_r = 4
N_e_c= 3
local_node_no = 4
total_element=12
total_nodes = total_element * local_node_no
no_of_qp = 9 # degree of polynomial**2
M = np.zeros((total_nodes, total_nodes))

xi =  [-0.7745967,          0,  0.7745967, -0.7745967,   0, 0.7745967, -0.7745967,         0, 0.7745967]
eta = [-0.7745967, -0.7745967, -0.7745967,          0,   0,         0,  0.7745967, 0.7745967, 0.7745967]
w_xi = [      5/9,        8/9,        5/9,        5/9, 8/9,       5/9,        5/9,       8/9,       5/9]
w_eta= [      5/9,        5/9,        5/9,        8/9, 8/9,       8/9,        5/9,       5/9,       5/9]
#-------------------- shape func, ddx and ddy ---------------------------------
shape_func = {1:lambda xi,eta: 1/4*(1-xi)*(1-eta),
              2:lambda xi,eta: 1/4*(1+xi)*(1-eta),
              3:lambda xi,eta: 1/4*(1-xi)*(1+eta),
              4:lambda xi,eta: 1/4*(1+xi)*(1+eta)}

ddx_shape_function = {1:lambda xi,eta: -1/4*(1-eta),
                      2:lambda xi,eta:  1/4*(1-eta),
                      3:lambda xi,eta: -1/4*(1+eta),
                      4:lambda xi,eta: -1/4*(1+eta)}

ddy_shape_function = {1:lambda xi,eta: -1/4*(1-xi),
                      2:lambda xi,eta: -1/4*(1+xi),
                      3:lambda xi,eta:  1/4*(1-xi),
                      4:lambda xi,eta:  1/4*(1+xi)}
#def func(xi,eta):
#    fn = shape_func[1](xi,eta)*shape_func[2](xi,eta)
#    return fn
#
#print(func(3,2))



#------------------------------ global node numbering -------------------------
def global_no(e,N_e_r):
    col=int(np.ceil(e/N_e_r))
    glob_no = np.array([(col-1)*2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*(e-1)+2, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+2])
    return glob_no


#global_no(7,4)
loc_to_glob_list=[]
for i in range(1,total_element+1):
    loc_to_glob_list.append(global_no(i,4))

#loc_to_glob_list[0][2]
#------------------------------ Main structure of the code --------------------
for element_no in range (total_element):
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
                
            
            def shape_function(xi, eta, iloc, jloc):
                shfunc = shape_func[iloc](xi,eta)*shape_func[jloc](xi,eta)
                return shfunc
            
            
            M[global_i-1,global_j-1] = gauss(shape_function(xi, eta, iloc, jloc),quadrature_points, weights)
            print(iloc, global_i, jloc, global_j)
#spy(M)                
            # here we should input the mass matrix function
                
            # here we should input the mass striffness function

#------------------------- global node number ---------------------------------






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