## 2D new method
#points, weights = np.polynomial.legendre.leggauss(3)
#points.shape = [points.shape[0], 1]
import numpy as np
import matplotlib.pyplot as plt
# quadrature_points, weights = np.polynomial.legendre.leggauss(3)
qpoint = np.array(([-np.sqrt(0.6), -np.sqrt(0.6)],[0, -np.sqrt(0.6)],[np.sqrt(0.6), -np.sqrt(0.6)],
                 [-np.sqrt(0.6), 0],              [0, 0],            [np.sqrt(0.6), 0],
                 [-np.sqrt(0.6), np.sqrt(0.6)],   [0, np.sqrt(0.6)], [np.sqrt(0.6), np.sqrt(0.6)]))
qweight = np.array(([25/81],[40/81],[25/81],
                   [40/81],[64/81],[40/81],
                   [25/81],[40/81],[25/81]))
C = 0.05
c = 0.1
dx = 0.125
dy = 0.1667
dt = C*dx*dy/(c*(dy+dx))
N_e_r = 4
N_e_c= 3
local_node_no = 4
total_element=12
total_sloc = 4
total_nodes = total_element * local_node_no
vol_qp = 9 # degree of polynomial**2
M = np.zeros((total_nodes, total_nodes))
K = np.zeros((total_nodes, total_nodes))
F = np.zeros((total_nodes, total_nodes))

# global node coordinates
def coordinates(e):            
    col =int(np.ceil((e+1)/N_e_r))
    row = int((e+1) - (col -1) * N_e_r)
    co_ordinates = np.array(([dx*(row-1), dy*(col-1)],
                             [dx*row  , dy*(col-1)],
                             [dx*(row-1), dy*(col)],
                             [dx*row  , dy*(col)]))
    return co_ordinates

#-------------------- shape func, ddx and ddy ---------------------------------
shape_func = {0:lambda xi,eta: 1/4*(1-xi)*(1-eta),
              1:lambda xi,eta: 1/4*(1+xi)*(1-eta),
              2:lambda xi,eta: 1/4*(1-xi)*(1+eta),
              3:lambda xi,eta: 1/4*(1+xi)*(1+eta)}

ddx_shape_func = {0:lambda eta: -1/4*(1-eta),
                  1:lambda eta:  1/4*(1-eta),
                  2:lambda eta: -1/4*(1+eta),
                  3:lambda eta:  1/4*(1+eta)}

ddy_shape_func = {0:lambda xi: -1/4*(1-xi),
                  1:lambda xi: -1/4*(1+xi),
                  2:lambda xi:  1/4*(1-xi),
                  3:lambda xi:  1/4*(1+xi)}

dx_dxi = lambda eta: 1/4*((eta-1)*coordinates(e)[0,0] +(1-eta)*coordinates(e)[1,0] -(1+eta)*coordinates(e)[2,0] +(1+eta)*coordinates(e)[3,0])
dy_dxi = lambda eta: 1/4*((eta-1)*coordinates(e)[0,1] +(1-eta)*coordinates(e)[1,1] -(1+eta)*coordinates(e)[2,1] +(1+eta)*coordinates(e)[3,1])
dx_deta = lambda xi: 1/4*((xi-1)*coordinates(e)[0,0] -(xi+1)*coordinates(e)[1,0] +(1-xi)*coordinates(e)[2,0] +(1+xi)*coordinates(e)[3,0])
dy_deta = lambda xi: 1/4*((xi-1)*coordinates(e)[0,1] -(xi+1)*coordinates(e)[1,1] +(1-xi)*coordinates(e)[2,1] +(1+xi)*coordinates(e)[3,1])

#------------------------------ global node numbering -------------------------
def global_no(e,N_e_r):
    col=int(np.ceil(e/N_e_r))
    glob_no = np.array([(col-1)*2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*(e-1)+2, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+2])
    return glob_no

loc_to_glob_list=[]
for i in range(1,total_element+1):
    loc_to_glob_list.append(global_no(i,N_e_r))

#------------------------------ Main structure of the code --------------------
for e in range (total_element):    # element numbers starts from 1 
    
    def jacobian(xi, eta, e):
        a = 1/4 * np.matrix(([-(1-eta), 1-eta, -(1+eta), 1+eta],
                                   [-(1-xi), -(1+xi), 1-xi, 1+xi])) 
        jacobi = a.dot(coordinates(e))
        return jacobi

    
    def det(jacobian):
        return np.linalg.det(jacobian)
    
    
    def inverse(jacobian):
        return np.linalg.inv(jacobian)
    
    
    for iloc in range(local_node_no):     
        global_i = loc_to_glob_list[e][iloc]
        for jloc in range(local_node_no):   
            global_j = loc_to_glob_list[e][jloc]           

#==============================================================================
            # M matrix Gaussian integration
            # loop over all qp
            M[global_i-1,global_j-1] = 0
            K[global_i-1,global_j-1] = 0
            det_jac = 0
            jac_matrix = np.zeros((2,2))
            for i in range(vol_qp):
                # Mass matrix
                jac_matrix = np.array(([dx_dxi, dy_dxi],[dx_deta, dy_deta]))
                det_jac = dx_dxi(qpoint[i,1])*dy_deta(qpoint[i,0]) - dy_dxi(qpoint[i,1]) * dx_deta(qpoint[i,0])
                inv_jac_matrix = 1/det_jac * np.array(([dy_deta, dy_dxi],[dx_deta, dx_dxi]))
                M[global_i-1,global_j-1] = M[global_i-1,global_j-1] + qweight[i] * shape_func[iloc](qpoint[i,0], qpoint[i,1]) * shape_func[jloc](qpoint[i,0], qpoint[i,1])* det_jac

#==============================================================================
                # Stiffness matrix
                K[global_i-1,global_j-1] =(K[global_i-1,global_j-1] + qweight[i] * c*dt*(shape_func[jloc](qpoint[i,0], qpoint[i,1]) * ddx_shape_func[iloc](qpoint[i,1]) * inverse(jacobian(qpoint[i,0], qpoint[i,1],e))[0,0]+
                                                                                         shape_func[jloc](qpoint[i,0], qpoint[i,1]) * ddy_shape_func[iloc](qpoint[i,0]) * inverse(jacobian(qpoint[i,0], qpoint[i,1],e))[1,1])*
                                                                                         det_jac)
            
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
L_xi = np.array([-0.5, 0.5,  1  ,  1  , 0.5, -0.5, -1  , -1  ])
L_eta = np.array([-1  ,-1  , -0.5,  0.5, 1  ,  1  ,  0.5, -0.5])
L_weights=np.array([1,1,1,1,1,1,1,1])
gp = len(L_xi)
sloc = 4
sdot = np.zeros(sloc)
n_hat = np.zeros(sloc)
domain_norm = [0,0,1]

for e in range(total_element):      # element numbering starts from 0
    
    # calculating the center of each element used to calculate the sign of normal
    def e_centre(e):
        x_centre =0
        y_centre =0
        for i in range(local_node_no):
            x_centre = x_centre + coordinates(e)[i,0]         
            y_centre = y_centre + coordinates(e)[i,1]

        z_centre = 0
        e_centre = [x_centre/local_node_no, y_centre/local_node_no, z_centre]
        return e_centre
    
    # =========================================================================
    for s in range(total_sloc):                     # number of surfaces = no.elements * no.faces in each e
        for siloc in range(sloc):                   # = 4, number of local surfaces
            
            # =================================================================
            # function giving global node numbers of 2 points of each face
            def s_glob_node(e,siloc):                
                global_nod ={0: [loc_to_glob_list[e][0] , loc_to_glob_list[e][1]],
                             1: [loc_to_glob_list[e][1] , loc_to_glob_list[e][3]],
                             2: [loc_to_glob_list[e][3] , loc_to_glob_list[e][2]],
                             3: [loc_to_glob_list[e][2] , loc_to_glob_list[e][0]]}
                return global_nod[siloc]
            
            # =================================================================
            # Jacobian
#            for gp in range(len(L_xi)):
#                jacamin =0
#                for sl in range(sloc):
#                    jacamin = jacamin + np.sqrt(dx_dxi(L_eta[gp])**2 + dy_dxi(L_eta[gp])**2)
#                print(jacamin)
            
            det_jac = {0: np.sqrt(dx_dxi(-1)**2 + dy_dxi(-1)**2),
                       1: np.sqrt(dx_deta(1)**2 + dy_deta(1)**2),
                       2: np.sqrt(dx_dxi(1)**2 + dy_dxi(1)**2),
                       3: np.sqrt(dx_deta(-1)**2 + dy_deta(-1)**2)}
            
            # =================================================================
            # normal to the boundary lines 
            snormal = {0: np.cross([dx_dxi(-1) , dy_dxi(-1)  ,0], domain_norm),                  # n_ds of the line (-1,-1) and (-1,1)
                       1: np.cross([dx_deta(1) , dy_deta(1)  ,0], domain_norm),                  # n_ds of the line (-1,1)  and (1,1)
                       2: np.cross([dx_dxi(1)  , dy_dxi(1)   ,0], domain_norm),      # n_ds of the line (1,1)   and (1,-1)
                       3: np.cross([dx_deta(-1), dy_deta(-1) ,0], domain_norm)}  # n_ds of the line (-1,1)  and (-1,-1)
            
            # vector from the e_centre to a node on a boundary line
            r = {0: np.subtract([coordinates(e)[0,0], coordinates(e)[0,1],0] , e_centre(e)),
                 1: np.subtract([coordinates(e)[1,0], coordinates(e)[1,1],0] , e_centre(e)),
                 2: np.subtract([coordinates(e)[3,0], coordinates(e)[3,1],0] , e_centre(e)),
                 3: np.subtract([coordinates(e)[2,0], coordinates(e)[2,1],0] , e_centre(e))}
            
            # dot product of Snormal and r 
            for s in range(sloc):
                sdot[s] = np.dot(snormal[s], r[s])
                if sdot[s] <= 0:
                    snormal[s] = snormal[s] * (-1)
            
            # sign of normal to the boundary
            n_hat = {0: np.sign(snormal[0][1]),
                     1: np.sign(snormal[1][0]),
                     2: np.sign(snormal[2][1]),
                     3: np.sign(snormal[3][0])}
                    
            # =================================================================
            sinod = s_glob_node(e,siloc)                # gives two nodes numbrs of each i-surface
            
            for sjloc in range(sloc):                   # = 4, number of local surfaces
                sjnod = s_glob_node(e,sjloc)            # gives two nodes numbrs of each j-surface
          
                # =============================================================
                # boundary Gaussian integration: loop over all qp
                flux = 0
                for gi in range(gp):
                    flux =  flux + L_weights[gi] * dt * c * n_hat[sjloc] * shape_func[siloc](L_xi[gi],L_eta[gi])*shape_func[sjloc](L_xi[gi],L_eta[gi]) *det_jac[sjloc]
                F[sinod[0]-1, sjnod[0]-1] = F[sinod[1]-1, sjnod[1]-1] = flux/2
                
                
                
#                def L_gauss(f, L_xi, L_eta, L_weights):  
#                    flux = 0
#                    for gi in range(gp):
#                        flux =  flux + L_weights[gi] * f(L_xi[gi], L_eta[gi]) * n_hat[sjloc] *det_jac[sjloc]
#                    return flux
#                # =============================================================
#                # values of 2 points of a line boundary are the same
#                F[s_glob_node(e,siloc)[0]-1, s_glob_node(e,sjloc)[0]-1] = F[s_glob_node(e,siloc)[1]-1, s_glob_node(e,sjloc)[1]-1] = (
#                        L_gauss(lambda xi,eta: dt * c * n_hat[sjloc] * shape_func[siloc](xi,eta)*shape_func[sjloc](xi,eta), L_xi, L_eta, L_weights))

# plt.spy(F)



    







            # vector from centre to one node on a boundary line
            # r = {0: np.subtract([coordinates(e)[0,0], coordinates(e)[0,1],0] , e_centre(e)),
            #      1: np.subtract([coordinates(e)[1,0], coordinates(e)[1,1],0] , e_centre(e)),
            #      2: np.subtract([coordinates(e)[3,0], coordinates(e)[3,1],0] , e_centre(e)),
            #      3: np.subtract([coordinates(e)[2,0], coordinates(e)[2,1],0] , e_centre(e))}
            
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
########################## line gaussian integration ##########################
#                def L_gauss(f, L_xi, L_eta, L_weights):  
#                    xi=eta = np.zeros(len(L_xi))
#                    for i in range(len(L_xi)):
#                        xi[i] =L_xi[i]
#                        eta[i] =L_eta[i]
#                    answer = 0
#                    for i in range(len(L_xi)):
#                        answer =  answer + L_weights[i] * f(xi[i], eta[i]) * n_hat[sjloc] *det_jac[siloc]
#                    return answer

############################## M vol Gaussian int ###############################
#            def gauss(f, quadrature_points, weights):
#                xi=eta = np.zeros(len(quadrature_points))
#                for i in range(len(quadrature_points)):
#                    xi[i] =quadrature_points[i]
#                    eta[i] =quadrature_points[i]
#                answer = 0
#                    
#                for i in range(len(quadrature_points)):
#                    for j in range(len(quadrature_points)):     
#                        answer =  answer + weights[i] * weights[j] * f(xi[i], eta[j]) * det(jacobian(xi[i], eta[j], element_no))
#                return answer

#loc_to_glob_list={}
#for i in range(1,total_element+1):
#    loc_to_glob_list[i]= (global_no(i,4))



#
#xi =  [-sqrt(0.6),          0,  sqrt(0.6), -sqrt(0.6),   0, sqrt(0.6), -sqrt(0.6),         0, sqrt(0.6)]
#eta = [-sqrt(0.6), -sqrt(0.6), -sqrt(0.6),          0,   0,         0,  sqrt(0.6), sqrt(0.6), sqrt(0.6)]
#w_xi = [      5/9,        8/9,        5/9,        5/9, 8/9,       5/9,        5/9,       8/9,       5/9]
#w_eta= [      5/9,        5/9,        5/9,        8/9, 8/9,       8/9,        5/9,       5/9,       5/9]