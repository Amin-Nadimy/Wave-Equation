import numpy as np
import matplotlib.pyplot as plt
qpoint = np.array(([-np.sqrt(0.6), -np.sqrt(0.6)],[0, -np.sqrt(0.6)],[np.sqrt(0.6), -np.sqrt(0.6)],
                 [-np.sqrt(0.6), 0],              [0, 0],            [np.sqrt(0.6), 0],
                 [-np.sqrt(0.6), np.sqrt(0.6)],   [0, np.sqrt(0.6)], [np.sqrt(0.6), np.sqrt(0.6)]))
qweight = np.array(([25/81],[40/81],[25/81],
                   [40/81],[64/81],[40/81],
                   [25/81],[40/81],[25/81]))
C = 0.05
c = 0.1
L = 0.5
N_e_r = 4
N_e_c= 3
nt = 10
dx = L/(N_e_r)
dy = L/(N_e_c)
dt = C*dx*dy/(c*(dy+dx))
local_node_no = 4
total_element=N_e_r * N_e_c

total_nodes = total_element * local_node_no
vol_qp = 9 # degree of polynomial**2
M = np.zeros((total_nodes, total_nodes))
K = np.zeros((total_nodes, total_nodes))
U = np.zeros(total_element*local_node_no)                          # Wave matrix
Un = np.zeros(total_element*local_node_no)                         # Dummy variable to save current components of U
U1=U2 = U                                # Dummy matrices to plot 3 time steps
#U_plot = np.ones((3,N*Np))                 # A matrix to save 3 time steps used for plotting the results

U[int(total_element*local_node_no*.3):int(total_element*local_node_no*.8)]=1              # Defining wave components

############### DG nodes and meshgrid #########################################
## DG x-coordinates
x = np.linspace(0, L, N_e_r+1)
x = [ele for ele in x for i in range(2)]
x=x[1:]
x=x[:-1]

## DG y-coordinates
y = np.linspace(0, L, N_e_c+1)
y = [ele for ele in y for i in range(2)]
y=y[1:]
y=y[:-1]

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

ddxi_shape_func = {0:lambda eta: -1/4*(1-eta),
                   1:lambda eta:  1/4*(1-eta),
                   2:lambda eta: -1/4*(1+eta),
                   3:lambda eta:  1/4*(1+eta)}

ddeta_shape_func = {0:lambda xi: -1/4*(1-xi),
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
            for i in range(vol_qp):
                det_jac = dx_dxi(qpoint[i,1])*dy_deta(qpoint[i,0]) - dy_dxi(qpoint[i,1]) * dx_deta(qpoint[i,0])
                
                # Mass matrix
                M[global_i-1,global_j-1] = M[global_i-1,global_j-1] + qweight[i] * shape_func[iloc](qpoint[i,0], qpoint[i,1]) * shape_func[jloc](qpoint[i,0], qpoint[i,1])* det_jac

#==============================================================================
                # Stiffness matrix
                ddx_shape_func = 1/det_jac * (dy_deta(qpoint[i,0]) * ddxi_shape_func[iloc](qpoint[i,1]) - 
                                    dy_dxi(qpoint[i,1]) * ddeta_shape_func[iloc](qpoint[i,0]))
                ddy_shape_func = 1/det_jac * (dx_dxi(qpoint[i,1]) * ddeta_shape_func[iloc](qpoint[i,0]) - 
                                    dx_deta(qpoint[i,0]) * ddxi_shape_func[iloc](qpoint[i,1]))
                
                K[global_i-1,global_j-1] =(K[global_i-1,global_j-1] + qweight[i] * c*dt*(shape_func[jloc](qpoint[i,0], qpoint[i,1]) * ddx_shape_func +
                                                                                         shape_func[jloc](qpoint[i,0], qpoint[i,1]) * ddy_shape_func)
                                                                                         * det_jac)

###############################################################################
###############################################################################
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
L_xi = [-1/np.sqrt(3), 1/np.sqrt(3),  1   ,    1,   1/np.sqrt(3), -1/np.sqrt(3), -1          , -1  ]
L_eta = [-1          , -1          ,  -1/np.sqrt(3),  1/np.sqrt(3),   1           ,  1           , 1/np.sqrt(3), -1/np.sqrt(3)]
L_qpoint = ([-1/np.sqrt(3), -1],[1/np.sqrt(3), -1],[1, -1/np.sqrt(3)],[1, 1/np.sqrt(3)],
                      [3**(-0.5), 1],[-3**(-0.5), 1],[-1, 3**(-0.5)],[-1, -3**(-0.5)])

s_ng = {0: [[-3**(-0.5), -1], [3**(-0.5), -1]],
        1: [[1, -3**(-0.5)],  [1, 3**(-0.5)]],
        3: [[3**(-0.5), 1],   [-3**(-0.5), 1]],
        2: [[-1, 3**(-0.5)],  [-1, -3**(-0.5)]]}

L_weights=[1,1]
gp = len(L_qpoint)
#sloc = 2
nsuf = 4
ng =2
total_nsuf = total_element * nsuf
sdot = np.zeros(nsuf)
n_hat = np.zeros(nsuf)
domain_norm = [0,0,1]
F = np.zeros((total_nodes, total_nodes))
#F = np.zeros((total_nodes+2*N_e_r, total_nodes+2*N_e_c))

# =============================================================================
# function giving global node numbers of 2 points of each face
def s_glob_node(e,suf):
    global_nod ={0: [loc_to_glob_list[e][0]-1, loc_to_glob_list[e][1]-1],
                 1: [loc_to_glob_list[e][1]-1, loc_to_glob_list[e][3]-1],
                 3: [loc_to_glob_list[e][3]-1, loc_to_glob_list[e][2]-1],
                 2: [loc_to_glob_list[e][2]-1, loc_to_glob_list[e][0]-1]}
    return global_nod[suf]


#==============================================================================
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
    
    
#------------------------------- Main stucture of the code --------------------    
for e in range(total_element):              # element numbering starts from 0
    for suf in range(nsuf):             # no.faces in each e
        for siloc in range(4):           # = 4, number of local surfaces
            sinod = s_glob_node(e,siloc)    # gives two nodes numbrs of each i-surface
            # =================================================================
            # normal to the boundary lines 
            snormal = {0: np.cross([dx_dxi(-1) , dy_dxi(-1)  ,0], domain_norm),                  # n_ds of the line (-1,-1) and (-1,1)
                       1: np.cross([dx_deta(1) , dy_deta(1)  ,0], domain_norm),                  # n_ds of the line (-1,1)  and (1,1)
                       2: np.cross([dx_dxi(1)  , dy_dxi(1)   ,0], domain_norm),      # n_ds of the line (1,1)   and (1,-1)
                       3: np.cross([dx_deta(-1), dy_deta(-1) ,0], domain_norm)}  # n_ds of the line (-1,1)  and (-1,-1)
            
            # vector from the centre of the element to a node on a boundary line
            r = {0: np.subtract([coordinates(e)[0,0], coordinates(e)[0,1],0] , e_centre(e)),
                 1: np.subtract([coordinates(e)[1,0], coordinates(e)[1,1],0] , e_centre(e)),
                 2: np.subtract([coordinates(e)[3,0], coordinates(e)[3,1],0] , e_centre(e)),
                 3: np.subtract([coordinates(e)[2,0], coordinates(e)[2,1],0] , e_centre(e))}
            
            # dot product of Snormal and r 
            for sl in range(nsuf):
                sdot[sl] = np.dot(snormal[sl], r[sl])
                if sdot[sl] <= 0:
                    snormal[sl] = snormal[sl] * (-1)
            
            # sign of normal for each surface
            n_hat = {0: np.sign(snormal[0][1]),
                     1: np.sign(snormal[1][0]),
                     3: np.sign(snormal[2][1]),
                     2: np.sign(snormal[3][0])}
                    
            # =================================================================
            for sjloc in range(4):                   # = 4, number of local surfaces
                sjnod = s_glob_node(e,sjloc)         # gives two nodes numbrs of each j-surface
          
                # =============================================================
                # boundary Gaussian integration: loop over all ng
                # cal det_jac
                L_det_jac =0
#                F[sinod[0], s_glob_node(e-N_e_r,3)[1]] = 0 #e11, global node 38 to 30
#                F[sinod[1], s_glob_node(e-N_e_r,3)[0]] = 0 #e11, global node 39 to 31
#                F[sinod[0], sinod[0]] = 0 # e11, global node 39 to 39
#                F[sinod[1], sinod[1]] = 0  # e11, global node 47 to 47
#                F[sinod[1], sinod[1]] = 0 # e11, global node 46 to 46
#                F[sinod[0], sinod[0]] = 0  # e11, global node 47 to 47
#                F[sinod[1], s_glob_node(e-1,1)[0]] = 0 # e11, global node 38 to 37
#                F[sinod[0], s_glob_node(e-1,1)[1]] = 0 # e11, global node 46 to 45
                flux=0
                for g in range(ng):
                    if suf==0 or suf==3:
                        eta=s_ng[sl][g][1]
                        L_det_jac = np.sqrt(dx_dxi(eta)**2 + dy_dxi(eta)**2)
                    elif suf==1 or suf==2:
                        xi = s_ng[sl][g][0]
                        L_det_jac = np.sqrt(dx_deta(xi)**2 + dy_deta(xi)**2)
                    #==========================================================
                    # cal flux
                    flux =  (flux+L_weights[g] * c*dt * n_hat[suf] 
                                    * shape_func[siloc](s_ng[siloc][g][0],s_ng[siloc][g][1])
                                    * shape_func[sjloc](s_ng[sjloc][g][0],s_ng[sjloc][g][1]) 
                                    * L_det_jac)
                    
                if suf ==0 and e>>N_e_r:
                    F[s_glob_node(e,0)[0], s_glob_node(e-N_e_r,3)[1]] = flux # 38 to 30
                    F[s_glob_node(e,0)[0], s_glob_node(e-N_e_r,3)[0]] = flux # 38 to 31
                    F[s_glob_node(e,0)[1], s_glob_node(e-N_e_r,3)[1]] = flux # 39 to 30
                    F[s_glob_node(e,0)[1], s_glob_node(e-N_e_r,3)[0]] = flux # 39 to 31
                elif suf==1:
                    F[s_glob_node(e,1)[0], s_glob_node(e,1)[0]] = flux # 39 to 39
                    F[s_glob_node(e,1)[1], s_glob_node(e,1)[1]] = flux # 47 to 47
                elif suf==3:
                    F[s_glob_node(e,3)[0], s_glob_node(e,3)[0]] = F[s_glob_node(e,1)[1], s_glob_node(e,1)[1]]+flux # 47 to 47
                    F[s_glob_node(e,3)[1], s_glob_node(e,3)[1]] = flux # 46 to 46
                elif suf==2 and e>>0 and e!=0 and e!=4 and e!=8:
                    F[s_glob_node(e,2)[0], s_glob_node(e-1,1)[1]] = flux # 46 to 45
                    F[s_glob_node(e,2)[0], s_glob_node(e-1,1)[0]] = flux # 46 to 37
                    F[s_glob_node(e,2)[1], s_glob_node(e-1,1)[1]] = flux # 38 to 45
                    F[s_glob_node(e,2)[1], s_glob_node(e-1,1)[0]] = flux # 38 to 37
                    
#plt.spy(F)                
########################### solving for U #####################################

RHS_cst = (M + K - F)
plt.spy(RHS_cst)
for n in range(nt):                 # Marching in time
    Un = U.copy()
    RHS = RHS_cst.dot(Un)           # saving U^t to be used at the next timestep calculation
    U=np.linalg.solve(M,RHS)        # solving for U(t+1)
    if n==1:                        # saving U at timestep 1 to plot
        U1=U
    elif n==nt//2:                  # saving U at half timestep to plot
        U2=U
    
     






















           
#            L_det_jac = {0: np.sqrt(dx_dxi(-1)**2 + dy_dxi(-1)**2),
#                         1: np.sqrt(dx_deta(1)**2 + dy_deta(1)**2),
#                         2: np.sqrt(dx_dxi(1)**2 + dy_dxi(1)**2),
#                         3: np.sqrt(dx_deta(-1)**2 + dy_deta(-1)**2)}


                
#                def L_gauss(f, L_xi, L_eta, L_weights):  
#                    flux = 0
#                    for gi in range(gp):
#                        flux =  flux + L_weights[gi] * f(L_xi[gi], L_eta[gi]) * n_hat[sjloc] *L_det_jac[sjloc]
#                    return flux
#                # =============================================================
#                # values of 2 points of a line boundary are the same
#                F[s_glob_node(e,siloc)[0]-1, s_glob_node(e,sjloc)[0]-1] = F[s_glob_node(e,siloc)[1]-1, s_glob_node(e,sjloc)[1]-1] = (
#                        L_gauss(lambda xi,eta: dt * c * n_hat[sjloc] * shape_func[siloc](xi,eta)*shape_func[sjloc](xi,eta), L_xi, L_eta, L_weights))




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
#                        answer =  answer + L_weights[i] * f(xi[i], eta[i]) * n_hat[sjloc] *L_det_jac[siloc]
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


#points, weights = np.polynomial.legendre.leggauss(3)
#points.shape = [points.shape[0], 1]
            
            
#loc_to_glob_list={}
#for i in range(1,total_element+1):
#    loc_to_glob_list[i]= (global_no(i,4))



#
#xi =  [-sqrt(0.6),          0,  sqrt(0.6), -sqrt(0.6),   0, sqrt(0.6), -sqrt(0.6),         0, sqrt(0.6)]
#eta = [-sqrt(0.6), -sqrt(0.6), -sqrt(0.6),          0,   0,         0,  sqrt(0.6), sqrt(0.6), sqrt(0.6)]
#w_xi = [      5/9,        8/9,        5/9,        5/9, 8/9,       5/9,        5/9,       8/9,       5/9]
#w_eta= [      5/9,        5/9,        5/9,        8/9, 8/9,       8/9,        5/9,       5/9,       5/9]