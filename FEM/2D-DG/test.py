import numpy as np
import matplotlib.pyplot as plt
qpoint = np.array(([-np.sqrt(0.6), -np.sqrt(0.6)],[0, -np.sqrt(0.6)],[np.sqrt(0.6), -np.sqrt(0.6)],
                   [-np.sqrt(0.6), 0],              [0, 0],            [np.sqrt(0.6), 0],
                   [-np.sqrt(0.6), np.sqrt(0.6)],   [0, np.sqrt(0.6)], [np.sqrt(0.6), np.sqrt(0.6)]))

qweight = np.array(([25/81],[40/81],[25/81],
                    [40/81],[64/81],[40/81],
                    [25/81],[40/81],[25/81]))

s_ng = {0: [[-3**(-0.5), -1], [3**(-0.5), -1]],
        1: [[1, -3**(-0.5)],  [1, 3**(-0.5)]],
        2: [[-1, 3**(-0.5)],  [-1, -3**(-0.5)]],
        3: [[3**(-0.5), 1],   [-3**(-0.5), 1]]}

L_weights=[1,1]

nsuf = 4
ng =2
gp = nsuf * ng
C = 0.05
c=np.array([0.1,0.1])
#c_x=0.1
#c_y=0.1
L = 0.5
N_e_r = 20
N_e_c= 10
nt = 1
domain_norm = [0,0,1]
dx = L/(N_e_r)
dy = L/(N_e_c)
dt = C/((c[0]/dx)+(c[1]/dy))
local_node_no = 4
total_element=N_e_r * N_e_c
total_nodes = total_element * local_node_no
total_nsuf = total_element * nsuf
vol_qp = 9 # degree of polynomial**2

# DG x-coordinates
x = np.linspace(0, L, N_e_r+1)
x = [ele for ele in x for i in range(2)]
x=x[1:]
x=x[:-1]

## DG y-coordinates
y = np.linspace(0, L, N_e_c+1)
y = [ele for ele in y for i in range(2)]
y=y[1:]
y=y[:-1]

sdot = np.zeros(nsuf)
n_hat = np.zeros(nsuf)
M = np.zeros((total_nodes, total_nodes))
K = np.zeros((total_nodes, total_nodes))
F = np.zeros((total_nodes, total_nodes))
U = np.zeros(total_nodes)                          

# initial con
for i in range(8):
    U[N_e_r*2*(i+4)+N_e_r//4+1:N_e_r*2*(i+4)+N_e_r//1+2]=1

#======================== Boundary elements in y-dir ==========================
    boundary_element_ydir=[]
n=0
while n< total_element:
    boundary_element_ydir.append(n)
    n=n+N_e_r
    
#==============================================================================
# global node coordinates
def coordinates(e):            
    col =int(np.ceil((e+1)/N_e_r))
    row = int((e+1) - (col -1) * N_e_r)
    co_ordinates = np.array(([dx*(row-1), dy*(col-1)],
                             [dx*row  , dy*(col-1)],
                             [dx*(row-1), dy*(col)],
                             [dx*row  , dy*(col)]))
    return co_ordinates


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


#------------------------------ global node numbering -------------------------
def global_no(e,N_e_r):
    col=int(np.ceil(e/N_e_r))
    glob_no = np.array([(col-1)*2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*(e-1)+2, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+2])
    return glob_no


loc_to_glob_list=[]
for i in range(1,total_element+1):
    loc_to_glob_list.append(global_no(i,N_e_r))

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
                
                K[global_i-1,global_j-1] =(K[global_i-1,global_j-1] + qweight[i]*dt*(c[0] * shape_func[jloc](qpoint[i,0], qpoint[i,1]) * ddx_shape_func +
                                                                                     c[1] * shape_func[jloc](qpoint[i,0], qpoint[i,1]) * ddy_shape_func)
                                                                                         * det_jac)
    #print(e)
M_inv = np.linalg.inv(M)
###############################################################################
############################ surface integration ##############################
Un_hat = np.zeros(total_element*nsuf)
for n in range(nt):
    Un=U.copy()
    r_cst = M_inv.dot(K)
    r_cst = r_cst.dot(Un)      # M_inv * K * Un
    sum = 0
    for e in range(total_element):
        for suf in range(nsuf):             # no.faces in each e
            sol = e*nsuf+suf#s_glob_node(e,suf)[0] # e*nsuf+suf
            for siloc in range(4):           # = 4, number of local surfaces
                sinod = s_glob_node(e,siloc)    # gives two nodes numbrs of each i-surface
                # =================================================================
                # normal to the boundary lines 
                snormal = {0: np.cross([dx_dxi(-1) , dy_dxi(-1)  ,0], domain_norm),      # n_ds of the line (-1,-1) and (-1,1)
                           1: np.cross([dx_deta(1) , dy_deta(1)  ,0], domain_norm),      # n_ds of the line (-1,1)  and (1,1)
                           2: np.cross([dx_dxi(1)  , dy_dxi(1)   ,0], domain_norm),      # n_ds of the line (1,1)   and (1,-1)
                           3: np.cross([dx_deta(-1), dy_deta(-1) ,0], domain_norm)}      # n_ds of the line (-1,1)  and (-1,-1)
                
                # vector from the centre of the element to a node on a boundary line
                r = {0: np.subtract([coordinates(e)[0,0], coordinates(e)[0,1],0] , e_centre(e)),
                     1: np.subtract([coordinates(e)[1,0], coordinates(e)[1,1],0] , e_centre(e)),
                     2: np.subtract([coordinates(e)[3,0], coordinates(e)[3,1],0] , e_centre(e)),
                     3: np.subtract([coordinates(e)[2,0], coordinates(e)[2,1],0] , e_centre(e))}
                
                # dot product of Snormal and r 
                for sl in range(nsuf):
                    sdot = np.dot(snormal[sl], r[sl])
                    if sdot <= 0:
                        snormal[sl] = snormal[sl] * (-1)
                
                # sign of normal for each surface
                n_hat = {0: [np.sign(snormal[0][0]), np.sign(snormal[0][1])],
                         1: [np.sign(snormal[1][0]), np.sign(snormal[1][1])],
                         3: [np.sign(snormal[2][0]), np.sign(snormal[2][1])],
                         2: [np.sign(snormal[3][0]), np.sign(snormal[3][1])]}
                        
                # =============================================================
                # cal det_jac
                L_det_jac =0
                flux =0
                #U[sol] = 0
                for g in range(ng):
                    eta=s_ng[siloc][g][1]
                    xi = s_ng[siloc][g][0]
                    L_det_jac = {0: np.sqrt(dx_dxi(eta)**2 + dy_dxi(eta)**2),
                                 1: np.sqrt(dx_deta(xi)**2 + dy_deta(xi)**2),
                                 2: np.sqrt(dx_deta(xi)**2 + dy_deta(xi)**2),
                                 3: np.sqrt(dx_dxi(eta)**2 + dy_dxi(eta)**2)}
                        
                    flux= L_weights[g] * L_det_jac[siloc] * dt *(c.dot(n_hat[siloc]) 
                                                            * shape_func[siloc](s_ng[siloc][g][0],s_ng[siloc][g][1]))
                    
                    if siloc==0:
                        Un_hat = Un[(sol-N_e_r*4)+3]
                    elif siloc==2:
                        Un_hat = Un[sol-5]
                    else:
                        Un_hat = Un[sol]
                        
                    j=0
                    while j <= total_nodes-1:
                        sum = sum + ( Un[sol] 
                                        +  r_cst[sol] 
                                        -  M_inv[sol, j] * flux * Un_hat )
                    
                        j+=1
                    U[sol]=sum

#for n in range(nt):
#    for j in range(total_nodes):
#        for i in range(4):
#            Un=U.copy
#            U[j] = ( Un[j] 
#                    + dt * M_inv[i,j] * K[i,j] * Un[j] 
#                    - dt * M_inv[i,j] * flux * Un[j] )
            
########################### solving for U #####################################
#RHS_cst = (M + K - F)
#for n in range(nt):                 # Marching in time
#    Un = U.copy()
#    RHS = RHS_cst.dot(Un)           # saving U^t to be used at the next timestep calculation
#    U=np.linalg.solve(M,RHS)        # solving for U(t+1)
#    U[0:N_e_r*2]=0#U[N_e_c*2:0]=0
#    if n==1:                        # saving U at timestep 1 to plot
#        U1=U
#    elif n==100:
#        U2=U
#    print(n)
#
########################### Plotting the results ################################
U1_plot=np.zeros(([len(y),len(x)]))
U2_plot=np.zeros(([len(y),len(x)]))

j=0
k=0
while j< len(y):
    i=0
    while i< len(x):
        #U1_plot[j,i]=U1[k]
        U2_plot[j,i]=U[k]
        i+=1
        k+=1
    j+=1
        
#X, Y =np.meshgsum(x,y)           # Creating a mesh gsum
#plt.figure(1)        
#ax = plt.gca(projection='3d')
#ax.plot_surface(X, Y, U2_plot , label='t=final')
##ax.plot_surface(X, Y, U1_plot , label='t=0')
#ax.set_ylabel('$y$')
#ax.set_xlabel('$x$')
#ax.set_zlabel('$U$')
#plt.legend()
#plt.show()