## 2D new method
#points, weights = np.polynomial.legendre.leggauss(3)
#points.shape = [points.shape[0], 1]
quadrature_points = [-sqrt(0.6),0, sqrt(0.6)]
weights = [5/9, 8/9, 5/9]

def gauss(f, quadrature_points, weights):
    y=x = np.zeros(len(quadrature_points))
    for i in range(len(quadrature_points)):
        x[i] =quadrature_points[i]
        y[i] =quadrature_points[i]
    answer = 0
    for i in range(len(quadrature_points)):
        for j in range(len(quadrature_points)):
            answer =  answer+weights[i] * weights[j] * f(x[i], y[j])
    return answer

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
def f(xi,eta):
    fn = shape_func[1](xi,eta)*shape_func[2](xi,eta)
    return fn

print(f(3,2))

#------------------------ cal J, det and coordinates --------------------------
class Element:
    def jacobian(self, xi, eta, element_no):
        self.xi = xi
        self.eta = eta
        self.coordinates = coordinates
        a = 1/4 * np.matrix(([-(1-eta), 1-eta, -(1+eta), 1+eta],
                               [-(1-xi), -(1+xi), 1-xi, 1+xi])) 
        jacobi = a.dot(coordinates(element_no))
        return jacobi
    
    def det(self, jacobian):
        self.jacobian = jacobian
        return np.linalg.det(jacobian)
    
    
    def coordinates(self, element_no):
        self.element_no = element_no
        col =int(ceil(element_no/N_e_r))
        row = int(element_no - (col -1) * N_e_r)
        coordinates = np.matrix(([dx*(row-1), dy*(col-1)],
                                [dx*row  , dy*(col-1)],
                                [dx*(row-1), dy*(col)],
                                [dx*row  , dy*(col)]))
        return coordinates


element = Element()
a= element.coordinates(12)
b=element.jacobian(-1,1,12)
element.det(element.jacobian(-1,1,12))
dx = 0.125
dy = 0.1667
print(coordinates(12))
#------------------------------ Main structure of the code --------------------
xi =  [-sqrt(0.6),          0,  sqrt(0.6), -sqrt(0.6),   0, sqrt(0.6), -sqrt(0.6),         0, sqrt(0.6)]
eta = [-sqrt(0.6), -sqrt(0.6), -sqrt(0.6),          0,   0,         0,  sqrt(0.6), sqrt(0.6), sqrt(0.6)]
w_xi = [      5/9,        8/9,        5/9,        5/9, 8/9,       5/9,        5/9,       8/9,       5/9]
w_eta= [      5/9,        5/9,        5/9,        8/9, 8/9,       8/9,        5/9,       5/9,       5/9]

for e in range (no_elements): # 12
    for i in range(local_node_in_row):      # 0 and 1
        global_i = 1 to 4             #?
        for j in range(local_node_in_column):   # 0 and 1
            global_j = 1 to 4         # ?
            for g in range(quadrature_points):  # 9
                int = int + w_xi[i] * w_eta[j] * f(xi[g],eta[g])
            
#------------------------- global node number cal -----------------------------
loc_2_global_list = {1:[1,2,9,10],     2:[3,4,11,12],    3:[5,6,13,14],    4:[7,8,15,16],
                     5:[17,18,25,26],  6:[19,20,27,28],  7:[21,22,29,30],  8:[23,24,31,32],
                     9:[33,34,41,42], 10:[35,36,43,44], 11:[37,38,45,46], 12:[39,40,47,48]}
def global_no(e,N_e_r):
    col=int(ceil(e/N_e_r))
    glob_no = np.array([(col-1)*2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*(e-1)+2, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+1, (col-1)*2*N_e_r+2*N_e_r+2*(e-1)+2])
    return glob_no


global_no(7,4)
total_e=12
loc_to_glob_list={}
for i in range(1,total_e+1):
    loc_to_glob_list[i]= (global_no(i,4))

  
