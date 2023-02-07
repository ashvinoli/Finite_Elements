import numpy as np

#units
m = 1.0
mm = 0.001*m
N = 1.0
kN = 1000.0*N
Mpa = N/(mm**2)
Gpa = 1000*Mpa


g_nodes = {}
elements = {}
g_dofs = 2

def node(node_tag,*coords):
    node = {}
    node["tag"] = node_tag
    node["coords"] = coords
    g_nodes[node_tag] = node
    g_nodes[node_tag]["load"] = None
    g_nodes[node_tag]["fixity"] = None
    g_nodes[node_tag]["disp"] = [0,0]

def fix(node_tag,*fixity):
    g_nodes[node_tag]["fixity"] = fixity

def load(node_tag,*loads):
    if g_nodes[node_tag]["load"] is None:
        g_nodes[node_tag]["load"]= np.array(loads)
    else: 
        g_nodes[node_tag]["load"]+= np.array(loads)

def B_matrix(nodes,element_type = "triangle"):
    coords = [g_nodes[i]['coords'] for i in nodes]
    (x1,y1),(x2,y2),(x3,y3)  = coords
    a_mat = np.array([[1.0,1.0,1.0],[x1,x2,x3],[y1,y2,y3]])
    mat_det = ((1/2)*np.linalg.det(a_mat))
    A  = abs(mat_det)
    inside_b = np.array([
        [y2-y3,0,y3-y1,0,y1-y2,0],
        [0,x3-x2,0,x1-x3,0,x2-x1],
        [x3-x2,y2-y3,x1-x3,y3-y1,x2-x1,y1-y2]
    ])
    b_mat = (1/(2*mat_det))*inside_b
    return b_mat,A,inside_b
    
def D_matrix(E,mu,stress_condition):
    coeff = None
    inside_d = None
    if stress_condition == "plane_stress":
        coeff = E/(1-mu**2)
        inside_d = np.array([
            [1,mu,0],
            [mu,1,0],
            [0,0,(1-mu)/2]
        ])
    else:
        coeff = E/((1+mu)*(1-2*mu))
        inside_d = np.array([
            [1-mu,mu,0],
            [mu,1-mu,0],
            [0,0,(1-2*mu)/2]
        ])
    return coeff*inside_d, inside_d, coeff

def element(ele_tag,nodes,E,mu,thickness,self_load = 0.0,stress_condition = "plane_stress"):
    D,_,_ = D_matrix(E,mu, stress_condition)
    B,A,_ = B_matrix(nodes)
    t = 1 if stress_condition=="plain_strain" else thickness
    total_load = A*t*self_load
    for node in nodes:
        load(node,0,-total_load/3)
    
    K = t*A*(np.dot(np.transpose(B),np.dot(D,B)))
    dofs = np.array([[2*i-1,2*i] for i in nodes])
    dofs = dofs.reshape(1,6)[0]

    mask_mat = {}
    for i,j in enumerate(dofs):
        mask_mat[j] = i

    my_element =  {"mat":K,"mask":mask_mat,"nodes":nodes,"Area":A,"B":B,"D":D,"t":t,"strain":None,"stress":None}
    elements[ele_tag] = my_element
    

def merge_k():
    dof = len(g_nodes.keys())*g_dofs
    big_k = np.zeros([dof,dof])
    for i in range(dof):
        for j in range(dof):
            for _,k in elements.items():
                if ((i+1) in k["mask"]) and ((j+1) in k["mask"]):
                    big_k[i][j] += k["mat"][k["mask"][i+1],k["mask"][j+1]]
    return big_k

def generate_final_dofs(type="disp"):
    """Generates dof vector or load vector. For dof vector 1 is fixed 0 is free

    Args:
        type (str, optional): . Defaults to "disp".

    Returns:
        list: dof vector or load vector
    """
    dofs = []
    all_nodes = sorted(g_nodes.keys())
    for item in all_nodes:
        if type == "disp":
            fixity = g_nodes[item]["fixity"]
        else:
            fixity = g_nodes[item]["load"]
        if fixity is not None:
            u1,v1 = fixity
            dofs += [u1,v1]
        else:
            dofs+=[0,0]
    return dofs


def reduce_k(merged_k,dofs:list):
    avail_dofs = dofs.count(0)
    reduced_k = np.zeros((avail_dofs,avail_dofs))
    zero_locs = [x for x,y in enumerate(dofs) if y==0]
    mask= {}
    for index,item in enumerate(zero_locs):
        mask[index] = item

    for i in range(avail_dofs):
        for j in range(avail_dofs):
            reduced_k[i,j] = merged_k[mask[i],mask[j]]
    
    return reduced_k
    
def reduce_f(f_vector:list,dofs):
    zero_locs = [x for x,y in enumerate(dofs) if y==0]
    reduced_f = [f_vector[i] for i in zero_locs]
    return np.array(reduced_f)


def analyse():
    dof = generate_final_dofs()
    K = merge_k()
    print(K)
    red_k = reduce_k(K,dof)
    F = generate_final_dofs(type="force")
    red_f = reduce_f(F,dof)
    zero_locs = [x for x,y in enumerate(dof) if y==0]
    disp = [0 for i in dof]
    f_disp =np.dot(np.linalg.inv(red_k),red_f)

    for i,item in enumerate(f_disp):
        disp[zero_locs[i]] = item

    for node in g_nodes:
        g_nodes[node]["disp"] = disp[0],disp[1]
        disp.pop(0)
        disp.pop(0)

    for _,element in elements.items():
        disp = []
        for node in element["nodes"]:
            disp += g_nodes[node]["disp"]
        disp = np.array(disp)

        element["strain"] = np.dot(element["B"],disp)
        element["stress"] = np.dot(element["D"],element["strain"])

def ele_load(element_tag:int,node_i:int,node_j:int,P_i:float,P_j:float):
    """Applies load to element. Convention for node_i and node_j is taken in such a way that while node_i is just below node_j, theta is zero. The face is vertical with load facing to right. And it rotates in anti-clockwise direction. So to apply load in bottom horizontal face, the load facing downwards, the node i will be the one in right and node j to the left. If opposite is taken, load will face up i.e. to face down the element will have to rotate by more than 180 degrees, vertical being zero degrees.

    Args:
        element_tag (int): _description_
        node_i (int): _description_
        node_j (int): _description_
        P_i (float): _description_
        P_j (float): _description_
    """
    x1,y1 = g_nodes[node_i]["coords"]
    x2,y2 = g_nodes[node_j]["coords"]
    req_element = elements[element_tag]
    del_x = x1-x2
    del_y = y2-y1
    if del_y==0:
        sign = -1 if del_x<0 else 1
        theta = np.arctan(sign*np.inf)
    else:
        theta = np.arctan(del_x/del_y)
    theta=abs(theta)
    L = ((x2-x1)**2+(y2-y1)**2)**(1/2)
    
    #signs:
    if del_x >=0 and del_y>=0:
        signx = 1
        signy = 1
    elif del_x>=0 and del_y<0:
        signx = -1
        signy = 1
    elif del_x<0 and del_y<0:
        signx = -1
        signy = -1
    else:
        signx = 1
        signy = -1

    load_xi = signx*(req_element["t"]*L/6)*(2*P_i*np.cos(theta)+P_j*np.cos(theta))
    load_yi = signy*(req_element["t"]*L/6)*(2*P_i*np.sin(theta)+P_j*np.sin(theta)) 
    load(node_i,
        load_xi,
        load_yi
    )
    
    load_xj =signx* (req_element["t"]*L/6)*(P_i*np.cos(theta)+2*P_j*np.cos(theta))
    load_yj = signy*(req_element["t"]*L/6)*(P_i*np.sin(theta)+2*P_j*np.sin(theta))
    load(node_j,
        load_xj,
        load_yj
    )
    