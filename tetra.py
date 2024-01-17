import numpy as np


class tetra_model():
    #units
    m = 1.0
    mm = 0.001*m
    N = 1.0
    kN = 1000.0*N
    Mpa = N/(mm**2)
    Gpa = 1000*Mpa

    def __init__(self) -> None:
        self.g_nodes = {}
        self.elements = {}
        self.g_dofs = 3

    def node(self,node_tag,*coords):
        node = {}
        node["tag"] = node_tag
        node["coords"] = coords
        self.g_nodes[node_tag] = node
        self.g_nodes[node_tag]["load"] = []
        self.g_nodes[node_tag]["fixity"] = None
        self.g_nodes[node_tag]["disp"] = []

    def fix(self,node_tag,*fixity):
        self.g_nodes[node_tag]["fixity"] = fixity

    def load(self,node_tag,*loads):
        self.g_nodes[node_tag]["load"].append(np.array(loads))
        

    def get_volume(self,nodes):
        coords = [self.g_nodes[i]['coords'] for i in nodes]
        (x1,y1,z1),(x2,y2,z2),(x3,y3,z3),(x4,y4,z4)  = coords
        vol_mat = np.array([
            [1, x1, y1, z1],
            [1, x2, y2, z2],
            [1, x3, y3, z3],
            [1, x4, y4, z4]
        ])

        V = np.linalg.det(vol_mat)/6
        return V

    def M_matrix(self,rho,V):
        M = np.zeros([12,12])
        for i in range(12):
            M[i,i] = 2

        i = 0
        for j in range(3,12):
            M[i,j] = 1
            M[j,i] = 1
            i+=1
        i=0
        for j in range(6,12):
            M[i,j] = 1
            M[j,i] = 1
            i+=1
        i=0
        for j in range(9,12):
            M[i,j] = 1
            M[j,i] = 1
            i+=1
        
        return rho*V*(1/20)*M


    def B_matrix(self,nodes):
        coords = [self.g_nodes[i]['coords'] for i in nodes]
        (x1,y1,z1),(x2,y2,z2),(x3,y3,z3),(x4,y4,z4)  = coords
        vol_mat = np.array([
            [1, x1, y1, z1],
            [1, x2, y2, z2],
            [1, x3, y3, z3],
            [1, x4, y4, z4]
        ])

        V = np.linalg.det(vol_mat)/6

        beta_1 = -1*np.linalg.det([
            [1,y2,z2],
            [1,y3,z3],
            [1,y4,z4]
        ])

        beta_2 = np.linalg.det([
            [1,y1,z1],
            [1,y3,z3],
            [1,y4,z4]
        ])

        
        beta_3 = -1*np.linalg.det([
            [1,y1,z1],
            [1,y2,z2],
            [1,y4,z4]
        ])

        
        beta_4 = np.linalg.det([
            [1,y1,z1],
            [1,y2,z2],
            [1,y3,z3]
        ])

        gamma_1 = np.linalg.det([
            [1,x2,z2],
            [1,x3,z3],
            [1,x4,z4]
        ])

        gamma_2 = -1*np.linalg.det([
            [1,x1,z1],
            [1,x3,z3],
            [1,x4,z4]
        ])

        gamma_3 = np.linalg.det([
            [1,x1,z1],
            [1,x2,z2],
            [1,x4,z4]
        ])

        gamma_4 = -1*np.linalg.det([
            [1,x1,z1],
            [1,x2,z2],
            [1,x3,z3]
        ])

        delta_1 = -1*np.linalg.det([
            [1,x2,y2],
            [1,x3,y3],
            [1,x4,y4]
        ])

        delta_2 = np.linalg.det([
            [1,x1,y1],
            [1,x3,y3],
            [1,x4,y4]
        ])

        delta_3 = -1*np.linalg.det([
            [1,x1,y1],
            [1,x2,y2],
            [1,x4,y4]
        ])

        delta_4 = np.linalg.det([
            [1,x1,y1],
            [1,x2,y2],
            [1,x3,y3]
            
        ])

        B = (1/(6*V))*np.array([
            [beta_1,0,0,beta_2,0,0,beta_3,0,0,beta_4,0,0],
            [0,gamma_1,0,0,gamma_2,0,0,gamma_3,0,0,gamma_4,0],
            [0,0,delta_1,0,0,delta_2,0,0,delta_3,0,0,delta_4],
            [gamma_1,beta_1,0,gamma_2,beta_2,0,gamma_3,beta_3,0,gamma_4,beta_4,0],
            [0,delta_1,gamma_1,0,delta_2,gamma_2,0,delta_3,gamma_3,0,delta_4,gamma_4],
            [delta_1,0,beta_1,delta_2,0,beta_2,delta_3,0,beta_3,delta_4,0,beta_4]
        ])

        return B,V
        
    def D_matrix(self, E,nu):
        coeff = E/((1+nu)*(1-2*nu))
        inside_d = np.array([
                [1-nu,nu,nu,0,0,0],
                [nu,1-nu, nu,0,0,0],
                [nu,nu,1-nu,0,0,0],
                [0,0,0,(1-2*nu)/2,0,0],
                [0,0,0,0,(1-2*nu)/2,0],
                [0,0,0,0,0,(1-2*nu)/2]
        ])
        return coeff*inside_d

        
    def element(self,ele_tag,nodes,E,nu,self_load = 0.0,rho = 0):
        D = self.D_matrix(E,nu)
        B,V = self.B_matrix(nodes)
        total_load = V*self_load
        #for node in nodes:
        #    self.load(node,0,0,-total_load/4)
        
        K = V*(np.dot(np.transpose(B),np.dot(D,B)))
        M =self.M_matrix(rho,V)
        dofs = np.array([[3*i-2,3*i-1,3*i] for i in nodes])
        dofs = dofs.reshape(1,12)[0]

        mask_mat = {}
        for i,j in enumerate(dofs):
            mask_mat[j] = i

        my_element =  {"mat":K,"mask":mask_mat,"mass":M,"nodes":nodes,"Volume":V,"B":B,"D":D,"stress":[],"strain":[]}
        self.elements[ele_tag] = my_element
        

    def merge_matrices(self,stiffness = True):
        dof = len(self.g_nodes.keys())*self.g_dofs
        big = np.zeros([dof,dof])
        for i in range(dof):
            for j in range(dof):
                for _,k in self.elements.items():
                    if ((i+1) in k["mask"]) and ((j+1) in k["mask"]):
                        if stiffness:
                            big[i][j] += k["mat"][k["mask"][i+1],k["mask"][j+1]]
                        else:
                            big[i][j] += k["mass"][k["mask"][i+1],k["mask"][j+1]]
        return big

    def merge_k(self):
        big_k = self.merge_matrices()
        return big_k

    def merge_M(self):
        big_m = self.merge_matrices(False)
        return big_m

    def generate_final_dofs(self,type="disp"):
        """Generates dof vector or load vector. For dof vector 1 is fixed 0 is free

        Args:
            type (str, optional): . Defaults to "disp".

        Returns:
            list: dof vector or load vector
        """
        dofs = []
        all_nodes = sorted(self.g_nodes.keys())
        free_node_tags = []
        free_nodes = 0
        
        for item in all_nodes:    
            fixity = self.g_nodes[item]["fixity"]
            if fixity is not None:
                u1,v1,w1 = fixity
                dofs += [u1,v1,w1]
            else:
                dofs+=[0,0,0]
                free_nodes += 1
                free_node_tags.append(self.g_nodes[item]["tag"])

        if type=="force": #directly returns the reduced F
            max_time_steps = 0
            for item in all_nodes:    
                time_length = len(self.g_nodes[item]["load"])
                if time_length > max_time_steps:
                    max_time_steps  = time_length
            avail_dofs = self.g_dofs * free_nodes
            force_vector = np.zeros([avail_dofs,max_time_steps])
            
            curr_dof = 0
            for item in free_node_tags:
                load = self.g_nodes[item]["load"]
                len_load = len(load)
                for i in range(len_load):
                    cur_load = load[i]
                    force_vector[curr_dof,i] = cur_load[0]
                    force_vector[curr_dof+1,i] = cur_load[1]
                    force_vector[curr_dof+2,i] = cur_load[2]
                curr_dof+=3
            return force_vector
                    

        return dofs


    def reduce_k(self,merged_k,dofs:list):
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
        


    def start_point(self):
        dof = self.generate_final_dofs()
        K = self.merge_k()
        M = self.merge_M()
        #print(K)
        red_k = self.reduce_k(K,dof)
        red_M = self.reduce_k(M,dof)
        red_f = self.generate_final_dofs(type="force")
        zero_locs = [x for x,y in enumerate(dof) if y==0]
        disp = [0 for i in dof]
        return dof, red_k, red_f,zero_locs,disp,red_M
        

    def update_g_nodes_disp(self,disp,f_disp,zero_locs,step=0):
        rows = len(f_disp)
        for i in range(rows):
            disp[zero_locs[i]] = f_disp[i,step]

        for node in self.g_nodes:
            self.g_nodes[node]["disp"].append([disp[0],disp[1],disp[2]])
            disp.pop(0)
            disp.pop(0)
            disp.pop(0)
        
    def update_elements(self,step=0):
        for _,element in self.elements.items():
            disp = []
            for node in element["nodes"]:
                disp += self.g_nodes[node]["disp"][step]
            disp = np.array(disp)
            disp = disp.reshape(len(disp),1)
            element["strain"].append(np.dot(element["B"],disp))
            element["stress"].append(np.dot(element["D"],element["strain"]))


    def analyse(self):
        dof, red_k, red_f, zero_locs,disp,_ = self.start_point()
        f_disp =np.dot(np.linalg.inv(red_k),red_f)
        self.update_g_nodes_disp(disp,f_disp,zero_locs)
        self.update_elements()


    def dyn_analyse(self,dt=0.0005):
        dof, red_k, red_f, zero_locs,disp,red_M = self.start_point()
        dofs,steps = red_f.shape
        n_time_inc = steps
        red_disp = np.zeros([dofs,n_time_inc+1])
        self.update_g_nodes_disp(disp.copy(),red_disp,zero_locs,0)
        self.update_elements()
        vel = np.zeros([dofs,n_time_inc+1])        
        disp_vector = self.get_dyn_disp(red_f,red_k,red_M,red_disp,vel,n_time_inc,dt,disp,zero_locs)
        return disp_vector
        

    def get_dyn_disp(self,F_vector,red_k,red_M,red_disp,vel,n_time_inc,dt,disp,zero_locs):
        accl_0 = np.dot(np.linalg.inv(red_M),F_vector[:,0]-np.dot(red_k,red_disp[:,0]))
        dispNdelta = red_disp[:,0]-dt*vel[:,0]+0.5*dt*dt*accl_0
        Feff = F_vector[:,0]+np.dot((2/(dt*dt))*red_M-red_k,red_disp[:,0])-np.dot((1/(dt*dt))*red_M,dispNdelta)
        red_disp[:,1] = np.dot(np.linalg.inv((1/(dt*dt))*red_M),Feff)
        self.update_g_nodes_disp(disp.copy(),red_disp,zero_locs,1)
        self.update_elements(1)
        
        for i in range(1,n_time_inc):
            Feff = F_vector[:,i]+np.dot((2/(dt*dt))*red_M-red_k,red_disp[:,i])-np.dot((1/(dt*dt))*red_M,red_disp[:,i-1])
            red_disp[:,i+1] = np.dot(np.linalg.inv((1/(dt*dt))*red_M),Feff)
            self.update_g_nodes_disp(disp.copy(),red_disp,zero_locs,i+1)
            self.update_elements(i+1)

        return red_disp
            