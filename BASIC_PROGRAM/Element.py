import numpy as np
from ElemUniv import ElemUniv

class Jacobian:
    def __init__(self):
        self.J = np.zeros((2,2))
        self.J1 = np.zeros((2,2)) # macierz odwrotna
        self.detJ = 0.0

    def calculate(self, dN_dksi_row, dN_deta_row, element_nodes_coords):
        self.J[0, 0] = np.dot(dN_dksi_row, element_nodes_coords[:, 0])
        self.J[0, 1] = np.dot(dN_dksi_row, element_nodes_coords[:, 1]) 

        self.J[1, 0] = np.dot(dN_deta_row, element_nodes_coords[:, 0])
        self.J[1, 1] = np.dot(dN_deta_row, element_nodes_coords[:, 1])

        self.detJ = np.linalg.det(self.J)
        self.J1 = np.linalg.inv(self.J)
    
    def __repr__(self):
        return f"\n--- JACOBIAN ---\n{self.J}\nDET: {self.detJ}\n\n--- JACOBIAN ODWROTNY ---\n{self.J1}\n"

# platy
class Element:
    def __init__(self, element_id, node_ids):
        self.id = int(element_id)
        self.node_ids = [int(nid) for nid in node_ids]
        self.jacobians = [] 
        self.H = np.zeros((4,4))
        self.Hbc = np.zeros((4,4))
        self.P = np.zeros((4))
        self.C = np.zeros((4,4))

    def calculate_jacobians(self, nodes_all, elem_univ:ElemUniv):
        self.jacobians.clear()
        # uwaga id zaczyna sie od 1
        element_nodes_coords = np.array([
            [nodes_all[node_id - 1].x, nodes_all[node_id - 1].y] 
            for node_id in self.node_ids
        ])

        for pc in range(elem_univ.npc):
            dN_dksi_row = elem_univ.dN_dksi[pc] 
            dN_deta_row = elem_univ.dN_deta[pc]
            
            jacobian = Jacobian()
            jacobian.calculate(dN_dksi_row, dN_deta_row, element_nodes_coords)
            self.jacobians.append(jacobian)

    def calculate_H(self, k, elem_univ, dens, spec_heat):
        self.H = np.zeros((4, 4))
        self.C = np.zeros((4, 4))

        for j_count,j in enumerate(self.jacobians):
            dN_dksi = elem_univ.dN_dksi[j_count, :]
            dN_deta = elem_univ.dN_deta[j_count, :] 

            dN_dx = j.J1[0, 0] * dN_dksi + j.J1[0, 1] * dN_deta
            dN_dy = j.J1[1, 0] * dN_dksi + j.J1[1, 1] * dN_deta
            # print(f"wiersz {j_count}: {dN_dx, dN_dy}")

            mala_h = elem_univ.w_pc_outer[j_count] * k * (np.outer(dN_dx, dN_dx) + np.outer(dN_dy, dN_dy)) * j.detJ
            # print(f"Mala macierz {j_count}: {mala_h}\n")
            self.H += mala_h

            self.C += elem_univ.w_pc_outer[j_count] * dens * spec_heat * (np.outer(elem_univ.N[j_count, :], elem_univ.N[j_count, :])) * j.detJ
    
    def calculate_Hbc(self, alfa, elem_univ, nodes_all, t_ot):
        self.Hbc = np.zeros((4,4))
        self.P = np.zeros((4))

        walls = [(0, 1), (1, 2), (2, 3), (3, 0)]
        # a zawsze beda 2 na sciane bo to czworokat
        for i,wall in enumerate(walls):
            node1 = nodes_all[self.node_ids[wall[0]] - 1]
            node2 = nodes_all[self.node_ids[wall[1]] - 1]

            if not node1.BC or not node2.BC:
                continue
            
            detJ = 0.5 * np.sqrt(np.pow((node1.x - node2.x), 2) + np.pow((node1.y - node2.y), 2))
            surface = elem_univ.surfaces[i]
            # print(surface)

            for j, w in enumerate(surface.w): # dla kazdego punktu calkowania
                self.Hbc += alfa * w * np.outer(surface.N[j], surface.N[j]) * detJ
                # if node1.x == 0 and node2.x == 0:
                #     self.P += alfa * w * surface.N[j] * 200 * detJ
                #     continue
                self.P += alfa * w * surface.N[j] * t_ot * detJ
 
        self.H += self.Hbc
        
    def __repr__(self):
        return f"\nElement(ID: {self.id}, Node IDs: {self.node_ids})\n" + f"{self.jacobians}\n\nH:\n{self.H}\n\nHbc:\n{self.Hbc}\n\n{self.P}\n"
