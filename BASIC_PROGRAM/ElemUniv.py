import numpy as np

def integration_scheme(npc_1d):
    if npc_1d == 2:
        x = [-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
        w = [1.0, 1.0]
    elif npc_1d == 3:
        x = [-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)]
        w = [5/9, 8/9, 5/9]
    elif npc_1d == 4:
        pass # dolozyc schemat
    
    return x, w

class Surface:
    def __init__(self, npc, position):
        self.N = np.zeros((npc, 4))
        self.npc_1d = npc
        self.position = position
        x_, w_ = integration_scheme(npc)
        self.x = x_.copy()
        self.w = w_.copy()
    
    def calculate_surface(self):
        for i, xi in enumerate(self.x):
            if self.position == 0:
                ksi = xi
                eta = -1
            elif self.position == 1:
                ksi = 1
                eta = xi
            elif self.position == 2:
                ksi = xi
                eta = 1
            elif self.position == 3:
                ksi = -1
                eta = xi

            self.N[i, 0] = 0.25 * (1 - ksi) * (1 - eta)
            self.N[i, 1] = 0.25 * (1 + ksi) * (1 - eta)
            self.N[i, 2] = 0.25 * (1 + ksi) * (1 + eta)
            self.N[i, 3] = 0.25 * (1 - ksi) * (1 + eta)
    
    def __repr__(self):
        return f"dN: {self.N}\n w = {self.w}\n"

# tabelka ksi(smieszne E) oraz eta(smieszne n)
class ElemUniv:
    def __init__(self, npc):
        self.npc_1d = npc
        self.npc = npc * npc 
        self.dN_dksi = np.zeros((self.npc, 4))
        self.dN_deta = np.zeros((self.npc, 4))
        self.N = np.zeros((self.npc, 4))

        self.surfaces = []

        x_, w_ = integration_scheme(self.npc_1d)
        self.x = x_.copy()
        self.w = w_.copy()
        self.w_pc_outer = np.outer(self.w, self.w).flatten()

    def calculate_elem_univ(self):
        pc_count = 0
        # siatka punktow calkowania
        for ksi in self.x:
            for eta in self.x:
                self.dN_dksi[pc_count] = [
                    -0.25 * (1 - ksi),
                    0.25 * (1 - ksi),
                    0.25 * (1 + ksi),
                    -0.25 * (1 + ksi)
                ]
                self.dN_deta[pc_count] = [
                    -0.25 * (1 - eta),
                    -0.25 * (1 + eta),
                    0.25 * (1 + eta),
                    0.25 * (1 - eta)
                ]

                self.N[pc_count, 0] = 0.25 * (1 - ksi) * (1 - eta)
                self.N[pc_count, 1] = 0.25 * (1 + ksi) * (1 - eta)
                self.N[pc_count, 2] = 0.25 * (1 + ksi) * (1 + eta)
                self.N[pc_count, 3] = 0.25 * (1 - ksi) * (1 + eta)

                pc_count += 1
        
        for i in range(4):
            self.surfaces.append(Surface(self.npc_1d, i))
            self.surfaces[i].calculate_surface()
    
    def __repr__(self):
        return f"\n--- ELEMENT UNIWERSALNY ---\nKSI:\n{self.dN_dksi}\nETA:\n{self.dN_deta}"