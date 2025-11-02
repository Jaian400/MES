import numpy as np
from integrate import V_ff
from visualize_grid import visualize_grid
from scipy.integrate import dblquad

# tabelka ksi(smieszne E) oraz eta(smieszne n)
class ElemUniv:
    def __init__(self, npc):
        self.npc_1d = npc
        self.npc = npc * npc 
        self.dN_dksi = np.zeros((self.npc, 4))
        self.dN_deta = np.zeros((self.npc, 4))

        self.w_pc_outer = []

    def calculate_elem_univ(self):
        if self.npc_1d == 2:
            x = [-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
            w = [1.0, 1.0]
            self.w_pc_outer = np.outer(w, w).flatten()
        elif self.npc_1d == 3:
            x = [-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)]
            w = [5/9, 8/9, 5/9]
            self.w_pc_outer = np.outer(w, w).flatten()
        
        # print(f"self.w_pc_outer {self.w_pc_outer}")
        # kazdy obrot to kolejny punkt calkowania; count bedzie 4 dla siatki 2 oraz 9 dla siatki 3
        pc_count = 0

        # siatka punktow calkowania
        for ksi in x:
            for eta in x:
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
                pc_count += 1
    
    def __repr__(self):
        return f"\n--- ELEMENT UNIWERSALNY ---\nKSI:\n{self.dN_dksi}\nETA:\n{self.dN_deta}"

#jedna kropka z położeniem
class Node:
    def __init__(self, node_id, x, y):
        self.id = int(node_id)
        self.x = float(x)
        self.y = float(y)

    def __repr__(self):
        return f"Node(ID: {self.id}, X: {self.x}, Y: {self.y})"

class Jacobian:
    def __init__(self):
        self.J = np.zeros((2,2))
        self.J1 = np.zeros((2,2)) # macierz odwrotna
        self.detJ = 0.0

    def calculate(self, dN_dksi_row, dN_deta_row, element_nodes_coords):
        # dot -> iloczyn skalarny; sum(element wiersza pc razy koordy po x(0) i potem y(1)) 
        self.J[0, 0] = np.dot(dN_dksi_row, element_nodes_coords[:, 0]) # lewy gorny - x po ksi
        self.J[0, 1] = np.dot(dN_dksi_row, element_nodes_coords[:, 1]) # prawy gorny - y po ksi

        self.J[1, 0] = np.dot(dN_deta_row, element_nodes_coords[:, 0]) # lewy dolny - x po eta
        self.J[1, 1] = np.dot(dN_deta_row, element_nodes_coords[:, 1]) # prawy dolny - y po eta

        self.detJ = np.linalg.det(self.J)
        self.J1 = np.linalg.inv(self.J)
    
    def __repr__(self):
        return f"\n--- JACOBIAN ---\n{self.J}\nDET: {self.detJ}\n\n--- JACOBIAN ODWROTNY ---\n{self.J1}\n"

# platy
class Element:
    def __init__(self, element_id, node_ids):
        self.id = int(element_id)
        self.node_ids = [int(nid) for nid in node_ids]
        self.jacobians = [] # jeden jakobian dla kazdego punktu calkowania
        # jacobian nie zalezy od wezla tylko od zestawu wezlow
        # macierz H - przewodnosci cieplnej
        self.H = np.zeros((4,4))

    def calculate_jacobians(self, nodes_all, elem_univ:ElemUniv):
        # uwaga id zaczyna sie od 1; tu mamy kolejne polozenia wezlow
        # mozna przekazac polozenia albo cale obiekty, ale latwiej chyba macierz coords
        element_nodes_coords = np.array([
            [nodes_all[node_id - 1].x, nodes_all[node_id - 1].y] 
            for node_id in self.node_ids
        ])

        for pc in range(elem_univ.npc):
            dN_dksi_row = elem_univ.dN_dksi[pc] # jeden wiersz to input do iloczynu w jaco, dane pc to dany wiersz
            dN_deta_row = elem_univ.dN_deta[pc]
            
            jacobian = Jacobian()
            jacobian.calculate(dN_dksi_row, dN_deta_row, element_nodes_coords)
            self.jacobians.append(jacobian)

    def calculate_H(self, k, elem_univ):
        # j_count bedzie tyle co npc
        for j_count,j in enumerate(self.jacobians):
            dN_dksi = elem_univ.dN_dksi[j_count, :]
            dN_deta = elem_univ.dN_deta[j_count, :] 

            # uzywamy J1 odwrotnej -> te dN_dcos to bedzie jeden rzad
            dN_dx = j.J1[0, 0] * dN_dksi + j.J1[0, 1] * dN_deta
            dN_dy = j.J1[1, 0] * dN_dksi + j.J1[1, 1] * dN_deta
            # print(f"wiersz {j_count}: {dN_dx, dN_dy}")

            mala_h = elem_univ.w_pc_outer[j_count] * k * (np.outer(dN_dx, dN_dx) + np.outer(dN_dy, dN_dy)) * j.detJ
            # print(f"Mala macierz {j_count}: {mala_h}\n")
            self.H += mala_h
        
    def __repr__(self):
        return f"\nElement(ID: {self.id}, Node IDs: {self.node_ids})\n" + f"{self.jacobians}\n\nH:\n{self.H}\n"

# pelna siatka - wszystko wszedzie naraz
class Grid:
    def __init__(self, nN, nE):
        self.nN = nN  
        self.nE = nE 
        self.nodes:Node = [] 
        self.elements:Element = [] 
        self.bc_nodes = [] # nunerki warunkow brzegowych

    def calculate_elements(self, elem_univ, cond):
        for e in self.elements:
            e.calculate_jacobians(self.nodes, elem_univ)
            e.calculate_H(cond, elem_univ)

# dane 
class GlobalData:
    def __init__(self):
        self.SimulationTime = 0
        self.SimulationStepTime = 0
        self.Conductivity = 0.0
        self.Alfa = 0.0
        self.Tot = 0.0
        self.InitialTemp = 0.0
        self.Density = 0.0
        self.SpecificHeat = 0.0

        self.npc = 2 # zakladamy ze odnosi sie do 1d

def load_data(file_path):
    try:     
        with open(file_path, 'r') as file:
            lines = file.readlines()

            gb = GlobalData()
            gb.SimulationTime = int(lines[0].split()[1])
            gb.SimulationStepTime = int(lines[1].split()[1])
            gb.Conductivity = float(lines[2].split()[1])
            gb.Alfa = float(lines[3].split()[1])
            gb.Tot = float(lines[4].split()[1])
            gb.InitialTemp = float(lines[5].split()[1])
            gb.Density = float(lines[6].split()[1])
            gb.SpecificHeat = float(lines[7].split()[1])

            nN = int(lines[8].split()[2]) 
            nE = int(lines[9].split()[2]) 
            
            grid = Grid(nN, nE)
            
            # Flagi jakby nie bylo po kolei
            reading_nodes = False
            reading_elements = False
            reading_bc = False

            for line in lines:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('*Node'):
                    reading_nodes = True
                    reading_elements = False
                    reading_bc = False
                    continue
                elif line.startswith('*Element'):
                    reading_nodes = False
                    reading_elements = True
                    reading_bc = False
                    continue
                elif line.startswith('*BC'):
                    reading_nodes = False
                    reading_elements = False
                    reading_bc = True
                    continue

                # Wczytywanie danych w zależności od aktywnej sekcji
                if reading_nodes:
                    parts = line.replace(',', ' ').split()
                    node_id, x, y = int(parts[0]), float(parts[1]), float(parts[2])
                    grid.nodes.append(Node(node_id, x, y)) 

                elif reading_elements:
                    parts = line.replace(',', ' ').split()
                    elem_id = int(parts[0])
                    
                    node_ids_for_element = [int(p) for p in parts[1:]]
                    grid.elements.append(Element(elem_id, node_ids_for_element))

                elif reading_bc:
                    # border condition
                    parts = line.replace(',', ' ').split()
                    grid.bc_nodes = [int(p) for p in parts]
                    reading_bc = False


            return gb, grid

    except FileNotFoundError:
        print(f"Błąd: Plik '{file_path}' nie został znaleziony.")
        return None, None
    except Exception as e:
        print(f"Wystąpił nieoczekiwany błąd: {e}")
        return None, None


if __name__ == "__main__":

    file_path = "Test1_4_4.txt"  
    # file_path = "Test2_4_4_MixGrid.txt"
    # file_path = "Test3_31_31_kwadrat.txt"

    global_data, grid_data = load_data(file_path)

    try:
        npc = int(input("Schemat calkowania (2 lub 3, domyślnie 2): "))
        if npc not in [2, 3]: npc = 2
    except ValueError:
        npc = 2

    global_data.npc = npc

    if global_data and grid_data:
        print("\n--- WSPÓŁRZĘDNE WSZYSTKICH WĘZŁÓW ---")
        for node_obj in grid_data.nodes:
            print(f"Węzeł ID: {node_obj.id}, Współrzędne (X, Y): ({node_obj.x}, {node_obj.y})")

        print("\n\n--- ID WĘZŁÓW DLA POSZCZEGÓLNYCH ELEMENTÓW ---")
        for elem_obj in grid_data.elements:
            print(f"Element ID: {elem_obj.id}, ID węzłów: {elem_obj.node_ids}")

        print("\n\n--- BOLDING CONDITION ---")
        for abc in grid_data.bc_nodes:
            print(f"{abc}")
    else:
        print("Nie otrzymano poprawnych danych z load_data\n")
        exit(-1)
    
    elem_univ = ElemUniv(global_data.npc)
    elem_univ.calculate_elem_univ()

    print(elem_univ)

    grid_data.calculate_elements(elem_univ, global_data.Conductivity)
    print(grid_data.elements)

    # visualize_grid(grid_data)
    exit(0)

    print("\n\nCALKOWANIE")
    calka = V_ff()

    f1 = lambda x:5*np.pow(x, 2) + 3*x + 6
    f2 = lambda y, x: 5 * np.power(x,2) * np.power(y,2) + 3*x*y + 6

    calka.quadrature(f2, -1, 1)
    print(calka.output)

    dblquad_, error = dblquad(f2, -100, 100, lambda x: -1, lambda x: 2)
    print(f"dblquad: {dblquad_}")
    print(f"roznica miedzy V_ff a dblquad: {100 * abs(calka.output - dblquad_)/calka.output} %")