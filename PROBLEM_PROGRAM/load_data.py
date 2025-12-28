from Global import global_data, grid_data, Node
from Element import Element

def load_data(file_path):
    try:     
        with open(file_path, 'r') as file:
            lines = file.readlines()

            global_data.SimulationTime = int(lines[0].split()[1])
            global_data.SimulationStepTime = int(lines[1].split()[1])
            global_data.Conductivity = float(lines[2].split()[1])
            global_data.Alfa = float(lines[3].split()[1])
            global_data.Tot = float(lines[4].split()[1])
            global_data.InitialTemp = float(lines[5].split()[1])
            global_data.Density = float(lines[6].split()[1])
            global_data.SpecificHeat = float(lines[7].split()[1])

            grid_data.nN = int(lines[8].split()[2]) 
            grid_data.nE = int(lines[9].split()[2]) 
            
            # Flagi jakby nie bylo po kolei
            reading_nodes = False
            reading_elements = False
            reading_bc = False
            reading_e_gen = False

            for line in lines:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('*Node'):
                    reading_nodes = True
                    reading_elements = False
                    reading_bc = False
                    reading_e_gen = False
                    continue
                elif line.startswith('*Element'):
                    reading_nodes = False
                    reading_elements = True
                    reading_bc = False
                    reading_e_gen = False
                    continue
                elif line.startswith('*BC'):
                    reading_nodes = False
                    reading_elements = False
                    reading_bc = True
                    reading_e_gen = False
                    continue
                elif line.startswith('*Elset'):
                    reading_nodes = False
                    reading_elements = False
                    reading_bc = False
                    reading_e_gen = True
                    continue

                # Wczytywanie danych w zależności od aktywnej sekcji
                if reading_nodes:
                    parts = line.replace(',', ' ').split()
                    node_id, x, y = int(parts[0]), float(parts[1]), float(parts[2])
                    grid_data.nodes.append(Node(node_id, x, y)) 

                elif reading_elements:
                    parts = line.replace(',', ' ').split()
                    elem_id = int(parts[0])
                    
                    node_ids_for_element = [int(p) for p in parts[1:]]
                    grid_data.elements.append(Element(elem_id, node_ids_for_element))

                elif reading_bc:
                    # border condition
                    parts = line.replace(',', ' ').split()
                    for p in parts:
                        grid_data.bc_nodes.append(int(p))
                
                elif reading_e_gen:
                    parts = line.replace(',', ' ').split()
                    for p in parts:
                        grid_data.elements[int(p) - 1].Q_gen = 100

            for n in grid_data.nodes:
                if n.id in grid_data.bc_nodes:
                    n.BC = True

            # return global_data, grid_data

    except FileNotFoundError:
        print(f"Błąd: Plik '{file_path}' nie został znaleziony.")
        return None, None
    except Exception as e:
        print(f"Wystąpił nieoczekiwany błąd: {e}")
        return None, None
