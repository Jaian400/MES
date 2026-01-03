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
            reading_wf = False
            reading_win = False
            reading_wout = False

            for line in lines:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('*Node'):
                    reading_nodes, reading_elements, reading_bc, reading_e_gen, reading_wf, reading_win, reading_wout = True, False, False, False, False, False, False
                    continue
                elif line.startswith('*Element'):
                    reading_nodes, reading_elements, reading_bc, reading_e_gen, reading_wf, reading_win, reading_wout = False, True, False, False, False, False, False
                    continue
                elif line.startswith('*BC'):
                    reading_nodes, reading_elements, reading_bc, reading_e_gen, reading_wf, reading_win, reading_wout = False, False, True, False, False, False, False
                    continue
                elif line.startswith('*Elset, elset=GEN'):
                    reading_nodes, reading_elements, reading_bc, reading_e_gen, reading_wf, reading_win, reading_wout = False, False, False, True, False, False, False
                    continue
                elif line.startswith('*Elset, elset=WALLS_FURNANCE'):
                    reading_nodes, reading_elements, reading_bc, reading_e_gen, reading_wf, reading_win, reading_wout = False, False, False, False, True, False, False
                    continue
                elif line.startswith('*Elset, elset=WALLS_INNER'):
                    reading_nodes, reading_elements, reading_bc, reading_e_gen, reading_wf, reading_win, reading_wout = False, False, False, False, False, True, False
                    continue
                elif line.startswith('*Elset, elset=WALLS_OUTTER'):
                    reading_nodes, reading_elements, reading_bc, reading_e_gen, reading_wf, reading_win, reading_wout = False, False, False, False, False, False, True
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
                    parts = line.replace(',', ' ').split()
                    for p in parts:
                        grid_data.bc_nodes.append(int(p))
                
                elif reading_e_gen:
                    parts = line.replace(',', ' ').split()
                    for p in parts:
                        element_id = int(p)
                        grid_data.gen_elements.append(element_id)
                        # if 0 < element_id <= len(grid_data.elements):
                        #     grid_data.elements[element_id - 1].Q_gen = 100

                elif reading_wf:
                    parts = line.replace(',', ' ').split()
                    for p in parts:
                        grid_data.w_f_elements.append(int(p))

                elif reading_win:
                    parts = line.replace(',', ' ').split()
                    for p in parts:
                        grid_data.w_in_elements.append(int(p))

                elif reading_wout:
                    parts = line.replace(',', ' ').split()
                    for p in parts:
                        grid_data.w_out_elements.append(int(p))

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
