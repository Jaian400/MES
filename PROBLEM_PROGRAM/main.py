import numpy as np
from integrate import V_ff
from visualize_grid import visualize_grid
from visualize_heatmap import visualize_heatmap
from scipy.integrate import dblquad
from scipy.sparse import csc_array, csr_array
from scipy.sparse.linalg import spsolve
from tabulate import tabulate

from Global import global_data, grid_data
from load_data import load_data
from ElemUniv import ElemUniv

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import time

#agregacja
class MatrixGlobal:
    def __init__(self):
        self.H = np.zeros((grid_data.nN, grid_data.nN))
        self.P = np.zeros((grid_data.nN))
        self.C = np.zeros((grid_data.nN, grid_data.nN))
    
    def calculate(self):
        for e in grid_data.elements:
            small_H = e.H
            small_P = e.P
            small_C = e.C
            for i in range(4):
                for j in range(4):
                    self.H[e.node_ids[i] - 1][e.node_ids[j] - 1] += small_H[i][j]
                    self.C[e.node_ids[i] - 1][e.node_ids[j] - 1] += small_C[i][j]
                
                self.P[e.node_ids[i] - 1] += small_P[i]

    def __repr__(self):
        result = "\nMacierz duza H:\n"
        result += tabulate(self.H, floatfmt=".3f", tablefmt="plain")
        result += f"\nMacierz duza P:\n{self.P}\n"
        result += "\nMacierz duza C:\n"
        result += tabulate(self.C, floatfmt=".3f", tablefmt="plain")

        return result

class SystemOfEquation:
    def __init__(self, matrix_H, nN):
        self.H = matrix_H.H
        self.P = matrix_H.P
        self.C = matrix_H.C
        self.t = np.zeros(nN)
    
    def solve(self, t0, delta_tau):
        # print(t0)
        new_C = self.C / delta_tau
        left = self.H + new_C
        # print(f"{left}")
        right = self.P + np.inner(t0, new_C)
        # print(f"{right}")
        self.t = np.linalg.solve(left, right)
        # self.t = spsolve(csc_array(left), csr_array(right))
    
    def __repr__(self):
        return f"\nt1: {self.t}\n"

def calculate_elements(elements, nodes, elem_univ):
    t_ot = global_data.Tot
    alfa = global_data.Alfa
    for e in elements:
        if e.material_type == 0: # POWIETRZE
            cond = 2.0 # uwazac z tym, standard: 0.026
            dens = 1.2      
            spec_heat = 1005 
        elif e.material_type == 1: # PIEC
            cond = 1.5     
            dens = 2200     
            spec_heat = 1000  
        elif e.material_type == 2: # SC WEWN
            cond = 0.8      
            dens = 1600     
            spec_heat = 850   
        elif e.material_type == 3: # SC ZEWN
            cond = 1.2    
            dens = 2100     
            spec_heat = 840 

        e.calculate_jacobians(nodes, elem_univ)
        e.calculate_H(cond, elem_univ, dens, spec_heat)
        e.calculate_Hbc(alfa, elem_univ, nodes, t_ot)

        if e.Q_gen != None and e.Q_gen != 0:
            e.calculate_Q_gen(elem_univ)

def calculate_state(matrixes, t0):
    system_of_equation = SystemOfEquation(matrixes, grid_data.nN)
    system_of_equation.solve(t0, global_data.SimulationStepTime)
    # print(system_of_equation)

    return system_of_equation.t

if __name__ == "__main__":
    # file_path = "Test_q_gen_big.txt"
    # file_path = "Test_q_gen.txt"
    file_path = "Problem_dom.txt"

    load_data(file_path)

    npc = 2
    global_data.npc = npc

    if global_data and grid_data:
        print("\n--- WSPÓŁRZĘDNE WSZYSTKICH WĘZŁÓW ---")
        for node_obj in grid_data.nodes:
            print(f"Węzeł ID: {node_obj.id}, Współrzędne (X, Y): ({node_obj.x}, {node_obj.y})")

        print("\n\n--- ID WĘZŁÓW DLA POSZCZEGÓLNYCH ELEMENTÓW ---")
        for elem_obj in grid_data.elements:
            print(f"Element ID: {elem_obj.id}, Q_gen: {elem_obj.Q_gen}, ID węzłów: {elem_obj.node_ids}")

        print("\n\n--- BORDER CONDITION ---")
        for abc in grid_data.bc_nodes:
            print(f"{abc}")
        
        # print(grid_data)
        # exit(0)
    else:
        print("Nie otrzymano poprawnych danych z load_data\n")
        exit(-1)
    
    elem_univ = ElemUniv(global_data.npc)
    elem_univ.calculate_elem_univ()

    # print(elem_univ)

    # OKRESLIC TYP MATERIALU, reszta 0 - powietrze
    for e_id in grid_data.w_f_elements:
        grid_data.elements[e_id - 1].material_type = 1 # 1 - sciana pieca
    for e_id in grid_data.w_in_elements:
        grid_data.elements[e_id - 1].material_type = 2 # 2 - sciana wewnetrzna
    for e_id in grid_data.w_out_elements:
        grid_data.elements[e_id - 1].material_type = 3 # 3 - sciana zewnetrzna

    # WIELKOSC GRZANIA
    for e_id in grid_data.gen_elements:
        grid_data.elements[e_id - 1].Q_gen = 5000
    
    calculate_elements(grid_data.elements, grid_data.nodes, elem_univ)
    # print(grid_data.elements)

    matrixes = MatrixGlobal()
    matrixes.calculate()
    # print(matrixes)

    # ZBIERANE DANE ------------------

    air_nodes_indices = set()
    for e in grid_data.elements:
        if e.material_type == 0: 
            for node_id in e.node_ids:
                air_nodes_indices.add(node_id - 1) 
    air_nodes_indices = list(air_nodes_indices)

    history_data = []

    # SYMULACJA ------------------

    taus = range(0, global_data.SimulationTime, global_data.SimulationStepTime)
    t0 = np.ones((grid_data.nN)) * global_data.InitialTemp
     
    start = time.time()
    for tau in taus:
        t0 = calculate_state(matrixes, t0)
        print(f"{tau+global_data.SimulationStepTime:02} s | Max: {t0.max():.5f} | Min: {t0.min():.5f} | Avg: {t0.mean():.5f}")

        avg_air_temp = t0[air_nodes_indices].mean()
        history_data.append({
            "Time": tau+global_data.SimulationStepTime,
            "Temp_Avg_Powietrze": avg_air_temp,
            "Temp_Max": t0.max(),
            "Temp_Min": t0.min()
        })

    end = time.time()

    print("Czas obliczeń:", end - start, "s")

    visualize_heatmap(grid_data, t0)

    df = pd.DataFrame(history_data)
    
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))
    
    # Wykres średniej temperatury powietrza
    plot = sns.lineplot(data=df, x="Time", y="Temp_Avg_Powietrze", color="blue", linewidth=2.5)

    plt.title("Zmiana temperatury powietrza w czasie", fontsize=15)
    plt.xlabel("Czas [s]", fontsize=12)
    plt.ylabel("Temperatura [st.C]", fontsize=12)
    plt.savefig(f'temp_in_time.png', dpi=300)

    # visualize_grid(grid_data)
    # print(t0)
    # print(global_data)
    exit(0)

    # print("\n\nCALKOWANIE")
    # calka = V_ff()

    # f1 = lambda x:5*np.pow(x, 2) + 3*x + 6
    # f2 = lambda y, x: 5 * np.power(x,2) * np.power(y,2) + 3*x*y + 6

    # calka.quadrature(f2, -1, 1)
    # print(calka.output)

    # dblquad_, error = dblquad(f2, -100, 100, lambda x: -1, lambda x: 2)
    # print(f"dblquad: {dblquad_}")
    # print(f"roznica miedzy V_ff a dblquad: {100 * abs(calka.output - dblquad_)/calka.output} %")