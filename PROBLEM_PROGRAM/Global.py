# GlobalData, Grid, ElemUniv
from dataclasses import dataclass, field
from typing import List
from Element import Element

class Node:
    def __init__(self, node_id, x, y):
        self.id = int(node_id)
        self.x = float(x)
        self.y = float(y)
        self.BC = False

    def __repr__(self):
        return f"Node(ID: {self.id}, X: {self.x}, Y: {self.y}) - BC:{self.BC}\n"

@dataclass
class GlobalData:
    SimulationTime: float = 0.0
    SimulationStepTime: float = 0.0
    Conductivity: float = 0.0    
    Alfa: float = 0.0            
    Tot: float = 0.0             
    InitialTemp: float = 0.0
    Density: float = 0.0         
    SpecificHeat: float = 0.0    
    npc: int = 2 # zakladamy odniesienie do 1d

@dataclass
class Grid:
    nN: int = 0
    nE: int = 0
    nodes: List[Node] = field(default_factory=list)
    elements: List[Element] = field(default_factory=list)
    bc_nodes: List[int] = field(default_factory=list)

    gen_elements: List[int] = field(default_factory=list)
    w_f_elements: List[int] = field(default_factory=list)
    w_in_elements: List[int] = field(default_factory=list)
    w_out_elements: List[int] = field(default_factory=list)

global_data = GlobalData()
grid_data = Grid()