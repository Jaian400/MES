import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# gemini
def visualize_grid(grid):
    fig, ax = plt.subplots(figsize=(8, 6))
    for elem in grid.elements:
        vertices = []
        for nid in elem.node_ids:
            node = grid.nodes[nid - 1]
            vertices.append((node.x, node.y))
        
        poly = Polygon(vertices, closed=True, edgecolor='blue', facecolor='cyan', alpha=0.3, linewidth=1.5)
        ax.add_patch(poly)
        
        center_x = sum([v[0] for v in vertices]) / 4
        center_y = sum([v[1] for v in vertices]) / 4
        ax.text(center_x, center_y, str(elem.id), color='blue', fontsize=10, ha='center', va='center', fontweight='bold')

    for node in grid.nodes:
        ax.plot(node.x, node.y, 'ko', markersize=4)
        ax.text(node.x, node.y, f" {node.id}", color='black', fontsize=8, va='bottom')

    if grid.bc_nodes:
        bc_x = []
        bc_y = []
        for nid in grid.bc_nodes:
            node = grid.nodes[nid - 1]
            bc_x.append(node.x)
            bc_y.append(node.y)
        ax.plot(bc_x, bc_y, 'ro', markersize=8, label='Warunek Brzegowy (BC)')

    ax.set_aspect('equal') 
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Wizualizacja Siatki MES')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend()
    
    ax.autoscale_view()
    
    plt.tight_layout()
    plt.savefig('visual.png')
    # plt.show()