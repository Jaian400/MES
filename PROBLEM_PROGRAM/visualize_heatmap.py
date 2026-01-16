import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.cm as cm
import matplotlib.patheffects as path_effects

def visualize_heatmap(grid, t):
    fig, ax = plt.subplots(figsize=(10, 8))

    verts = []
    poly_colors = []

    for elem in grid.elements:
        vertices = []
        elem_temps = []
        
        for nid in elem.node_ids:
            node = grid.nodes[nid - 1] 
            vertices.append((node.x, node.y))
            elem_temps.append(t[nid - 1])
        
        verts.append(vertices)
        
        avg_temp = sum(elem_temps) / len(elem_temps)
        poly_colors.append(avg_temp)

    p = PolyCollection(verts, cmap=cm.jet, edgecolor='black', linewidths=0.5)
    
    p.set_array(poly_colors)
    
    p.set_clim(vmin=0, vmax=40) 
    # -----------------------------
    
    ax.add_collection(p)

    cbar = fig.colorbar(p, ax=ax, extend='max') 
    cbar.set_label('Temperatura [$^\circ$C]', rotation=270, labelpad=15)

    ax.set_aspect('equal')
    ax.autoscale_view()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Rozkład temperatury')
    
    plt.tight_layout()
    plt.savefig(f'heatmap.png', dpi=300)
    # plt.show()