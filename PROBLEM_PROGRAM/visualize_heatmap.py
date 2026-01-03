import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.cm as cm
    
import matplotlib.patheffects as path_effects

def visualize_heatmap(grid, t):
    fig, ax = plt.subplots(figsize=(10, 8))

    verts = []
    poly_colors = []

    for elem in grid.elements:
        # 1. Pobierz współrzędne węzłów elementu
        vertices = []
        elem_temps = []
        
        for nid in elem.node_ids:
            node = grid.nodes[nid - 1] # Zakładamy indeksowanie od 1
            vertices.append((node.x, node.y))
            
            # 2. Pobierz temperaturę dla tego węzła z listy t
            elem_temps.append(t[nid - 1])
        
        verts.append(vertices)
        
        # 3. Oblicz średnią temperaturę elementu (kolor powierzchni)
        avg_temp = sum(elem_temps) / len(elem_temps)
        poly_colors.append(avg_temp)

    # 4. Tworzenie kolekcji wielokątów (efektywne rysowanie)
    # cmap='jet' to klasyczna skala (niebieski-zimny, czerwony-gorący)
    p = PolyCollection(verts, cmap=cm.jet, edgecolor='black', linewidths=0.5)
    
    # Przypisanie kolorów na podstawie średnich temperatur
    p.set_array(poly_colors)
    
    # Dodanie kolekcji do wykresu
    ax.add_collection(p)

    # 5. Dodanie paska kolorów (Colorbar)
    cbar = fig.colorbar(p, ax=ax)
    cbar.set_label('Temperatura [$^\circ$C]', rotation=270, labelpad=15)

    # 6. Opcjonalne: Rysowanie numerów węzłów (tylko jeśli siatka jest mała)
    if len(grid.nodes) < 100:
        for node in grid.nodes:
            ax.text(node.x, node.y, f"{t[node.id-1]:.1f}", 
                    fontsize=7, ha='center', va='center', color='white',
                    path_effects=[path_effects.withStroke(linewidth=2, foreground='black')])

    # Ustawienia widoku
    ax.set_aspect('equal')
    ax.autoscale_view()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Rozkład temperatury (MES Heatmap)')
    
    plt.tight_layout()
    plt.savefig(f'heatmap{0}.png', dpi=300)
    # plt.show()