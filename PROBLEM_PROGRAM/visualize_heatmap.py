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
            node = grid.nodes[nid - 1] 
            vertices.append((node.x, node.y))
            elem_temps.append(t[nid - 1])
        
        verts.append(vertices)
        
        # 3. Oblicz średnią temperaturę elementu
        avg_temp = sum(elem_temps) / len(elem_temps)
        poly_colors.append(avg_temp)

    # 4. Tworzenie kolekcji
    p = PolyCollection(verts, cmap=cm.jet, edgecolor='black', linewidths=0.5)
    
    # Przypisanie kolorów
    p.set_array(poly_colors)
    
    # --- KLUCZOWA ZMIANA TUTAJ ---
    # Ustawiamy sztywne granice skali kolorów od 0 do 25
    p.set_clim(vmin=0, vmax=40) 
    # -----------------------------
    
    # Dodanie kolekcji do wykresu
    ax.add_collection(p)

    # 5. Dodanie paska kolorów
    # Opcja extend='max' dodaje strzałkę na górze paska,
    # sygnalizując, że są wartości wykraczające poza skalę (powyżej 25 stopni)
    cbar = fig.colorbar(p, ax=ax, extend='max') 
    cbar.set_label('Temperatura [$^\circ$C]', rotation=270, labelpad=15)

    # Ustawienia widoku
    ax.set_aspect('equal')
    ax.autoscale_view()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Rozkład temperatury')
    
    plt.tight_layout()
    plt.savefig(f'heatmap.png', dpi=300)
    # plt.show()