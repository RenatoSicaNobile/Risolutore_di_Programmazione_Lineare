import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.ticker import MaxNLocator
import logging
import traceback

def plot_solution_to_axes(ax, c, A, b, variables, show_all_quadrants, solution=None, show_solution=True):
    """Grafico avanzato: regione ammissibile, vincoli, intersezioni, gradiente, curve di livello, soluzione ottima, tutti i quadranti."""
    try:
        # 1. Limiti e aspetto grafico (sempre quadrato).
        if show_all_quadrants:
            x_min, x_max = -2, 12
            y_min, y_max = -2, 12
        else:
            x_min, x_max = 0, 12
            y_min, y_max = 0, 12
        ax.set_aspect('equal')

        # 2. Vincoli e punti di intersezione
        x = np.linspace(x_min, x_max, 1000)
        intersection_points = []

        constraint_colors = ["#1f77b4", "#1f77b4", "#1f77b4", "#ff7f0e"] # L'ultimo arancione per enfasi
        constraint_labels = [
            "Vincolo 1",
            "Vincolo 2",
            "Vincolo 3",
            "Vincolo 4"
        ]

        for i in range(len(b)):
            color = constraint_colors[i % len(constraint_colors)]
            label = constraint_labels[i] if i < len(constraint_labels) else f"Vincolo {i+1}"
            if abs(A[i, 1]) > 1e-10:
                y = (b[i] - A[i, 0] * x) / A[i, 1]
                valid = (y >= y_min) & (y <= y_max)
                if np.any(valid):
                    ax.plot(x[valid], y[valid], color=color, lw=2, label=label if i < 4 else None, zorder=3)
                # Intersezioni con assi
                try:
                    y0 = b[i] / A[i, 1]
                    if x_min <= 0 <= x_max and y_min <= y0 <= y_max:
                        intersection_points.append([0, y0])
                except ZeroDivisionError:
                    pass
                try:
                    x0 = b[i] / A[i, 0]
                    if y_min <= 0 <= y_max and x_min <= x0 <= x_max:
                        intersection_points.append([x0, 0])
                except ZeroDivisionError:
                    pass
            else:
                xval = b[i] / A[i, 0] if abs(A[i, 0]) > 1e-10 else 0
                if x_min <= xval <= x_max:
                    ax.axvline(xval, color=color, lw=2, label=label if i < 4 else None, zorder=3)
                    intersection_points.append([xval, 0])

        # 3. Intersezioni tra vincoli
        for i in range(len(b)):
            for j in range(i+1, len(b)):
                try:
                    pt = np.linalg.solve(A[[i,j], :2], b[[i,j]])
                    if (x_min-1e-6 <= pt[0] <= x_max+1e-6) and (y_min-1e-6 <= pt[1] <= y_max+1e-6):
                        intersection_points.append(pt)
                except Exception:
                    pass

        # 4. Vertici della regione ammissibile (e etichette)
        vertices = find_feasible_vertices(A, b, show_all_quadrants)
        if vertices:
            polygon = Polygon(vertices, closed=True, alpha=0.4, label='Regione ammissibile', facecolor='#a6cee3', edgecolor='none', zorder=2)
            ax.add_patch(polygon)
            for idx, (xpt, ypt) in enumerate(vertices):
                ax.scatter(xpt, ypt, color='blue', s=40, zorder=5)
                # Etichetta con offset alternato
                if idx % 2 == 0:
                    ax.annotate(f"({xpt:.1f}, {ypt:.1f})", (xpt, ypt), textcoords="offset points", xytext=(10, 5), ha='left', va='bottom',
                                fontsize=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'), zorder=6)
                else:
                    ax.annotate(f"({xpt:.1f}, {ypt:.1f})", (xpt, ypt), textcoords="offset points", xytext=(-55, -10), ha='left', va='top',
                                fontsize=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'), zorder=6)

        # 5. Gradiente e curve di livello
        if c is not None and len(c) >= 2 and np.any(c):
            # Curve di livello (più visibili)
            if solution is not None and abs(c[1]) > 1e-10:
                zval = np.dot(c, solution)
                y_level = (zval - c[0]*x) / c[1]
                valid = (y_level >= y_min) & (y_level <= y_max)
                ax.plot(x[valid], y_level[valid], 'g--', lw=2, label='Direzione ottimizzazione', zorder=9)
            # Diverse curve di livello (ben visibili)
            if solution is not None:
                offsets = np.linspace(-8, 8, 9)
                for off in offsets:
                    zval = np.dot(c, solution) + off
                    if abs(c[1]) > 1e-10:
                        y_level = (zval - c[0]*x) / c[1]
                        valid = (y_level >= y_min) & (y_level <= y_max)
                        ax.plot(x[valid], y_level[valid], color='gray', ls=(0, (4, 6)), lw=2, alpha=0.7, zorder=1)
            # Vettore gradiente (sempre dall'origine)
            gradvec = c / (np.linalg.norm(c) + 1e-10)
            arr_scale = 3.5 if not show_all_quadrants else 6
            arrow = FancyArrowPatch((0,0), gradvec*arr_scale, color='red', arrowstyle='->', mutation_scale=25, lw=3, zorder=20, label='Vettore gradiente')
            ax.add_patch(arrow)

        # 6. Soluzione ottima
        if show_solution and solution is not None and len(solution) >= 2:
            ax.plot(solution[0], solution[1], 'o', color='red', markersize=11, label='Soluzione ottima', zorder=11)

        # 7. Assi cartesiani, griglia, legenda
        ax.axhline(0, color='black', linewidth=1.5, zorder=0)
        ax.axvline(0, color='black', linewidth=1.5, zorder=0)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel(variables[0], fontsize=12)
        ax.set_ylabel(variables[1], fontsize=12)
        ax.set_title("Regione Ammissibile", fontsize=14)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(True, which='both', linestyle='--', alpha=0.6, linewidth=1.1, zorder=0)

        # Legenda custom
        handles = []
        labels = []
        for i in range(min(len(b), len(constraint_labels))):
            line = plt.Line2D([0], [0], color=constraint_colors[i], lw=2)
            handles.append(line)
            labels.append(constraint_labels[i])
        handles.append(plt.Rectangle((0,0),1,1, color='#a6cee3', alpha=0.4, ec='none'))
        labels.append("Regione ammissibile")
        handles.append(plt.Line2D([0], [0], color='red', lw=0, marker='o', markersize=11))
        labels.append("Soluzione ottima")
        handles.append(plt.Line2D([0], [0], color='green', ls='--', lw=2))
        labels.append("Direzione ottimizzazione")
        handles.append(FancyArrowPatch((0,0), (0.5,0.5), color='red', arrowstyle='->', mutation_scale=22, lw=3))
        labels.append("Vettore gradiente")
        ax.legend(handles, labels, loc='upper right', fontsize=12, frameon=True)

    except Exception as e:
        logging.error(f"Errore in plot_solution_to_axes: {e}")
        logging.error(traceback.format_exc())

def find_feasible_vertices(A, b, show_all_quadrants):
    """Trova e ordina i vertici della regione ammissibile per problemi 2D."""
    n = len(b)
    points = []
    for i in range(n):
        try:
            if abs(A[i,0]) > 1e-10:
                x0 = b[i]/A[i,0]
                if show_all_quadrants or x0 >= -1e-10:
                    points.append([x0, 0])
        except Exception:
            pass
        try:
            if abs(A[i,1]) > 1e-10:
                y0 = b[i]/A[i,1]
                if show_all_quadrants or y0 >= -1e-10:
                    points.append([0, y0])
        except Exception:
            pass

    for i in range(n):
        for j in range(i+1, n):
            try:
                pt = np.linalg.solve(A[[i,j], :2], b[[i,j]])
                if show_all_quadrants or (pt[0] >= -1e-8 and pt[1] >= -1e-8):
                    points.append(pt)
            except Exception:
                pass

    feasible = []
    for pt in points:
        if np.all(np.dot(A, pt) <= b + 1e-8):
            feasible.append(tuple(pt))
    feasible = list(dict.fromkeys(feasible))
    if len(feasible) > 2:
        center = np.mean(feasible, axis=0)
        feasible.sort(key=lambda p: np.arctan2(p[1]-center[1], p[0]-center[0]))
    return feasible