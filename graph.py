import numpy as np
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.ticker import MaxNLocator

def plot_solution_to_axes(ax, c, A, b, variables, show_all_quadrants=False, solution=None, show_solution=True):
    if len(variables) != 2 or A.shape[1] != 2:
        ax.text(0.5, 0.5, "Grafico disponibile solo per problemi con 2 variabili",
                ha="center", va="center")
        return

    # Set plot limits
    if show_all_quadrants:
        # Show -10 to +10 on both axes, origin at center
        x_min, x_max = -10, 10
        y_min, y_max = -10, 10
    else:
        x_min, x_max = 0, 10
        y_min, y_max = 0, 10

    x = np.linspace(x_min, x_max, 400)
    intersection_points = []

    for i in range(len(b)):
        if np.abs(A[i, 1]) > 1e-12:
            y = (b[i] - A[i, 0] * x) / A[i, 1]
            ax.plot(x, y, label=f"Vincolo {i+1}")
            y_intercept = b[i] / A[i, 1] if np.abs(A[i, 1]) > 1e-12 else 0
            intersection_points.append((0, y_intercept))
            if np.abs(A[i, 0]) > 1e-12:
                x_intercept = b[i] / A[i, 0]
                intersection_points.append((x_intercept, 0))
        else:
            x_val = b[i] / A[i, 0] if np.abs(A[i, 0]) > 1e-12 else 0
            ax.axvline(x_val, label=f"Vincolo {i+1}")
            intersection_points.append((x_val, 0))

    for i in range(len(b)):
        for j in range(i+1, len(b)):
            Ai, Aj = A[i, :2], A[j, :2]
            bi, bj = b[i], b[j]
            det = Ai[0]*Aj[1] - Ai[1]*Aj[0]
            if np.abs(det) > 1e-10:
                point = np.linalg.solve(np.vstack([Ai, Aj]), np.array([bi, bj]))
                intersection_points.append(point)

    vertices = find_feasible_vertices(A, b, show_all_quadrants=show_all_quadrants)
    if vertices:
        polygon = Polygon(vertices, closed=True, alpha=0.3, label='Regione ammissibile')
        ax.add_patch(polygon)

    if show_solution and solution is not None and len(solution) == 2:
        ax.plot(solution[0], solution[1], 'ro', markersize=8, label='Soluzione ottima')
        if c is not None and len(c) == 2 and np.linalg.norm(c) > 1e-12 and c[1] != 0:
            obj_slope = -c[0] / c[1]
            obj_line_x = np.array([x_min, x_max])
            obj_line_y = obj_slope * (obj_line_x - solution[0]) + solution[1]
            ax.plot(obj_line_x, obj_line_y, 'g--', linewidth=2, label='Direzione ottimizzazione')
            for level in np.linspace(0.2, 1.0, 4):
                offset = level * 2
                ax.plot(obj_line_x, obj_line_y + offset, 'g:', linewidth=0.5, alpha=0.5)
                ax.plot(obj_line_x, obj_line_y - offset, 'g:', linewidth=0.5, alpha=0.5)
            vec_length = min(x_max-x_min, y_max-y_min) * 0.2
            vec_x = vec_length * c[0] / np.linalg.norm(c)
            vec_y = vec_length * c[1] / np.linalg.norm(c)
            arrow = FancyArrowPatch((solution[0], solution[1]),
                                   (solution[0] + vec_x, solution[1] + vec_y),
                                   color='red', arrowstyle='->', mutation_scale=15,
                                   linewidth=2, label='Vettore gradiente')
            ax.add_patch(arrow)

    if intersection_points:
        points = np.array(intersection_points)
        ax.scatter(points[:, 0], points[:, 1], color='blue', zorder=5)
        for point in points:
            ax.text(point[0], point[1], f"({point[0]:.1f}, {point[1]:.1f})",
                    fontsize=8, ha='right', va='bottom',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel(variables[0])
    ax.set_ylabel(variables[1])
    ax.set_title("Regione Ammissibile")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc='upper right')

def find_feasible_vertices(A, b, show_all_quadrants=False):
    if A.shape[1] != 2:
        return None
    n = len(b)
    vertices = []
    min_val = -10 if show_all_quadrants else 0
    vertices.append((min_val, min_val))
    for i in range(n):
        for j in range(i+1, n):
            a1, a2 = A[i, 0], A[i, 1]
            b1 = b[i]
            c1, c2 = A[j, 0], A[j, 1]
            b2 = b[j]
            det = a1 * c2 - a2 * c1
            if abs(det) > 1e-10:
                x = (c2 * b1 - a2 * b2) / det
                y = (a1 * b2 - c1 * b1) / det
                if show_all_quadrants or (x >= -1e-10 and y >= -1e-10):
                    vertices.append((x, y))
    feasible_vertices = []
    for x, y in vertices:
        feasible = True
        for i in range(n):
            if A[i, 0]*x + A[i, 1]*y > b[i] + 1e-10:
                feasible = False
                break
        if feasible:
            feasible_vertices.append((x, y))
    if len(feasible_vertices) > 2:
        center = np.mean(feasible_vertices, axis=0)
        feasible_vertices.sort(key=lambda p: np.arctan2(p[1]-center[1], p[0]-center[0]))
    return feasible_vertices