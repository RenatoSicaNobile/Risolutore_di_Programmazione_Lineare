import numpy as np
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.ticker import MaxNLocator

def intersect_lines(a1, b1, a2, b2):
    """Return intersection point of a1[0]*x + a1[1]*y = b1 and a2[0]*x + a2[1]*y = b2, or None if parallel."""
    A = np.array([a1, a2])
    b = np.array([b1, b2])
    try:
        if np.linalg.matrix_rank(A) < 2:
            return None
        pt = np.linalg.solve(A, b)
        return pt
    except Exception:
        return None

def get_polygon_vertices(A, b, show_all_quadrants):
    """Find all feasible intersections and sort them for polygon plotting."""
    n = len(b)
    verts = []
    for i in range(n):
        for j in range(i+1, n):
            pt = intersect_lines(A[i], b[i], A[j], b[j])
            if pt is not None and np.all(np.dot(A, pt) <= b + 1e-8):
                # For all quadrants, allow negative; else only x>=0, y>=0
                if show_all_quadrants or (pt[0] >= -1e-8 and pt[1] >= -1e-8):
                    verts.append(tuple(np.round(pt, 8)))
    verts = list(set(verts))
    if len(verts) > 2:
        center = np.mean(verts, axis=0)
        verts.sort(key=lambda p: np.arctan2(p[1]-center[1], p[0]-center[0]))
    return verts

def plot_constraint_segments(ax, A, b, xlims, ylims, colors):
    """Plot constraint lines only within axis bounds, including vertical/horizontal lines."""
    n = len(b)
    for i in range(n):
        a = A[i]
        bi = b[i]
        color = colors[i % len(colors)]
        # Try to find endpoints inside plot window
        points = []
        # Intersect with left/right borders (x=xmin/xmax)
        for x in xlims:
            if abs(a[1]) > 1e-12:
                y = (bi - a[0]*x) / a[1]
                if ylims[0]-1e-8 <= y <= ylims[1]+1e-8:
                    points.append((x, y))
        # Intersect with bottom/top borders (y=ymin/ymax)
        for y in ylims:
            if abs(a[0]) > 1e-12:
                x = (bi - a[1]*y) / a[0]
                if xlims[0]-1e-8 <= x <= xlims[1]+1e-8:
                    points.append((x, y))
        # Remove duplicates
        points = [tuple(np.round(p, 8)) for p in points]
        points = list(set(points))
        # Plot only when two distinct points are found
        if len(points) >= 2:
            points = sorted(points, key=lambda p: (p[0], p[1]))
            ax.plot([p[0] for p in points], [p[1] for p in points], color=color, linewidth=2, label=f"Vincolo {i+1}")

def plot_solution_to_axes(ax, c, A, b, variables, show_all_quadrants=False, solution=None, show_solution=True):
    if len(variables) != 2 or A.shape[1] != 2:
        ax.text(0.5, 0.5, "Grafico disponibile solo per problemi con 2 variabili",
                ha="center", va="center")
        return

    # Plot limits
    if show_all_quadrants:
        x_min, x_max = -10, 10
        y_min, y_max = -10, 10
    else:
        x_min, x_max = 0, 10
        y_min, y_max = 0, 10
    xlims = (x_min, x_max)
    ylims = (y_min, y_max)
    constraint_colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:cyan', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

    plot_constraint_segments(ax, A, b, xlims, ylims, constraint_colors)

    # Feasible region
    verts = get_polygon_vertices(A, b, show_all_quadrants)
    if verts and len(verts) > 2:
        polygon = Polygon(verts, closed=True, alpha=0.3, color='skyblue', label='Regione ammissibile')
        ax.add_patch(polygon)
        for vx, vy in verts:
            ax.plot(vx, vy, 'bo')
            ax.text(vx, vy, f"({vx:.1f}, {vy:.1f})", fontsize=8, ha='right', va='bottom',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # Highlight vertical segment at x1=3 if inside plot window
    x1_segments = []
    for pt in verts:
        if abs(pt[0] - 3) < 1e-6:
            x1_segments.append(pt)
    if len(x1_segments) == 2:
        x1_segments = sorted(x1_segments, key=lambda p: p[1])
        ax.plot([3, 3], [x1_segments[0][1], x1_segments[1][1]], color='red', linewidth=3, label='Segmento x₁=3')

    # Solution and direction
    if show_solution and solution is not None and len(solution) == 2:
        ax.plot(solution[0], solution[1], 'ro', markersize=8, label='Soluzione ottima')
        # Objective direction line (isocost)
        if c is not None and len(c) == 2 and np.linalg.norm(c) > 1e-12:
            # For max: direction is -c; for min: direction is +c
            direction = -c / np.linalg.norm(c)
            ortho = np.array([direction[1], -direction[0]])
            line_len = max(x_max - x_min, y_max - y_min) * 1.5
            p1 = solution + ortho * line_len
            p2 = solution - ortho * line_len
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'g--', linewidth=2, label='Direzione ottimizzazione')
            # Gradient vector: from origin, direction of increase
            vec_length = min(x_max-x_min, y_max-y_min)*0.3
            grad_start = np.array([0, 0])
            grad_end = grad_start + direction * vec_length
            arrow = FancyArrowPatch((grad_start[0], grad_start[1]),
                                   (grad_end[0], grad_end[1]),
                                   color='red', arrowstyle='->', mutation_scale=15,
                                   linewidth=2, label='Vettore gradiente')
            ax.add_patch(arrow)

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
        seen = set()
        new_handles = []
        new_labels = []
        for h, l in zip(handles, labels):
            if l not in seen:
                new_handles.append(h)
                new_labels.append(l)
                seen.add(l)
        ax.legend(new_handles, new_labels, loc='upper right')