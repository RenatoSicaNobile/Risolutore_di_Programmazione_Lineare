import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull

def get_feasible_points(A, b, x1_min=0, x1_max=10, x2_min=0, x2_max=10, show_all_quadrants=False):
    m, n = A.shape
    points = []
    # Collect all intersections of pairs of constraints
    for i in range(m):
        for j in range(i+1, m):
            AA = np.array([A[i], A[j]])
            if np.linalg.matrix_rank(AA) < 2:
                continue
            try:
                pt = np.linalg.solve(AA, np.array([b[i], b[j]]))
                if not np.isfinite(pt).all():
                    continue
                points.append(pt)
            except Exception:
                continue
    # Add axis intercepts
    for i in range(m):
        a, bb = A[i], b[i]
        if abs(a[0]) > 1e-8:
            points.append(np.array([bb / a[0], 0]))
        if abs(a[1]) > 1e-8:
            points.append(np.array([0, bb / a[1]]))
    points.append(np.array([0, 0]))
    # Filter feasible points
    feasible_points = []
    for pt in points:
        if not np.isfinite(pt).all():
            continue
        if np.all(A @ pt - b <= 1e-7) and \
           (show_all_quadrants or (pt[0] >= -1e-7 and pt[1] >= -1e-7)):
            feasible_points.append(pt)
    return np.unique(np.round(feasible_points, 8), axis=0)

def find_optimal_segment(c, feasible_points, z_opt):
    """Find segment(s) of feasible_points where c.x == z_opt (for alternate optima)"""
    eps = 1e-7
    optimal_points = [pt for pt in feasible_points if abs(np.dot(c, pt) - z_opt) <= eps]
    # Order segment endpoints (for display)
    if len(optimal_points) > 2:
        hull = ConvexHull(optimal_points)
        optimal_points = [optimal_points[i] for i in hull.vertices]
    return np.array(optimal_points)

def plot_solution_to_axes(ax, c, A, b, variables, show_all_quadrants=False, solution=None, show_solution=True):
    # Only for 2 variables
    if len(c) != 2 or A.shape[1] != 2:
        ax.text(0.5, 0.5, "Plotting only for 2 variables", ha="center", va="center", fontsize=16)
        return

    # Set up plot range
    x1_min, x1_max = 0, 10
    x2_min, x2_max = 0, 10
    if show_all_quadrants:
        x1_min, x2_min = -10, -10

    x1 = np.linspace(x1_min, x1_max, 500)
    x2 = np.linspace(x2_min, x2_max, 500)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

    # Plot constraints
    for i, (a, bb) in enumerate(zip(A, b)):
        if abs(a[1]) > 1e-8:
            y = (bb - a[0]*x1) / a[1]
            ax.plot(x1, y, color=colors[i % len(colors)], label=f'Vincolo {i+1}', linewidth=2)
        elif abs(a[0]) > 1e-8:
            x_vert = bb / a[0]
            ax.axvline(x_vert, color=colors[i % len(colors)], label=f'Vincolo {i+1}', linewidth=2)

    # Feasible region
    feasible_points = get_feasible_points(A, b, x1_min, x1_max, x2_min, x2_max, show_all_quadrants)
    if len(feasible_points) > 2:
        hull = ConvexHull(feasible_points)
        region_pts = feasible_points[hull.vertices]
        poly = Polygon(region_pts, closed=True, facecolor='cornflowerblue', alpha=0.2, label="Regione ammissibile")
        ax.add_patch(poly)
    else:
        region_pts = feasible_points

    # Draw intersection points with labels
    for pt in feasible_points:
        ax.plot(pt[0], pt[1], 'o', color='blue')
        # Offset label to avoid overlap
        offset_x = 0.15 if pt[0] >= (x1_max-x1_min)*0.7 else -0.7
        offset_y = 0.18 if pt[1] >= (x2_max-x2_min)*0.7 else -0.35
        ax.text(pt[0]+offset_x, pt[1]+offset_y, f"({pt[0]:.1f}, {pt[1]:.1f})", fontsize=10,
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # Plot level curves (grey dotted)
    if solution is not None:
        z_opt = float(np.dot(c, solution))
        levels = np.linspace(z_opt - 6, z_opt + 6, 7)
        for lev in levels:
            if abs(c[1]) > 1e-8:
                y = (lev - c[0]*x1) / c[1]
                ax.plot(x1, y, ls=":", color="grey", alpha=0.5, lw=1, zorder=1)
            elif abs(c[0]) > 1e-8:
                x = lev / c[0]
                ax.axvline(x, ls=":", color="grey", alpha=0.5, lw=1, zorder=1)

    # Plot gradient vector from origin
    grad = np.array(c)
    grad_norm = np.linalg.norm(grad)
    if grad_norm > 1e-8:
        grad_unit = grad / grad_norm
        grad_scale = min(x1_max, x2_max) * 0.5
        grad_vec = grad_unit * grad_scale
        ax.arrow(0, 0, grad_vec[0], grad_vec[1],
                 head_width=0.25, head_length=0.4, fc='red', ec='red', linewidth=2, length_includes_head=True, zorder=3)

    # Plot optimal segment (if multiple optima)
    if solution is not None:
        z_opt = float(np.dot(c, solution))
        opt_pts = find_optimal_segment(c, feasible_points, z_opt)
        if len(opt_pts) >= 2:
            order = np.argsort(opt_pts[:, 0])  # Sort by x for visual clarity
            opt_pts = opt_pts[order]
            ax.plot(opt_pts[:, 0], opt_pts[:, 1], 'g-', linewidth=4, label="Direzione ottimizzazione", zorder=4)
        else:
            # Plot optimization direction line (old style, dashed)
            dir_vec = grad / grad_norm if grad_norm > 1e-8 else np.array([1, 0])
            t = np.linspace(-10, 10, 100)
            line_x = solution[0] + dir_vec[0] * t
            line_y = solution[1] + dir_vec[1] * t
            ax.plot(line_x, line_y, "g--", linewidth=2, label="Direzione ottimizzazione", zorder=2)

    # Plot optimal solution
    if show_solution and solution is not None:
        ax.plot(solution[0], solution[1], 'o', color='red', markersize=10, label="Soluzione ottima", zorder=5)

    # Axes/labels/legend
    ax.set_xlabel(variables[0])
    ax.set_ylabel(variables[1])
    ax.set_xlim(x1_min, x1_max)
    ax.set_ylim(x2_min, x2_max)
    ax.set_title("Regione Ammissibile")
    # Custom legend
    legend_elements = []
    for i in range(len(A)):
        legend_elements.append(Line2D([0], [0], color=colors[i % len(colors)], lw=2, label=f'Vincolo {i+1}'))
    legend_elements.append(Line2D([0], [0], color='cornflowerblue', lw=8, alpha=0.2, label='Regione ammissibile'))
    legend_elements.append(Line2D([0], [0], marker='o', color='red', lw=0, markersize=10, label='Soluzione ottima'))
    legend_elements.append(Line2D([0], [0], color='green', lw=4, label='Direzione ottimizzazione'))
    legend_elements.append(Line2D([0], [0], color='red', lw=2, label='Vettore gradiente'))
    ax.legend(handles=legend_elements, loc='best')

    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_aspect('auto')