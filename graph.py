import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull

def get_feasible_points(A, b, x1_min, x1_max, x2_min, x2_max, show_all_quadrants):
    """Compute all feasible intersection points given constraints."""
    m, n = A.shape
    points = []
    # All pairwise intersections
    for i in range(m):
        for j in range(i+1, m):
            Ai = np.array([A[i], A[j]])
            if np.linalg.matrix_rank(Ai) < 2:
                continue
            try:
                pt = np.linalg.solve(Ai, np.array([b[i], b[j]]))
                if not np.isfinite(pt).all():
                    continue
                points.append(pt)
            except Exception:
                continue
    # Intersections with axes
    for i in range(m):
        a, bb = A[i], b[i]
        if abs(a[0]) > 1e-8:
            points.append(np.array([bb / a[0], 0]))
        if abs(a[1]) > 1e-8:
            points.append(np.array([0, bb / a[1]]))
    # Origin
    points.append(np.array([0, 0]))
    # Filter feasible points
    feasible_points = []
    for pt in points:
        if not np.isfinite(pt).all():
            continue
        if np.all(A @ pt - b <= 1e-7) and (show_all_quadrants or (pt[0] >= -1e-7 and pt[1] >= -1e-7)):
            feasible_points.append(pt)
    return np.unique(np.round(feasible_points, 8), axis=0)

def find_optimal_segment(c, feasible_points, z_opt):
    """Return segment endpoints of optimal solutions (if multiple) for alternate optima."""
    eps = 1e-7
    optimal_points = [pt for pt in feasible_points if abs(np.dot(c, pt) - z_opt) <= eps]
    if len(optimal_points) > 2:
        hull = ConvexHull(optimal_points)
        optimal_points = [optimal_points[i] for i in hull.vertices]
    return np.array(optimal_points)

def place_label(ax, pt, placed_labels, min_dist=0.8):
    """Place label for pt, avoiding overlap with others in placed_labels."""
    offsets = [
        (0.25, 0.25), (-0.95, 0.25), (0.25, -0.55), (-0.95, -0.55),
        (0.0, 0.7), (0.7, 0.0), (0.0, -0.7), (-0.7, 0.0),
        (0.55, 0.55), (-0.55, 0.55), (0.55, -0.55), (-0.55, -0.55),
    ]
    for dx, dy in offsets:
        new_pt = (pt[0] + dx, pt[1] + dy)
        overlap = False
        for lab in placed_labels:
            if np.linalg.norm(np.array(new_pt) - np.array(lab)) < min_dist:
                overlap = True
                break
        if not overlap:
            placed_labels.append(new_pt)
            ax.text(new_pt[0], new_pt[1], f"({pt[0]:.1f}, {pt[1]:.1f})",
                    fontsize=10, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            return
    # fallback
    ax.text(pt[0], pt[1], f"({pt[0]:.1f}, {pt[1]:.1f})",
            fontsize=10, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

def plot_solution_to_axes(
    ax, c, A, b, variables, show_all_quadrants=False, solution=None, show_solution=True
):
    """
    Fully rewritten: Plots constraints, feasible region, all intersections, level curves (with legend),
    optimal solution, optimal segment/direction (dotted), gradient (with legend), axes (if show_all_quadrants).
    """
    if len(c) != 2 or A.shape[1] != 2:
        ax.text(0.5, 0.5, "Plotting only for 2 variables", ha="center", va="center", fontsize=16)
        return

    # Plot range
    if show_all_quadrants:
        x1_min, x1_max = -10, 10
        x2_min, x2_max = -10, 10
    else:
        x1_min, x1_max = 0, 10
        x2_min, x2_max = 0, 10

    x1 = np.linspace(x1_min, x1_max, 500)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

    # ----------- Plot constraints ----------------
    for i, (a, bb) in enumerate(zip(A, b)):
        if abs(a[1]) > 1e-8:
            y = (bb - a[0]*x1) / a[1]
            ax.plot(x1, y, color=colors[i % len(colors)], label=f'Vincolo {i+1}', linewidth=2)
        elif abs(a[0]) > 1e-8:
            x_vert = bb / a[0]
            ax.axvline(x_vert, color=colors[i % len(colors)], label=f'Vincolo {i+1}', linewidth=2)

    # ----------- Feasible region -----------------
    feasible_points = get_feasible_points(A, b, x1_min, x1_max, x2_min, x2_max, show_all_quadrants)
    if len(feasible_points) > 2:
        hull = ConvexHull(feasible_points)
        region_pts = feasible_points[hull.vertices]
        poly = Polygon(region_pts, closed=True, facecolor='cornflowerblue', alpha=0.2, label="Regione ammissibile")
        ax.add_patch(poly)
    else:
        region_pts = feasible_points

    # ----------- Intersection points and labels (no overlap) -----------------
    placed_labels = []
    for pt in feasible_points:
        ax.plot(pt[0], pt[1], 'o', color='blue')
        place_label(ax, pt, placed_labels, min_dist=0.85)

    # ----------- Level curves (grey dotted), legend entry -------------
    level_curve_handle = None
    if solution is not None:
        z_opt = float(np.dot(c, solution))
        levels = np.linspace(z_opt - 6, z_opt + 6, 7)
        for i, lev in enumerate(levels):
            if abs(c[1]) > 1e-8:
                y = (lev - c[0]*x1) / c[1]
                l, = ax.plot(x1, y, ls=":", color="grey", alpha=0.7, lw=1, zorder=1)
                if i == 3:  # Use middle for legend handle
                    level_curve_handle = l
            elif abs(c[0]) > 1e-8:
                x = lev / c[0]
                l = ax.axvline(x, ls=":", color="grey", alpha=0.7, lw=1, zorder=1)
                if i == 3:
                    level_curve_handle = l

    # ----------- Gradient vector from origin -----------------
    grad = np.array(c)
    grad_norm = np.linalg.norm(grad)
    if grad_norm > 1e-8:
        grad_unit = grad / grad_norm
        grad_scale = min(x1_max, x2_max) * 0.5
        grad_vec = grad_unit * grad_scale
        ax.arrow(0, 0, grad_vec[0], grad_vec[1],
                 head_width=0.25, head_length=0.4, fc='red', ec='red', linewidth=2, length_includes_head=True, zorder=3)

    # ----------- Optimal segment (dotted) or direction (dotted) -------------
    if solution is not None:
        z_opt = float(np.dot(c, solution))
        opt_pts = find_optimal_segment(c, feasible_points, z_opt)
        if len(opt_pts) >= 2:
            order = np.argsort(opt_pts[:, 0])
            opt_pts = opt_pts[order]
            ax.plot(opt_pts[:, 0], opt_pts[:, 1], 'g:', linewidth=4, label="Direzione ottimizzazione", zorder=4)
        else:
            # Show optimization direction as dotted green line
            dir_vec = grad / grad_norm if grad_norm > 1e-8 else np.array([1, 0])
            t = np.linspace(-10, 10, 100)
            line_x = solution[0] + dir_vec[0] * t
            line_y = solution[1] + dir_vec[1] * t
            ax.plot(line_x, line_y, "g:", linewidth=2, label="Direzione ottimizzazione", zorder=2)

    # ----------- Optimal solution point ---------------
    if show_solution and solution is not None:
        ax.plot(solution[0], solution[1], 'o', color='red', markersize=10, label="Soluzione ottima", zorder=5)

    # ----------- Draw axes for quadrants -----------------
    if show_all_quadrants:
        ax.axhline(0, color='black', linewidth=1.3, linestyle='-', alpha=0.7, zorder=0)
        ax.axvline(0, color='black', linewidth=1.3, linestyle='-', alpha=0.7, zorder=0)

    # ----------- Axes/labels/title ----------------------
    ax.set_xlabel(variables[0])
    ax.set_ylabel(variables[1])
    ax.set_xlim(x1_min, x1_max)
    ax.set_ylim(x2_min, x2_max)
    ax.set_title("Regione Ammissibile")

    # ----------- Custom legend (with level curves) --------------
    legend_elements = []
    for i in range(len(A)):
        legend_elements.append(Line2D([0], [0], color=colors[i % len(colors)], lw=2, label=f'Vincolo {i+1}'))
    legend_elements.append(Line2D([0], [0], color='cornflowerblue', lw=8, alpha=0.2, label='Regione ammissibile'))
    legend_elements.append(Line2D([0], [0], marker='o', color='red', lw=0, markersize=10, label='Soluzione ottima'))
    legend_elements.append(Line2D([0], [0], color='green', lw=4, ls=':', label='Direzione ottimizzazione'))
    legend_elements.append(Line2D([0], [0], color='red', lw=2, label='Vettore gradiente'))
    if level_curve_handle is not None:
        legend_elements.append(Line2D([0], [0], color='grey', lw=2, ls=':', label='Curve di livello'))
    ax.legend(handles=legend_elements, loc='best')

    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_aspect('auto')