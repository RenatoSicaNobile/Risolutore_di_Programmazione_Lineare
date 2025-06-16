import numpy as np

def format_tableau(tableau, n_vars, n_slacks):
    """
    Nicely formats a simplex tableau as text.
    Args:
        tableau: The simplex tableau as a numpy array.
        n_vars: Number of decision variables.
        n_slacks: Number of slack variables.
    Returns:
        A string representation of the tableau.
    """
    lines = []
    header = ["Base"] + [f"x{i+1}" for i in range(n_vars)] + [f"s{i+1}" for i in range(n_slacks)] + ["RHS"]
    lines.append("  ".join(f"{h:>8}" for h in header))
    for i in range(n_slacks):
        row = [f"s{i+1}"]
        for j in range(n_vars + n_slacks + 1):
            row.append(f"{tableau[i, j]:8.4f}")
        lines.append("  ".join(row))
    cost_row = ["Cost"]
    for j in range(n_vars + n_slacks + 1):
        cost_row.append(f"{tableau[-1, j]:8.4f}")
    lines.append("  ".join(cost_row))
    return "\n".join(lines) + "\n"