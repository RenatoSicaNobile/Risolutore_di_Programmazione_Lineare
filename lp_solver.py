import numpy as np

def preprocess_problem(obj_type, obj_coeffs, constraints, nn_signs, nn_values):
    """Convert input into canonical form for simplex and plotting. Returns c, A, b (numpy arrays)."""
    c = []
    for sign, coeff in obj_coeffs:
        coeff = float(coeff)
        c.append(coeff if sign == "+" else -coeff)
    if obj_type == "max":
        c = [-x for x in c]
    A = []
    b = []
    ineqs = []
    n_vars = len(c)
    for cons in constraints:
        row = []
        for sign, val in zip(cons["signs"], cons["values"]):
            coeff = float(val)
            row.append(coeff if sign == "+" else -coeff)
        A.append(row)
        b.append(float(cons["rhs"]))
        ineqs.append(cons["ineq"])
    # Non-negativity constraints
    for i, (sign, val) in enumerate(zip(nn_signs, nn_values)):
        val = float(val)
        if sign in (">=", "≥"):
            row = [0.0] * n_vars
            row[i] = -1.0
            A.append(row)
            b.append(-val)
            ineqs.append("≤")
        elif sign in ("<=", "≤"):
            row = [0.0] * n_vars
            row[i] = 1.0
            A.append(row)
            b.append(val)
            ineqs.append("≤")
    # Canonical form: all ≤
    A_std = []
    b_std = []
    for ai, bi, op in zip(A, b, ineqs):
        if op == "≤":
            A_std.append(ai)
            b_std.append(bi)
        elif op == "≥":
            A_std.append([-x for x in ai])
            b_std.append(-bi)
        elif op == "=":
            A_std.append(ai)
            b_std.append(bi)
            A_std.append([-x for x in ai])
            b_std.append(-bi)
    return np.array(c), np.array(A_std), np.array(b_std)

def simplex_solve(c, A, b):
    """Standard primal simplex with tableau, no integer constraints."""
    m, n = A.shape
    tableau = np.zeros((m+1, n+m+1))
    tableau[:-1, :n] = A
    tableau[:-1, n:n+m] = np.eye(m)
    tableau[:-1, -1] = b
    tableau[-1, :n] = c
    max_iter = 100 * (n + m)
    for _ in range(max_iter):
        if np.all(tableau[-1, :-1] >= -1e-10):
            break
        candidates = np.where(tableau[-1, :-1] < -1e-10)[0]
        if len(candidates) == 0:
            break
        entering = candidates[0]
        ratios = np.where(
            tableau[:-1, entering] > 1e-10,
            tableau[:-1, -1] / tableau[:-1, entering],
            np.inf
        )
        if np.all(ratios == np.inf):
            return None, None
        leaving = np.argmin(ratios)
        pivot = tableau[leaving, entering]
        tableau[leaving, :] /= pivot
        for i in range(m + 1):
            if i != leaving:
                tableau[i, :] -= tableau[i, entering] * tableau[leaving, :]
    solution = np.zeros(n)
    for j in range(n):
        col = tableau[:-1, j]
        if np.count_nonzero(np.abs(col - 1) < 1e-8) == 1 and np.count_nonzero(np.abs(col) < 1e-8) == m - 1:
            row = np.where(np.abs(col - 1) < 1e-8)[0][0]
            solution[j] = tableau[row, -1]
    optimal_value = tableau[-1, -1]
    return solution, optimal_value