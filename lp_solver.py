import numpy as np
import utils

def run_simplex(
    c, A, b, obj_type, original_constraint_count, constraints, steps_text=None
):
    """
    Implementation of the simplex algorithm with detailed steps.
    Parameters:
        c: 1D array, coefficients of the objective function.
        A: 2D array, constraint matrix.
        b: 1D array, right-hand side values.
        obj_type: 'max' or 'min'
        original_constraint_count: int, number of original constraints.
        constraints: list of tuples, each representing a constraint.
        steps_text: Tkinter Text widget or None; outputs steps if provided.
    Returns:
        solution: optimal variable values
        slack_values: slack variable values
        optimal_value: optimal value of the objective function
        final_dual: dual variable values (shadow prices)
        problem_history: list of (description, tableau) tuples
    """
    m, n = A.shape

    # Standard form: always minimize, so negate c if maximizing
    c = np.array(c, dtype=float)
    if obj_type == "max":
        c = -c

    # Build initial tableau
    tableau = np.zeros((m + 1, n + m + 1))
    tableau[:-1, :n] = A
    tableau[:-1, n:n + m] = np.eye(m)
    tableau[:-1, -1] = b
    tableau[-1, :n] = c

    problem_history = [("Initial Tableau", tableau.copy())]

    if steps_text is not None:
        steps_text.insert("end", "\n=== INITIAL TABLEAU ===\n")
        steps_text.insert("end", utils.format_tableau(tableau, n, m))

    max_iter = 100 * (n + m)
    iter_count = 0
    while iter_count < max_iter:
        # Check for optimality
        if np.all(tableau[-1, :-1] >= -1e-10):
            if steps_text is not None:
                steps_text.insert("end", "\n=== OPTIMALITY CONDITION SATISFIED ===\n")
            break

        # Pick entering variable (Bland's rule: smallest index with negative cost)
        candidates = np.where(tableau[-1, :-1] < -1e-10)[0]
        if len(candidates) == 0:
            break
        entering = candidates[0]

        # Minimum ratio test
        ratios = np.where(
            tableau[:-1, entering] > 1e-10,
            tableau[:-1, -1] / tableau[:-1, entering],
            np.inf
        )
        if np.all(ratios == np.inf):
            if steps_text is not None:
                steps_text.insert("end", "\n=== UNBOUNDED PROBLEM ===\n")
            return None, None, None, None, problem_history
        leaving = np.argmin(ratios)

        # Pivot
        pivot = tableau[leaving, entering]
        tableau[leaving, :] /= pivot
        for i in range(m + 1):
            if i != leaving:
                tableau[i, :] -= tableau[i, entering] * tableau[leaving, :]

        iter_count += 1
        problem_history.append((f"Iteration {iter_count}", tableau.copy()))
        if steps_text is not None:
            steps_text.insert("end", f"\n=== ITERATION {iter_count} ===\n")
            steps_text.insert("end", utils.format_tableau(tableau, n, m))

    # Extract solution
    solution = np.zeros(n)
    for j in range(n):
        col = tableau[:-1, j]
        if np.count_nonzero(np.abs(col - 1) < 1e-8) == 1 and np.count_nonzero(np.abs(col) < 1e-8) == m - 1:
            row = np.where(np.abs(col - 1) < 1e-8)[0][0]
            solution[j] = tableau[row, -1]
    slack_values = tableau[:-1, -1][n:]
    optimal_value = tableau[-1, -1]
    if obj_type == "max":
        optimal_value = -optimal_value

    # Duals: last row, coefficients of the slack variables
    dual_solution = tableau[-1, n:n + m]
    if original_constraint_count is not None and m > original_constraint_count:
        dual_solution = dual_solution[:original_constraint_count]

    return solution, slack_values, optimal_value, dual_solution, problem_history