import numpy as np
import logging
import traceback
import utils  # Absolute import
import tkinter as tk  # Import tkinter

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
        constraints: list of tuples, each representing a constraint (see GUI code).
        steps_text: Tkinter Text widget or None; outputs steps if provided.
    Returns:
        solution: optimal variable values
        slack_values: slack variable values
        optimal_value: optimal value of the objective function
        final_dual: dual variable values (shadow prices)
        problem_history: list of (description, tableau) tuples
    """
    try:
        m, n = A.shape

        # Record initial tableau
        if steps_text is not None:
            steps_text.insert(tk.END, "\n=== INITIAL TABLEAU ===\n")

        # Add slack variables
        tableau = np.zeros((m + 1, n + m + 1))
        tableau[:-1, :n] = A
        tableau[:-1, n:n + m] = np.eye(m)
        tableau[:-1, -1] = b
        tableau[-1, :n] = c  # We're always minimizing in standard form

        # Display initial tableau
        if steps_text is not None:
            steps_text.insert(tk.END, utils.format_tableau(tableau, n, m))
        problem_history = [("Initial Tableau", tableau.copy())]

        # Bland's rule for anti-cycling
        iteration = 1
        while True:
            # Check for optimality
            if np.all(tableau[-1, :-1] >= -1e-6):
                if steps_text is not None:
                    steps_text.insert(tk.END, "\n=== OPTIMALITY CONDITION SATISFIED ===\n")
                break

            # Select entering variable (smallest index with negative reduced cost)
            entering = np.argmin(tableau[-1, :-1])
            entering_var = f"x{entering + 1}" if entering < n else f"s{entering - n + 1}"

            if steps_text is not None:
                steps_text.insert(tk.END, f"\n=== ITERATION {iteration} ===\n")
                steps_text.insert(tk.END, f"Entering variable: {entering_var} (column {entering + 1})\n")

            if tableau[-1, entering] >= -1e-10:
                break  # Optimal solution found

            # Select leaving variable using minimum ratio test
            ratios = np.where(tableau[:-1, entering] > 1e-10,
                              tableau[:-1, -1] / tableau[:-1, entering], np.inf)

            if np.all(ratios == np.inf):
                if steps_text is not None:
                    steps_text.insert(tk.END, "\n=== UNBOUNDED PROBLEM ===\n")
                return None, None, None, None, problem_history  # Indicate unbounded

            leaving = np.argmin(ratios)
            leaving_var = f"s{leaving + 1}"

            if steps_text is not None:
                steps_text.insert(tk.END, f"Leaving variable: {leaving_var} (row {leaving + 1})\n")

            # Pivot
            pivot = tableau[leaving, entering]
            tableau[leaving] /= pivot

            for i in range(m + 1):
                if i != leaving:
                    tableau[i] -= tableau[i, entering] * tableau[leaving]

            # Update basis variables (for display)
            basis_vars = []
            for col in range(n):
                col_vals = tableau[:-1, col]
                if np.sum(np.isclose(col_vals, 1)) == 1 and np.sum(np.isclose(col_vals, 0)) == m - 1:
                    basis_vars.append(f"x{col + 1}")
            for col in range(m):
                col_vals = tableau[:-1, n + col]
                if np.sum(np.isclose(col_vals, 1)) == 1 and np.sum(np.isclose(col_vals, 0)) == m - 1:
                    basis_vars.append(f"s{col + 1}")

            if steps_text is not None:
                steps_text.insert(tk.END, f"Basis variables: {', '.join(basis_vars)}\n")
            problem_history.append((f"Iteration {iteration} - After pivot", tableau.copy()))

            iteration += 1

        # Extract solution
        solution = np.zeros(n)
        basis = []

        for col in range(n):
            col_vals = tableau[:-1, col]
            if np.sum(np.isclose(col_vals, 1)) == 1 and np.sum(np.isclose(col_vals, 0)) == m - 1:
                row = np.argmax(col_vals)
                solution[col] = tableau[row, -1]
                basis.append((row, col))

        slack_values = b - np.dot(A, solution)

        # Calculate optimal value
        optimal_value = tableau[-1, -1]
        if obj_type == "max":
            optimal_value = -optimal_value

        # Extract dual variables (shadow prices)
        dual_solution = np.zeros(m)
        for row in range(m):
            if n + row < tableau.shape[1] - 1:
                dual_solution[row] = tableau[-1, n + row]

        # Adjust dual variables for converted constraints
        final_dual = np.zeros(original_constraint_count)
        current_idx = 0
        for i, op in enumerate([x[3].get() if hasattr(x[3], "get") else x[3] for x in constraints]):
            if op == "≤":
                final_dual[i] = dual_solution[current_idx]
                current_idx += 1
            elif op == "≥":
                final_dual[i] = -dual_solution[current_idx]
                current_idx += 1
            elif op == "=":
                final_dual[i] = dual_solution[current_idx] - dual_solution[current_idx + 1]
                current_idx += 2

        return solution, slack_values, optimal_value, final_dual, problem_history

    except Exception as e:
        logging.error(f"Error in run_simplex: {e}")
        logging.error(traceback.format_exc())
        return None, None, None, None, []  # Indicate failure