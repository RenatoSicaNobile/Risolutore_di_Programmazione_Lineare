import numpy as np

def pretty_var(name, idx):
    sub = "₀₁₂₃₄₅₆₇₈₉"
    return f"{name}{''.join(sub[int(d)] for d in str(idx+1))}"

def canonical_problem(obj_type, obj_coeffs, constraints, nn_signs, nn_values):
    lines = []
    n = len(obj_coeffs)
    c = np.zeros(n)
    obj_terms = []
    for i, (sign, coeff) in enumerate(obj_coeffs):
        coeff = float(coeff)
        real_coeff = coeff if sign == "+" else -coeff
        c[i] = real_coeff
        if abs(real_coeff) == 1:
            term = f"{('-' if real_coeff<0 else '')}{pretty_var('x',i)}"
        else:
            term = f"{('-' if real_coeff<0 else '')}{abs(real_coeff):g}{pretty_var('x',i)}"
        obj_terms.append(term)
    obj_line = f"max z = {' + '.join(obj_terms).replace('+ -','- ')}" if obj_type=="max" else f"min z = {' + '.join(obj_terms).replace('+ -','- ')}"
    lines.append("Objective:")
    lines.append(obj_line)
    lines.append("")
    constr_lines = []
    A = []
    b = []
    ineqs = []
    for cons in constraints:
        lhs = []
        row = []
        for i, (sign, val) in enumerate(zip(cons["signs"], cons["values"])):
            coeff = float(val)
            real_coeff = coeff if sign == "+" else -coeff
            row.append(real_coeff)
            if abs(real_coeff) == 1:
                term = f"{('-' if real_coeff<0 else '')}{pretty_var('x',i)}"
            else:
                term = f"{('-' if real_coeff<0 else '')}{abs(real_coeff):g}{pretty_var('x',i)}"
            lhs.append(term)
        ineq = cons["ineq"]
        rhs = float(cons["rhs"])
        A.append(row)
        b.append(rhs)
        ineqs.append(ineq)
        constr_lines.append(f"{' + '.join(lhs).replace('+ -','- '):<22} {ineq} {rhs:g}")
    for i, (sign, val) in enumerate(zip(nn_signs, nn_values)):
        if sign in (">=", "≥"):
            constr_lines.append(f"{pretty_var('x',i)} ≥ 0")
        elif sign in ("<=", "≤"):
            constr_lines.append(f"{pretty_var('x',i)} ≤ 0")
    lines.append("Constraints:")
    for cl in constr_lines:
        lines.append(cl)
    lines.append("")
    return c, np.array(A), np.array(b), ineqs, "\n".join(lines)

def add_slacks(A, ineqs):
    m, n = A.shape
    n_slack = 0
    slack_names = []
    new_A = []
    for i, row in enumerate(A):
        new_row = list(row)
        slack_cols = [0]*m
        if ineqs[i] == "≤":
            slack_cols[i] = 1
            slack_names.append(f"S{i+1}")
        elif ineqs[i] == "≥":
            slack_cols[i] = -1
            slack_names.append(f"S{i+1}")
        elif ineqs[i] == "=":
            slack_cols[i] = 0
            slack_names.append(f"S{i+1}")
        new_row.extend(slack_cols)
        new_A.append(new_row)
    return np.array(new_A, dtype=float), slack_names

def build_tableau(A, b, c):
    m, n = A.shape
    tableau = np.zeros((m+1, n+1))
    tableau[:m, :n] = A
    tableau[:m, -1] = b
    tableau[-1, :n] = -c  # for maximization
    return tableau

def aligned_tableau(tableau, all_vars, basis=None):
    m, n = tableau.shape
    col_width = 8
    header = ["Base"] + [f"{v:^{col_width}}" for v in all_vars] + ["RHS"]
    lines = []
    lines.append(" ".join(header))
    if basis is None:
        basis = [len(all_vars) - m + 1 + i for i in range(m - 1)]
    for i in range(m - 1):
        row = [f"{all_vars[basis[i]]:<4}"] + [f"{tableau[i, j]:{col_width}.2f}" for j in range(n - 1)] + [f"{tableau[i, -1]:{col_width}.2f}"]
        lines.append(" ".join(row))
    zrow = ["Z   "] + [f"{tableau[-1, j]:{col_width}.2f}" for j in range(n - 1)] + [f"{tableau[-1, -1]:{col_width}.2f}"]
    lines.append(" ".join(zrow))
    return "\n".join(lines)

def simplex_steps(obj_type, obj_coeffs, constraints, nn_signs, nn_values):
    c, A, b, ineqs, cano_text = canonical_problem(obj_type, obj_coeffs, constraints, nn_signs, nn_values)
    steps = [cano_text]
    nvars = len(obj_coeffs)
    A_slack, slack_names = add_slacks(A, ineqs)
    all_vars = [pretty_var("x",i) for i in range(nvars)] + slack_names

    steps.append("Slack variable introduction:")
    for i, row in enumerate(A_slack):
        lhs = []
        for j, aij in enumerate(row):
            if abs(aij) == 1 and aij != 0:
                lhs.append(f"{('-' if aij<0 else '')}{all_vars[j]}")
            elif aij != 0:
                lhs.append(f"{('-' if aij<0 else '')}{abs(aij):g}{all_vars[j]}")
        steps.append(f"{' + '.join(lhs).replace('+ -','- '):<30} = {b[i]:g}")
    steps.append("")
    steps.append("Non-negativity: " + ", ".join([f"{v} ≥ 0" for v in all_vars]))
    steps.append("")

    # Pad c for slack variables
    c_full = np.concatenate([c, np.zeros(len(slack_names))])
    tableau = build_tableau(A_slack, b, c_full)
    steps.append("Initial simplex tableau:")
    steps.append(aligned_tableau(tableau, all_vars))
    steps.append("")

    m, n = tableau.shape
    basis = [nvars + i for i in range(len(slack_names))]
    max_iter = 100
    solution_found = False
    for it in range(max_iter):
        obj_row = tableau[-1, :-1]
        if np.all(obj_row >= -1e-10):
            steps.append("All Z-row coefficients are non-negative. Optimal solution found.\n")
            solution_found = True
            break
        entering = np.argmin(obj_row)
        ratios = []
        for i in range(m - 1):
            if tableau[i, entering] > 1e-10:
                ratios.append(tableau[i, -1] / tableau[i, entering])
            else:
                ratios.append(np.inf)
        if np.all(np.isinf(ratios)):
            steps.append("Problem is unbounded: no finite solution.")
            return steps, None, None, None, None
        leaving = np.argmin(ratios)
        steps.append(f"Pivot: Enter {all_vars[entering]}, Leave {all_vars[basis[leaving]]} (row {leaving+1})")
        tableau[leaving, :] /= tableau[leaving, entering]
        for i in range(m):
            if i != leaving:
                tableau[i, :] -= tableau[i, entering] * tableau[leaving, :]
        basis[leaving] = entering
        steps.append("New tableau:")
        steps.append(aligned_tableau(tableau, all_vars, basis))
        steps.append("")
    sol = np.zeros(len(all_vars))
    for i, var in enumerate(basis):
        sol[var] = tableau[i, -1]
    optval = tableau[-1, -1]
    steps.append("Optimal solution found:")
    steps.append("  " + ", ".join([f"{all_vars[i]} = {sol[i]:.4f}" for i in range(len(all_vars))]))
    steps.append(f"  Optimal value: {optval:.4f}")
    steps.append("")

    # Check for alternative optima: if any non-basic variable in Z-row has zero reduced cost
    zero_reduced = []
    for j in range(len(all_vars)):
        if j not in basis:
            reduced_cost = tableau[-1, j]
            if np.isclose(reduced_cost, 0):
                zero_reduced.append(j)
    if zero_reduced:
        altvars = ", ".join([all_vars[j] for j in zero_reduced])
        steps.append("Multiple optimal solutions detected!")
        steps.append(f"  Variables with zero reduced cost: {altvars}")
        steps.append("  All points on the edge(s) defined by these variables (with others fixed) are optimal.")
        # For 2 variables, explain the segment
        xopt = sol[:nvars]
        eqs = []
        for j in zero_reduced:
            if j < nvars:
                s = [f"{all_vars[i]} = {xopt[i]:.4f}" for i in range(nvars) if i != j]
                s.append(f"{all_vars[j]} in [{0:.4f}, {xopt[j]:.4f}]")
                eqs.append(", ".join(s))
        if eqs:
            steps.append("  Segment of optimal solutions (for variables):")
            for eq in eqs:
                steps.append("    " + eq)
        steps.append("")
    return steps, sol[:nvars], optval, tableau, all_vars

def simplex_full(obj_type, obj_coeffs, constraints, nn_signs, nn_values):
    steps, xsol, optval, tableau, varnames = simplex_steps(obj_type, obj_coeffs, constraints, nn_signs, nn_values)
    return steps, xsol, optval, tableau, varnames