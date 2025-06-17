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
    lines.append("Problema in forma canonica:")
    lines.append(obj_line)
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
        constr_lines.append(f"{' + '.join(lhs).replace('+ -','- ')} {ineq} {rhs:g}")
    for i, (sign, val) in enumerate(zip(nn_signs, nn_values)):
        if sign in (">=", "≥"):
            constr_lines.append(f"{pretty_var('x',i)} ≥ 0")
        elif sign in ("<=", "≤"):
            constr_lines.append(f"{pretty_var('x',i)} ≤ 0")
    lines.append("Vincoli:")
    for cl in constr_lines:
        lines.append(cl)
    return c, np.array(A), np.array(b), ineqs, "\n".join(lines)

def add_slacks(A, ineqs):
    # Adds slack variables to each constraint, returns new A, names
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

def simplex_steps(obj_type, obj_coeffs, constraints, nn_signs, nn_values):
    c, A, b, ineqs, cano_text = canonical_problem(obj_type, obj_coeffs, constraints, nn_signs, nn_values)
    steps = [cano_text, ""]
    nvars = len(obj_coeffs)
    A_slack, slack_names = add_slacks(A, ineqs)
    all_vars = [pretty_var("x",i) for i in range(nvars)] + slack_names

    steps.append("Introduzione delle variabili di slack:")
    for i, row in enumerate(A_slack):
        lhs = []
        for j, aij in enumerate(row):
            if abs(aij) == 1 and aij != 0:
                lhs.append(f"{('-' if aij<0 else '')}{all_vars[j]}")
            elif aij != 0:
                lhs.append(f"{('-' if aij<0 else '')}{abs(aij):g}{all_vars[j]}")
        steps.append(f"{' + '.join(lhs).replace('+ -','- ')} = {b[i]:g}")
    steps.append("")
    steps.append("Condizioni: " + ", ".join([f"{v} ≥ 0" for v in all_vars]))

    # Pad c with zeros for slack variables
    c_full = np.concatenate([c, np.zeros(len(slack_names))])

    tableau = build_tableau(A_slack, b, c_full)
    steps.append("\nTableau iniziale del simplesso:")
    steps.append(format_tableau(tableau, all_vars))

    m, n = tableau.shape
    basis = [nvars + i for i in range(len(slack_names))]
    max_iter = 100
    for it in range(max_iter):
        obj_row = tableau[-1, :-1]
        if np.all(obj_row >= -1e-10):
            steps.append("\nTutti i coefficienti nella riga z sono non negativi. Soluzione ottima trovata.")
            break
        entering = np.argmin(obj_row)
        ratios = []
        for i in range(m - 1):
            if tableau[i, entering] > 1e-10:
                ratios.append(tableau[i, -1] / tableau[i, entering])
            else:
                ratios.append(np.inf)
        if np.all(np.isinf(ratios)):
            steps.append("Il problema è illimitato: nessuna soluzione finita.")
            return steps, None
        leaving = np.argmin(ratios)
        steps.append(f"\nPivot: entra {all_vars[entering]}, esce {all_vars[basis[leaving]]} (riga {leaving+1})")
        tableau[leaving, :] /= tableau[leaving, entering]
        for i in range(m):
            if i != leaving:
                tableau[i, :] -= tableau[i, entering]*tableau[leaving, :]
        basis[leaving] = entering
        steps.append("Nuovo tableau:")
        steps.append(format_tableau(tableau, all_vars, basis))
    sol = np.zeros(len(all_vars))
    for i, var in enumerate(basis):
        sol[var] = tableau[i, -1]
    optval = tableau[-1, -1]
    steps.append(f"\nSoluzione ottima trovata: " + ", ".join([f"{all_vars[i]} = {sol[i]:.4f}" for i in range(len(all_vars))]))
    steps.append(f"Valore ottimo: {optval:.4f}")
    xsol = sol[:nvars]
    return steps, (xsol, optval)

def format_tableau(tableau, all_vars, basis=None):
    m, n = tableau.shape
    lines = []
    header = ["Base"] + [f"{v}" for v in all_vars] + ["RHS"]
    lines.append(" | ".join(header))
    if basis is None:
        basis = [len(all_vars) - m + 1 + i for i in range(m - 1)]
    for i in range(m - 1):
        row = [f"{all_vars[basis[i]]}"] + [f"{tableau[i, j]:7.2f}" for j in range(n - 1)] + [f"{tableau[i, -1]:7.2f}"]
        lines.append(" | ".join(row))
    zrow = ["Z"] + [f"{tableau[-1, j]:7.2f}" for j in range(n - 1)] + [f"{tableau[-1, -1]:7.2f}"]
    lines.append(" | ".join(zrow))
    return "\n".join(lines)

def simplex_full(obj_type, obj_coeffs, constraints, nn_signs, nn_values):
    steps, result = simplex_steps(obj_type, obj_coeffs, constraints, nn_signs, nn_values)
    if result is not None:
        xsol, optval = result
        return steps, xsol, optval
    else:
        return steps, None, None