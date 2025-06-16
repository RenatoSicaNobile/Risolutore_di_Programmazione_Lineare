import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import json
from PIL import Image, ImageTk
import graph

def subscript_number(n):
    subscript_digits = "₀₁₂₃₄₅₆₇₈₉"
    return ''.join(subscript_digits[int(d)] for d in str(n))

class LinearProgrammingGUI:
    def __init__(self, master):
        self.master = master
        master.title("Risolutore di Programmazione Lineare")
        master.geometry("1200x800")

        # Logo and Title
        logo_title_frame = tk.Frame(master)
        logo_title_frame.pack(side=tk.TOP, fill=tk.X, pady=5)
        try:
            self.logo_img = Image.open("rock.jpg")
            self.logo_img = self.logo_img.resize((90, 90), Image.LANCZOS)
            self.logo_photo = ImageTk.PhotoImage(self.logo_img)
            logo_label = tk.Label(logo_title_frame, image=self.logo_photo)
            logo_label.pack(side=tk.LEFT, padx=12)
        except Exception:
            logo_label = tk.Label(logo_title_frame, text="LOGO")
            logo_label.pack(side=tk.LEFT, padx=12)

        title_label = tk.Label(
            logo_title_frame,
            text="RISOLUTORE DI PROGRAMMAZIONE LINEARE",
            font=("Arial", 22, "bold"),
            anchor="center",
            justify="center"
        )
        title_label.pack(side=tk.LEFT, expand=True, padx=20)

        main_frame = tk.Frame(master)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # LEFT: Problem definition
        self.problem_frame = ttk.LabelFrame(main_frame, text="Definizione del Problema")
        self.problem_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10, anchor="n")

        # Objective Function Type
        self.obj_type_label = ttk.Label(self.problem_frame, text="Tipo funzione obiettivo:")
        self.obj_type_label.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.obj_type = tk.StringVar(value="max")
        self.obj_type_combo = ttk.Combobox(
            self.problem_frame, textvariable=self.obj_type, values=["max", "min"], state="readonly"
        )
        self.obj_type_combo.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W)
        self.obj_type_combo.bind('<<ComboboxSelected>>', lambda e: self.show_objective_function())

        # Number of Variables
        self.num_vars_label = ttk.Label(self.problem_frame, text="Numero di variabili:")
        self.num_vars_label.grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.num_vars = tk.IntVar(value=2)
        self.num_vars_entry = ttk.Entry(self.problem_frame, textvariable=self.num_vars)
        self.num_vars_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W)
        self.add_var_button = ttk.Button(
            self.problem_frame, text="AGGIORNA VARIABILI", command=self.update_variables
        )
        self.add_var_button.grid(row=1, column=2, padx=5, pady=5)

        # Objective Function Display
        self.obj_func_label = ttk.Label(self.problem_frame, text="Funzione obiettivo:")
        self.obj_func_label.grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.obj_func_value = ttk.Label(
            self.problem_frame, text="", font=("Arial", 12, "bold")
        )
        self.obj_func_value.grid(row=2, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)

        self.obj_coeffs = []
        self.obj_signs = []
        self.obj_coeffs_frame = ttk.Frame(self.problem_frame)
        self.obj_coeffs_frame.grid(row=3, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)

        # Constraints
        self.constraints_label = ttk.Label(self.problem_frame, text="Vincoli:")
        self.constraints_label.grid(row=4, column=0, padx=5, pady=5, sticky=tk.W)
        self.constraints_frame = ttk.Frame(self.problem_frame)
        self.constraints_frame.grid(row=4, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)
        self.constraints = []
        self.add_constraint_button = ttk.Button(
            self.problem_frame, text="AGGIUNGI VINCOLO", command=self.add_constraint
        )
        self.add_constraint_button.grid(row=5, column=0, padx=5, pady=5)
        self.remove_constraint_button = ttk.Button(
            self.problem_frame, text="RIMUOVI VINCOLO", command=self.remove_constraint
        )
        self.remove_constraint_button.grid(row=5, column=1, padx=5, pady=5)

        # Non-Negativity Constraints
        self.nn_constraints_label = ttk.Label(self.problem_frame, text="Vincoli di non negatività:")
        self.nn_constraints_label.grid(row=6, column=0, padx=5, pady=5, sticky=tk.W)
        self.nn_constraints_frame = ttk.Frame(self.problem_frame)
        self.nn_constraints_frame.grid(row=6, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)
        self.nn_signs = []
        self.nn_values = []

        # Integer Variables
        self.integer_vars_label = ttk.Label(self.problem_frame, text="Variabili intere:")
        self.integer_vars_label.grid(row=7, column=0, padx=5, pady=5, sticky=tk.W)
        self.integer_vars = []
        self.integer_vars_frame = ttk.Frame(self.problem_frame)
        self.integer_vars_frame.grid(row=7, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)

        self.reset_button = ttk.Button(self.problem_frame, text="RESETTA", command=self.reset_fields)
        self.reset_button.grid(row=8, column=0, padx=5, pady=5)
        self.load_button = ttk.Button(self.problem_frame, text="CARICA PROBLEMA", command=self.load_problem)
        self.load_button.grid(row=8, column=1, padx=5, pady=5)
        self.save_button = ttk.Button(self.problem_frame, text="SALVA PROBLEMA", command=self.save_problem)
        self.save_button.grid(row=8, column=2, padx=5, pady=5)

        # RIGHT: Solution description/process only
        right_frame = tk.Frame(main_frame)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.solution_frame = ttk.LabelFrame(right_frame, text="Processo risolutivo")
        self.solution_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))
        self.solution_text = tk.Text(self.solution_frame, height=30, wrap=tk.WORD)
        self.solution_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Button row
        button_row = tk.Frame(self.solution_frame)
        button_row.pack(fill=tk.X, padx=5, pady=2)
        self.solve_button = ttk.Button(button_row, text="RISOLVI", command=self.solve)
        self.solve_button.pack(side=tk.LEFT, padx=5)
        self.export_button = ttk.Button(button_row, text="ESPORTA IN PDF", command=self.export_pdf)
        self.export_button.pack(side=tk.LEFT, padx=5)
        self.copy_button = ttk.Button(button_row, text="COPIA", command=self.copy_solution)
        self.copy_button.pack(side=tk.LEFT, padx=5)

        self.show_all_quadrants = tk.BooleanVar(value=False)
        self.show_all_quadrants_check = ttk.Checkbutton(
            button_row, text="MOSTRA QUADRANTI", variable=self.show_all_quadrants, command=self.plot_graph
        )
        self.show_all_quadrants_check.pack(side=tk.LEFT, padx=5)
        self.graph_button = ttk.Button(button_row, text="GRAFICO", command=self.open_graph_window)
        self.graph_button.pack(side=tk.LEFT, padx=5)

        # Store last solution for graphing
        self.last_solution = None
        self.last_A = None
        self.last_b = None
        self.last_c = None

        self.update_variables()

    def copy_solution(self):
        try:
            self.master.clipboard_clear()
            self.master.clipboard_append(self.solution_text.get("1.0", tk.END))
            messagebox.showinfo("Copiato", "Soluzione copiata negli appunti!")
        except Exception:
            messagebox.showerror("Errore", "Impossibile copiare la soluzione negli appunti.")

    def update_variables(self):
        num_vars = self.num_vars.get()
        for widget in self.obj_coeffs_frame.winfo_children():
            widget.destroy()
        self.obj_coeffs = []
        self.obj_signs = []
        for i in range(num_vars):
            sign_var = tk.StringVar(value="+")
            sign_combo = ttk.Combobox(
                self.obj_coeffs_frame, textvariable=sign_var, values=["+", "-"], width=3, state="readonly"
            )
            sign_combo.grid(row=0, column=i * 3, padx=2, pady=2)
            self.obj_signs.append(sign_var)
            sign_var.trace_add('write', lambda *args: self.show_objective_function())
            coeff_entry = ttk.Entry(self.obj_coeffs_frame, width=5)
            coeff_entry.insert(0, "1.0")
            coeff_entry.grid(row=0, column=i * 3 + 1, padx=2, pady=2)
            coeff_entry.bind('<KeyRelease>', lambda e: self.show_objective_function())
            self.obj_coeffs.append(coeff_entry)
            var_label = ttk.Label(self.obj_coeffs_frame, text="x" + subscript_number(i + 1))
            var_label.grid(row=0, column=i * 3 + 2, padx=2, pady=2)
        self.show_objective_function()

        for widget in self.nn_constraints_frame.winfo_children():
            widget.destroy()
        self.nn_signs = []
        self.nn_values = []
        for i in range(num_vars):
            var_label = ttk.Label(self.nn_constraints_frame, text="x" + subscript_number(i + 1) + " ≥")
            var_label.grid(row=0, column=i * 2, padx=2, pady=2)
            self.nn_signs.append(">=")
            value_entry = ttk.Entry(self.nn_constraints_frame, width=5)
            value_entry.insert(0, "0")
            value_entry.grid(row=0, column=i * 2 + 1, padx=2, pady=2)
            self.nn_values.append(value_entry)

        for widget in self.integer_vars_frame.winfo_children():
            widget.destroy()
        self.integer_vars = []
        for i in range(num_vars):
            var_check = ttk.Checkbutton(self.integer_vars_frame, text="x" + subscript_number(i + 1))
            var_check.grid(row=0, column=i, padx=2, pady=2)
            self.integer_vars.append(var_check)

        # Update constraints entries to match variable count
        for idx, (frame, signs, entries, ineq, rhs_entry) in enumerate(self.constraints):
            for widget in frame.winfo_children():
                widget.destroy()
            signs.clear()
            entries.clear()
            for j in range(num_vars):
                sign_var = tk.StringVar(value="+")
                sign_combo = ttk.Combobox(
                    frame, textvariable=sign_var, values=["+", "-"], width=3, state="readonly"
                )
                sign_combo.grid(row=0, column=j * 3, padx=2, pady=2)
                signs.append(sign_var)
                entry = ttk.Entry(frame, width=5)
                entry.insert(0, "1.0")
                entry.grid(row=0, column=j * 3 + 1, padx=2, pady=2)
                entries.append(entry)
                var_label = ttk.Label(frame, text="x" + subscript_number(j + 1))
                var_label.grid(row=0, column=j * 3 + 2, padx=2, pady=2)
            ineq = tk.StringVar(value="≤")
            ineq_combo = ttk.Combobox(
                frame, textvariable=ineq, values=["≤", "≥", "="], width=3, state="readonly"
            )
            ineq_combo.grid(row=0, column=num_vars * 3, padx=2, pady=2)
            new_rhs_entry = ttk.Entry(frame, width=5)
            new_rhs_entry.insert(0, "0.0")
            new_rhs_entry.grid(row=0, column=num_vars * 3 + 1, padx=2, pady=2)
            self.constraints[idx] = (frame, signs, entries, ineq, new_rhs_entry)

    def show_objective_function(self):
        num_vars = self.num_vars.get()
        terms = []
        for i in range(num_vars):
            sign = self.obj_signs[i].get() if i < len(self.obj_signs) else "+"
            coeff = self.obj_coeffs[i].get() if i < len(self.obj_coeffs) else "1"
            var = "x" + subscript_number(i + 1)
            try:
                coeff_val = float(coeff)
            except Exception:
                coeff_val = 0
            if coeff_val == 0:
                continue
            # Fix: Avoid duplicate operators
            sign_str = "+" if sign == "+" else "-"
            # Always show sign except for first positive term
            if i == 0:
                if sign == "-":
                    term = f"-{abs(coeff_val):g}{var}" if abs(coeff_val) != 1 else f"-{var}"
                else:
                    term = f"{abs(coeff_val):g}{var}" if abs(coeff_val) != 1 else f"{var}"
            else:
                if sign == "-":
                    term = f" - {abs(coeff_val):g}{var}" if abs(coeff_val) != 1 else f" - {var}"
                else:
                    term = f" + {abs(coeff_val):g}{var}" if abs(coeff_val) != 1 else f" + {var}"
            terms.append(term)
        if not terms:
            terms = ["0"]
        obj_type = self.obj_type.get()
        obj_str = ("max" if obj_type == "max" else "min") + " z = " + "".join(terms)
        self.obj_func_value.config(text=obj_str)

    def add_constraint(self):
        num_vars = self.num_vars.get()
        frame = ttk.Frame(self.constraints_frame)
        frame.pack(pady=2, padx=2)
        signs = []
        entries = []
        for i in range(num_vars):
            sign_var = tk.StringVar(value="+")
            sign_combo = ttk.Combobox(
                frame, textvariable=sign_var, values=["+", "-"], width=3, state="readonly"
            )
            sign_combo.grid(row=0, column=i * 3, padx=2, pady=2)
            signs.append(sign_var)
            entry = ttk.Entry(frame, width=5)
            entry.insert(0, "1.0")
            entry.grid(row=0, column=i * 3 + 1, padx=2, pady=2)
            entries.append(entry)
            var_label = ttk.Label(frame, text="x" + subscript_number(i + 1))
            var_label.grid(row=0, column=i * 3 + 2, padx=2, pady=2)
        ineq = tk.StringVar(value="≤")
        ineq_combo = ttk.Combobox(
            frame, textvariable=ineq, values=["≤", "≥", "="], width=3, state="readonly"
        )
        ineq_combo.grid(row=0, column=num_vars * 3, padx=2, pady=2)
        rhs_entry = ttk.Entry(frame, width=5)
        rhs_entry.insert(0, "0.0")
        rhs_entry.grid(row=0, column=num_vars * 3 + 1, padx=2, pady=2)
        constraint_tuple = (frame, signs, entries, ineq, rhs_entry)
        self.constraints.append(constraint_tuple)

    def remove_constraint(self):
        if self.constraints:
            frame, *_ = self.constraints[-1]
            frame.destroy()
            self.constraints.pop()

    def reset_fields(self):
        self.num_vars.set(2)
        self.obj_type.set("max")
        self.constraints.clear()
        for widget in self.constraints_frame.winfo_children():
            widget.destroy()
        self.update_variables()
        self.solution_text.delete(1.0, tk.END)
        self.last_solution = None
        self.last_A = None
        self.last_b = None
        self.last_c = None

    def save_problem(self):
        data = self._get_problem_data()
        filename = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("File JSON", "*.json")],
            title="Salva problema"
        )
        if filename:
            with open(filename, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
            messagebox.showinfo("Salva problema", "Problema salvato con successo!")

    def load_problem(self):
        filename = filedialog.askopenfilename(
            filetypes=[("File JSON", "*.json")],
            title="Carica problema"
        )
        if filename:
            with open(filename, "r", encoding="utf-8") as f:
                data = json.load(f)
            self._set_problem_data(data)
            messagebox.showinfo("Carica problema", "Problema caricato con successo!")

    def _get_problem_data(self):
        num_vars = self.num_vars.get()
        obj_type = self.obj_type.get()
        obj_coeffs = [entry.get() for entry in self.obj_coeffs]
        obj_signs = [sign.get() for sign in self.obj_signs]
        nn_values = [entry.get() for entry in self.nn_values]
        integer_vars = [var.instate(['selected']) for var in self.integer_vars]
        constraints = []
        for frame, signs, entries, ineq, rhs_entry in self.constraints:
            constraint = {
                "signs": [s.get() for s in signs],
                "values": [e.get() for e in entries],
                "ineq": ineq.get(),
                "rhs": rhs_entry.get()
            }
            constraints.append(constraint)
        return {
            "num_vars": num_vars,
            "obj_type": obj_type,
            "obj_coeffs": obj_coeffs,
            "obj_signs": obj_signs,
            "nn_values": nn_values,
            "integer_vars": integer_vars,
            "constraints": constraints
        }

    def _set_problem_data(self, data):
        self.num_vars.set(data.get("num_vars", 2))
        self.obj_type.set(data.get("obj_type", "max"))
        self.update_variables()
        for i, (c, s) in enumerate(zip(data.get("obj_coeffs", []), data.get("obj_signs", []))):
            if i < len(self.obj_coeffs):
                self.obj_coeffs[i].delete(0, tk.END)
                self.obj_coeffs[i].insert(0, c)
            if i < len(self.obj_signs):
                self.obj_signs[i].set(s)
        for i, v in enumerate(data.get("nn_values", [])):
            if i < len(self.nn_values):
                self.nn_values[i].delete(0, tk.END)
                self.nn_values[i].insert(0, v)
        for i, selected in enumerate(data.get("integer_vars", [])):
            if i < len(self.integer_vars):
                if selected:
                    self.integer_vars[i].state(['selected'])
                else:
                    self.integer_vars[i].state(['!selected'])
        for widget in self.constraints_frame.winfo_children():
            widget.destroy()
        self.constraints = []
        for cons in data.get("constraints", []):
            self.add_constraint()
            frame, signs, entries, ineq, rhs_entry = self.constraints[-1]
            for j, s in enumerate(cons.get("signs", [])):
                if j < len(signs):
                    signs[j].set(s)
            for j, val in enumerate(cons.get("values", [])):
                if j < len(entries):
                    entries[j].delete(0, tk.END)
                    entries[j].insert(0, val)
            ineq.set(cons.get("ineq", "≤"))
            rhs_entry.delete(0, tk.END)
            rhs_entry.insert(0, cons.get("rhs", "0.0"))
        self.show_objective_function()

    def solve(self):
        try:
            num_vars = self.num_vars.get()
            obj_type = self.obj_type.get()
            obj_coeffs = []
            for i in range(num_vars):
                sign = self.obj_signs[i].get()
                coeff = float(self.obj_coeffs[i].get())
                obj_coeffs.append(coeff if sign == "+" else -coeff)
            A = []
            b = []
            ineq_ops = []
            for frame, signs, entries, ineq, rhs_entry in self.constraints:
                row = []
                for i in range(num_vars):
                    sign = signs[i].get()
                    coeff = float(entries[i].get())
                    row.append(coeff if sign == "+" else -coeff)
                A.append(row)
                b.append(float(rhs_entry.get()))
                ineq_ops.append(ineq.get())
            # Non-negativity constraints
            for i, entry in enumerate(self.nn_values):
                val = float(entry.get())
                row = [0.0] * num_vars
                if self.nn_signs[i] == ">=" or self.nn_signs[i] == "≥":
                    row[i] = -1.0
                    A.append(row)
                    b.append(-val)
                    ineq_ops.append("≤")
                elif self.nn_signs[i] == "<=" or self.nn_signs[i] == "≤":
                    row[i] = 1.0
                    A.append(row)
                    b.append(val)
                    ineq_ops.append("≤")
            # Standard form for simplex and plotting
            A_std = []
            b_std = []
            for ai, bi, op in zip(A, b, ineq_ops):
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
            A = np.array(A_std)
            b = np.array(b_std)
            c = -np.array(obj_coeffs) if obj_type == "max" else np.array(obj_coeffs)
            self.solution_text.delete("1.0", tk.END)
            sol, optimal = self.simplex_solve(c, A, b)
            self.last_solution = sol
            self.last_A = A
            self.last_b = b
            self.last_c = c
            if sol is not None:
                self.solution_text.insert(tk.END, "=== SOLUZIONE ===\n\n")
                self.solution_text.insert(tk.END, f"Valore ottimo ({obj_type}): z = {optimal:.4f}\n\n")
                self.solution_text.insert(tk.END, "Variabili decisionali:\n")
                for i, val in enumerate(sol):
                    self.solution_text.insert(tk.END, f"  x{subscript_number(i + 1)} = {val:.4f}\n")
                self.solution_text.insert(tk.END, "\n=== FINE SOLUZIONE ===\n")
            else:
                self.solution_text.insert(tk.END, "Il problema non ha una soluzione ammissibile.\n")
            self.plot_graph()
        except Exception as e:
            messagebox.showerror("Errore", f"Si è verificato un errore: {e}")

    def simplex_solve(self, c, A, b):
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
        if self.obj_type.get() == "max":
            optimal_value = -optimal_value
        return solution, optimal_value

    def plot_graph(self):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        if hasattr(self, "canvas") and self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
        fig, ax = plt.subplots(figsize=(7, 7))
        variables = [f"x{subscript_number(i + 1)}" for i in range(self.num_vars.get())]
        if self.last_c is not None and self.last_A is not None and self.last_b is not None:
            graph.plot_solution_to_axes(
                ax, self.last_c, self.last_A, self.last_b, variables,
                show_all_quadrants=self.show_all_quadrants.get(),
                solution=self.last_solution, show_solution=True
            )
        self.canvas = FigureCanvasTkAgg(fig, self.solution_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.draw()

    def open_graph_window(self):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        popup = tk.Toplevel(self.master)
        popup.title("Grafico - Regione Ammissibile")
        popup.geometry("820x820")
        fig, ax = plt.subplots(figsize=(8, 8))
        variables = [f"x{subscript_number(i + 1)}" for i in range(self.num_vars.get())]
        if self.last_c is not None and self.last_A is not None and self.last_b is not None:
            graph.plot_solution_to_axes(
                ax, self.last_c, self.last_A, self.last_b, variables,
                show_all_quadrants=self.show_all_quadrants.get(),
                solution=self.last_solution, show_solution=True
            )
        canvas = FigureCanvasTkAgg(fig, master=popup)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, popup)
        toolbar.update()

    def export_pdf(self):
        # Implement PDF export logic here if needed
        pass

if __name__ == "__main__":
    root = tk.Tk()
    gui = LinearProgrammingGUI(root)
    root.mainloop()