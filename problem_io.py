import json
import tkinter as tk
from tkinter import filedialog, messagebox

def save_problem(variables, obj_type, obj_coeffs, obj_signs, constraints, nn_signs, nn_values, integer_vars):
    """Save the current problem to a file."""
    filename = filedialog.asksaveasfilename(
        defaultextension=".json",
        filetypes=[("JSON Files", "*.json")],
        title="Salva problema"
    )

    if not filename:
        return

    try:
        problem_data = {
            "variables": variables,
            "obj_type": obj_type.get(),
            "obj_coeffs": [entry.get() for entry in obj_coeffs],
            "obj_signs": [sign.get() for sign in obj_signs],
            "constraints": [
                {
                    "signs": [sign.get() for sign in constraint[1]],
                    "coeffs": [entry.get() for entry in constraint[2]],
                    "ineq": constraint[3].get(),
                    "rhs": constraint[4].get()
                }
                for constraint in constraints
            ],
            "nn_signs": [sign.get() for sign in nn_signs],
            "nn_values": [val.get() for val in nn_values],
            "integer_vars": [var.get() for var in integer_vars]
        }

        with open(filename, "w") as f:
            json.dump(problem_data, f, indent=4)

        messagebox.showinfo("Successo", f"Problema salvato con successo:\n{filename}")
        return True  # Indicate success
    except Exception as e:
        messagebox.showerror("Errore", f"Impossibile salvare il problema:\n{str(e)}")
        return False  # Indicate failure

def load_problem(obj_type, obj_coeffs, obj_signs, constraints, nn_signs, nn_values, integer_vars):
    """Load a problem from a file."""
    filename = filedialog.askopenfilename(
        filetypes=[("JSON Files", "*.json")],
        title="Carica problema"
    )

    if not filename:
        return None  # Indicate no file selected

    try:
        with open(filename, "r") as f:
            problem_data = json.load(f)

        return problem_data  # Return the loaded data
    except Exception as e:
        messagebox.showerror("Errore", f"Impossibile caricare il problema:\n{str(e)}")
        return None  # Indicate failure