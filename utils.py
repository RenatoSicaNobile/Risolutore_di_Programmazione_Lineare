import numpy as np
import tkinter as tk

def format_tableau(tableau, n, m):
    """Formats the tableau as a string for display."""
    s = ""
    for i in range(tableau.shape[0]):
        for j in range(tableau.shape[1]):
            s += "{:8.4f} ".format(tableau[i, j])
        s += "\n"
    return s

def is_number(s):
    """Checks if a string is a number (int or float)."""
    try:
        float(s)
        return True
    except ValueError:
        return False