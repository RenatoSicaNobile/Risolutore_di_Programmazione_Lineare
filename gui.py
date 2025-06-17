import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import json
import numpy as np
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import io

import lp_solver
import graph
import pdf_exporter

def pedice_numero(n):
    pedice = "₀₁₂₃₄₅₆₇₈₉"
    return ''.join(pedice[int(d)] for d in str(n))

class RisolutorePLGUI:
    def __init__(self, master):
        self.master = master
        master.title("Risolutore di Programmazione Lineare")
        master.geometry("1200x800")

        # Logo e Titolo
        frame_logo_titolo = tk.Frame(master)
        frame_logo_titolo.pack(side=tk.TOP, fill=tk.X, pady=5)
        try:
            self.logo_img = Image.open("rock.jpg")
            self.logo_img = self.logo_img.resize((90, 90), Image.LANCZOS)
            self.logo_photo = ImageTk.PhotoImage(self.logo_img)
            label_logo = tk.Label(frame_logo_titolo, image=self.logo_photo)
            label_logo.pack(side=tk.LEFT, padx=12)
        except Exception:
            label_logo = tk.Label(frame_logo_titolo, text="LOGO")
            label_logo.pack(side=tk.LEFT, padx=12)

        label_titolo = tk.Label(
            frame_logo_titolo,
            text="RISOLUTORE DI PROGRAMMAZIONE LINEARE",
            font=("Arial", 22, "bold"),
            anchor="center",
            justify="center"
        )
        label_titolo.pack(side=tk.LEFT, expand=True, padx=20)

        frame_principale = tk.Frame(master)
        frame_principale.pack(fill=tk.BOTH, expand=True)

        # SINISTRA: Definizione problema
        self.frame_problema = ttk.LabelFrame(frame_principale, text="Definizione del Problema")
        self.frame_problema.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10, anchor="n")

        # Tipo funzione obiettivo
        self.label_tipo_fo = ttk.Label(self.frame_problema, text="Tipo funzione obiettivo:")
        self.label_tipo_fo.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.tipo_fo = tk.StringVar(value="max")
        self.combo_tipo_fo = ttk.Combobox(
            self.frame_problema, textvariable=self.tipo_fo, values=["max", "min"], state="readonly"
        )
        self.combo_tipo_fo.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W)
        self.combo_tipo_fo.bind('<<ComboboxSelected>>', lambda e: self.mostra_fo())

        # Numero variabili
        self.label_num_var = ttk.Label(self.frame_problema, text="Numero di variabili:")
        self.label_num_var.grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.num_var = tk.IntVar(value=2)
        self.entry_num_var = ttk.Entry(self.frame_problema, textvariable=self.num_var)
        self.entry_num_var.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W)
        self.bottone_aggiorna_var = ttk.Button(
            self.frame_problema, text="AGGIORNA VARIABILI", command=self.aggiorna_variabili
        )
        self.bottone_aggiorna_var.grid(row=1, column=2, padx=5, pady=5)

        # Funzione obiettivo
        self.label_fo = ttk.Label(self.frame_problema, text="Funzione obiettivo:")
        self.label_fo.grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.label_fo_valore = ttk.Label(
            self.frame_problema, text="", font=("Arial", 12, "bold")
        )
        self.label_fo_valore.grid(row=2, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)

        self.fo_coeff = []
        self.fo_segn = []
        self.frame_fo_coeff = ttk.Frame(self.frame_problema)
        self.frame_fo_coeff.grid(row=3, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)

        # Vincoli
        self.label_vincoli = ttk.Label(self.frame_problema, text="Vincoli:")
        self.label_vincoli.grid(row=4, column=0, padx=5, pady=5, sticky=tk.W)
        self.frame_vincoli = ttk.Frame(self.frame_problema)
        self.frame_vincoli.grid(row=4, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)
        self.vincoli = []
        self.bottone_aggiungi_vincolo = ttk.Button(
            self.frame_problema, text="AGGIUNGI VINCOLO", command=self.aggiungi_vincolo
        )
        self.bottone_aggiungi_vincolo.grid(row=5, column=0, padx=5, pady=5)
        self.bottone_rimuovi_vincolo = ttk.Button(
            self.frame_problema, text="RIMUOVI VINCOLO", command=self.rimuovi_vincolo
        )
        self.bottone_rimuovi_vincolo.grid(row=5, column=1, padx=5, pady=5)

        # Vincoli di non negatività
        self.label_nn = ttk.Label(self.frame_problema, text="Vincoli di non negatività:")
        self.label_nn.grid(row=6, column=0, padx=5, pady=5, sticky=tk.W)
        self.frame_nn = ttk.Frame(self.frame_problema)
        self.frame_nn.grid(row=6, column=1, columnspan=2, padx=5, pady=5, sticky=tk.W)
        self.nn_segn = []
        self.nn_val = []

        # Pulsanti gestione
        self.bottone_reset = ttk.Button(self.frame_problema, text="RESETTA", command=self.resetta)
        self.bottone_reset.grid(row=7, column=0, padx=5, pady=5)
        self.bottone_carica = ttk.Button(self.frame_problema, text="CARICA PROBLEMA", command=self.carica_problema)
        self.bottone_carica.grid(row=7, column=1, padx=5, pady=5)
        self.bottone_salva = ttk.Button(self.frame_problema, text="SALVA PROBLEMA", command=self.salva_problema)
        self.bottone_salva.grid(row=7, column=2, padx=5, pady=5)

        # DESTRA: Output soluzione
        frame_destra = tk.Frame(frame_principale)
        frame_destra.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.frame_soluzione = ttk.LabelFrame(frame_destra, text="Procedura risolutiva")
        self.frame_soluzione.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))
        self.text_soluzione = tk.Text(self.frame_soluzione, height=30, wrap=tk.WORD, font=("Consolas", 12))
        self.text_soluzione.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        frame_bottoni = tk.Frame(self.frame_soluzione)
        frame_bottoni.pack(fill=tk.X, padx=5, pady=2)
        self.bottone_risolvi = ttk.Button(frame_bottoni, text="RISOLVI", command=self.risolvi)
        self.bottone_risolvi.pack(side=tk.LEFT, padx=5)
        self.bottone_esporta = ttk.Button(frame_bottoni, text="ESPORTA IN PDF", command=self.esporta_pdf)
        self.bottone_esporta.pack(side=tk.LEFT, padx=5)
        self.bottone_copia = ttk.Button(frame_bottoni, text="COPIA", command=self.copia_soluzione)
        self.bottone_copia.pack(side=tk.LEFT, padx=5)
        self.var_mostra_quadranti = tk.BooleanVar(value=False)
        self.check_quadranti = ttk.Checkbutton(
            frame_bottoni, text="MOSTRA QUADRANTI", variable=self.var_mostra_quadranti
        )
        self.check_quadranti.pack(side=tk.LEFT, padx=5)
        self.bottone_grafico = ttk.Button(frame_bottoni, text="GRAFICO", command=self.mostra_grafico)
        self.bottone_grafico.pack(side=tk.LEFT, padx=5)
        self.last_soluzione = None
        self.last_A = None
        self.last_b = None
        self.last_c = None
        self.last_tabella = None
        self.last_nomi_var = None
        self.last_passi = None
        self.last_pivot = None
        self.last_img_grafico = None  # PIL Image
        self.aggiorna_variabili()

    def aggiorna_variabili(self):
        num_var = self.num_var.get()
        # Coefficienti FO
        for widget in self.frame_fo_coeff.winfo_children():
            widget.destroy()
        self.fo_coeff = []
        self.fo_segn = []
        for i in range(num_var):
            segno = tk.StringVar(value="+")
            combo_segno = ttk.Combobox(
                self.frame_fo_coeff, textvariable=segno, values=["+", "-"], width=3, state="readonly"
            )
            combo_segno.grid(row=0, column=i * 3, padx=2, pady=2)
            self.fo_segn.append(segno)
            segno.trace_add('write', lambda *args: self.mostra_fo())
            coeff_entry = ttk.Entry(self.frame_fo_coeff, width=5)
            coeff_entry.insert(0, "1.0")
            coeff_entry.grid(row=0, column=i * 3 + 1, padx=2, pady=2)
            coeff_entry.bind('<KeyRelease>', lambda e: self.mostra_fo())
            self.fo_coeff.append(coeff_entry)
            label_var = ttk.Label(self.frame_fo_coeff, text="x" + pedice_numero(i + 1))
            label_var.grid(row=0, column=i * 3 + 2, padx=2, pady=2)
        self.mostra_fo()

        # Vincoli di non negatività
        for widget in self.frame_nn.winfo_children():
            widget.destroy()
        self.nn_segn = []
        self.nn_val = []
        for i in range(num_var):
            label_var = ttk.Label(self.frame_nn, text="x" + pedice_numero(i + 1) + " ≥")
            label_var.grid(row=0, column=i * 2, padx=2, pady=2)
            self.nn_segn.append(">=")
            entry_val = ttk.Entry(self.frame_nn, width=5)
            entry_val.insert(0, "0")
            entry_val.grid(row=0, column=i * 2 + 1, padx=2, pady=2)
            self.nn_val.append(entry_val)

        # Vincoli
        for frame, segni, entrate, ineq, rhs_entry in self.vincoli:
            for widget in frame.winfo_children():
                widget.destroy()
            segni.clear()
            entrate.clear()
            for j in range(num_var):
                segno = tk.StringVar(value="+")
                combo_segno = ttk.Combobox(
                    frame, textvariable=segno, values=["+", "-"], width=3, state="readonly"
                )
                combo_segno.grid(row=0, column=j * 3, padx=2, pady=2)
                segni.append(segno)
                entrata = ttk.Entry(frame, width=5)
                entrata.insert(0, "1.0")
                entrata.grid(row=0, column=j * 3 + 1, padx=2, pady=2)
                entrate.append(entrata)
                label_var = ttk.Label(frame, text="x" + pedice_numero(j + 1))
                label_var.grid(row=0, column=j * 3 + 2, padx=2, pady=2)
            combo_ineq = ttk.Combobox(
                frame, textvariable=ineq, values=["≤", "≥", "="], width=3, state="readonly"
            )
            combo_ineq.grid(row=0, column=num_var * 3, padx=2, pady=2)
            rhs_entry.grid(row=0, column=num_var * 3 + 1, padx=2, pady=2)

    def mostra_fo(self):
        num_var = self.num_var.get()
        termini = []
        for i in range(num_var):
            segno = self.fo_segn[i].get() if i < len(self.fo_segn) else "+"
            coeff = self.fo_coeff[i].get() if i < len(self.fo_coeff) else "1"
            var = "x" + pedice_numero(i + 1)
            try:
                coeff_val = float(coeff)
            except Exception:
                coeff_val = 0
            if coeff_val == 0:
                continue
            if i == 0:
                termine = f"{'-' if segno == '-' else ''}{abs(coeff_val):g}{var}" if abs(coeff_val) != 1 else f"{'-' if segno == '-' else ''}{var}"
            else:
                termine = f" {'-' if segno == '-' else '+'} {abs(coeff_val):g}{var}" if abs(coeff_val) != 1 else f" {'-' if segno == '-' else ''}{var}"
            termini.append(termine)
        if not termini:
            termini = ["0"]
        tipo = self.tipo_fo.get()
        fo_str = ("max" if tipo == "max" else "min") + " z = " + "".join(termini)
        self.label_fo_valore.config(text=fo_str)

    def aggiungi_vincolo(self):
        num_var = self.num_var.get()
        frame = ttk.Frame(self.frame_vincoli)
        frame.pack(pady=2, padx=2)
        segni = []
        entrate = []
        for i in range(num_var):
            segno = tk.StringVar(value="+")
            combo_segno = ttk.Combobox(
                frame, textvariable=segno, values=["+", "-"], width=3, state="readonly"
            )
            combo_segno.grid(row=0, column=i * 3, padx=2, pady=2)
            segni.append(segno)
            entrata = ttk.Entry(frame, width=5)
            entrata.insert(0, "1.0")
            entrata.grid(row=0, column=i * 3 + 1, padx=2, pady=2)
            entrate.append(entrata)
            label_var = ttk.Label(frame, text="x" + pedice_numero(i + 1))
            label_var.grid(row=0, column=i * 3 + 2, padx=2, pady=2)
        ineq = tk.StringVar(value="≤")
        combo_ineq = ttk.Combobox(
            frame, textvariable=ineq, values=["≤", "≥", "="], width=3, state="readonly"
        )
        combo_ineq.grid(row=0, column=num_var * 3, padx=2, pady=2)
        rhs_entry = ttk.Entry(frame, width=5)
        rhs_entry.insert(0, "0.0")
        rhs_entry.grid(row=0, column=num_var * 3 + 1, padx=2, pady=2)
        self.vincoli.append((frame, segni, entrate, ineq, rhs_entry))

    def rimuovi_vincolo(self):
        if self.vincoli:
            frame, *_ = self.vincoli[-1]
            frame.destroy()
            self.vincoli.pop()

    def resetta(self):
        self.num_var.set(2)
        self.tipo_fo.set("max")
        for widget in self.frame_vincoli.winfo_children():
            widget.destroy()
        self.vincoli = []
        self.aggiorna_variabili()
        self.text_soluzione.delete(1.0, tk.END)
        self.last_soluzione = None
        self.last_A = None
        self.last_b = None
        self.last_c = None

    def salva_problema(self):
        dati = self._ottieni_dati_problema()
        filename = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("File JSON", "*.json")],
            title="Salva problema"
        )
        if filename:
            with open(filename, "w", encoding="utf-8") as f:
                json.dump(dati, f, indent=2)
            messagebox.showinfo("Salva problema", "Problema salvato!")

    def carica_problema(self):
        filename = filedialog.askopenfilename(
            filetypes=[("File JSON", "*.json")],
            title="Carica problema"
        )
        if filename:
            with open(filename, "r", encoding="utf-8") as f:
                dati = json.load(f)
            self._imposta_dati_problema(dati)
            messagebox.showinfo("Carica problema", "Problema caricato!")

    def _ottieni_dati_problema(self):
        num_var = self.num_var.get()
        tipo = self.tipo_fo.get()
        fo_coeffs = [(segn.get(), entry.get()) for segn, entry in zip(self.fo_segn, self.fo_coeff)]
        nn_segn = self.nn_segn
        nn_val = [entry.get() for entry in self.nn_val]
        vincoli = []
        for _, segni, entrate, ineq, rhs_entry in self.vincoli:
            vincolo = {
                "segni": [s.get() for s in segni],
                "valori": [e.get() for e in entrate],
                "ineq": ineq.get(),
                "rhs": rhs_entry.get()
            }
            vincoli.append(vincolo)
        return {
            "num_var": num_var,
            "tipo": tipo,
            "fo_coeffs": fo_coeffs,
            "nn_segn": nn_segn,
            "nn_val": nn_val,
            "vincoli": vincoli
        }

    def _imposta_dati_problema(self, dati):
        self.num_var.set(dati.get("num_var", 2))
        self.tipo_fo.set(dati.get("tipo", "max"))
        self.aggiorna_variabili()
        for i, coeff_pair in enumerate(dati.get("fo_coeffs", [])):
            if isinstance(coeff_pair, (list, tuple)) and len(coeff_pair) == 2:
                s, c = coeff_pair
            else:
                s, c = "+", coeff_pair
            if i < len(self.fo_segn):
                self.fo_segn[i].set(s)
            if i < len(self.fo_coeff):
                self.fo_coeff[i].delete(0, tk.END)
                self.fo_coeff[i].insert(0, c)
        for i, v in enumerate(dati.get("nn_val", [])):
            if i < len(self.nn_val):
                self.nn_val[i].delete(0, tk.END)
                self.nn_val[i].insert(0, v)
        for widget in self.frame_vincoli.winfo_children():
            widget.destroy()
        self.vincoli = []
        for vincolo in dati.get("vincoli", []):
            self.aggiungi_vincolo()
            frame, segni, entrate, ineq, rhs_entry = self.vincoli[-1]
            for j, s in enumerate(vincolo.get("segni", [])):
                if j < len(segni):
                    segni[j].set(s)
            for j, val in enumerate(vincolo.get("valori", [])):
                if j < len(entrate):
                    entrate[j].delete(0, tk.END)
                    entrate[j].insert(0, val)
            ineq.set(vincolo.get("ineq", "≤"))
            rhs_entry.delete(0, tk.END)
            rhs_entry.insert(0, vincolo.get("rhs", "0.0"))
        self.mostra_fo()

    def risolvi(self):
        try:
            self.text_soluzione.delete("1.0", tk.END)
            num_var = self.num_var.get()
            tipo = self.tipo_fo.get()
            fo_coeffs = [(self.fo_segn[i].get(), self.fo_coeff[i].get()) for i in range(num_var)]
            nn_segn = self.nn_segn
            nn_val = [entry.get() for entry in self.nn_val]
            vincoli = []
            for _, segni, entrate, ineq, rhs_entry in self.vincoli:
                vincolo = {
                    "segni": [s.get() for s in segni],
                    "valori": [e.get() for e in entrate],
                    "ineq": ineq.get(),
                    "rhs": rhs_entry.get()
                }
                vincoli.append(vincolo)
            # FIX: rimuovi 'return_pivots', chiama simplex_full col numero corretto di parametri
            passi, sol, ottimo, tabella, nomi_var, pivot = lp_solver.simplex_full(tipo, fo_coeffs, vincoli, nn_segn, nn_val)
            for passo in passi:
                self.text_soluzione.insert(tk.END, passo + "\n\n")
            self.last_soluzione = sol
            self.last_tabella = tabella
            self.last_nomi_var = nomi_var
            self.last_passi = passi
            self.last_pivot = pivot if len(pivot) > 0 else None
            self.last_c = [float(c[1]) if c[0] == "+" else -float(c[1]) for c in fo_coeffs]
            self.last_A = []
            self.last_b = []
            for v in vincoli:
                riga = []
                for s, val in zip(v["segni"], v["valori"]):
                    riga.append(float(val) if s == "+" else -float(val))
                self.last_A.append(riga)
                self.last_b.append(float(v["rhs"]))
        except Exception as e:
            messagebox.showerror("Errore", f"Errore durante la risoluzione: {e}")

    def mostra_grafico(self):
        if self.last_A is None or self.last_b is None or self.last_c is None:
            messagebox.showwarning("Grafico", "Risolvere il problema prima di generare il grafico.")
            return
        popup = tk.Toplevel(self.master)
        popup.title("Grafico - Regione Ammissibile")
        popup.geometry("820x820")
        fig, ax = plt.subplots(figsize=(8, 8))
        variabili = [f"x{pedice_numero(i + 1)}" for i in range(self.num_var.get())]
        if len(self.last_c) == 2:
            graph.plot_solution_to_axes(
                ax, np.array(self.last_c), np.array(self.last_A), np.array(self.last_b), variabili,
                show_all_quadrants=self.var_mostra_quadranti.get(),
                solution=self.last_soluzione, show_solution=True
            )
            canvas = FigureCanvasTkAgg(fig, master=popup)
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, popup)
            toolbar.update()
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=150)
            buf.seek(0)
            from PIL import Image as PILImage
            self.last_img_grafico = PILImage.open(buf).copy()
            buf.close()
            plt.close(fig)
        else:
            ax.text(0.5, 0.5, "Grafico disponibile solo per 2 variabili.", ha="center", va="center", fontsize=16)
            fig.tight_layout()
            canvas = FigureCanvasTkAgg(fig, master=popup)
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            canvas.draw()

    def copia_soluzione(self):
        try:
            self.master.clipboard_clear()
            self.master.clipboard_append(self.text_soluzione.get("1.0", tk.END))
            messagebox.showinfo("Copiato", "Soluzione copiata negli appunti!")
        except Exception:
            messagebox.showerror("Errore", "Impossibile copiare la soluzione negli appunti.")

    def esporta_pdf(self):
        try:
            filename = filedialog.asksaveasfilename(
                defaultextension=".pdf",
                filetypes=[("File PDF", "*.pdf")],
                title="Salva Soluzione in PDF"
            )
            if not filename:
                return
            vincoli_export = []
            for _, segni, entrate, ineq, rhs_entry in self.vincoli:
                vincoli_export.append({
                    "segni": [s.get() for s in segni],
                    "valori": [e.get() for e in entrate],
                    "ineq": ineq.get(),
                    "rhs": rhs_entry.get()
                })
            variabili = [f"x{pedice_numero(i + 1)}" for i in range(self.num_var.get())]
            if self.last_img_grafico is None and self.last_A is not None and self.last_b is not None and self.last_c is not None:
                fig, ax = plt.subplots(figsize=(8, 8))
                graph.plot_solution_to_axes(
                    ax, np.array(self.last_c), np.array(self.last_A), np.array(self.last_b), variabili,
                    show_all_quadrants=self.var_mostra_quadranti.get(),
                    solution=self.last_soluzione, show_solution=True
                )
                buf = io.BytesIO()
                fig.savefig(buf, format='png', dpi=150)
                buf.seek(0)
                from PIL import Image as PILImage
                self.last_img_grafico = PILImage.open(buf).copy()
                buf.close()
                plt.close(fig)
            pdf_exporter.esporta_in_pdf(
                filename=filename,
                variabili=variabili,
                tipo_fo=self.tipo_fo.get(),
                fo_coeff=[entry.get() for entry in self.fo_coeff],
                fo_segn=[sign.get() for sign in self.fo_segn],
                vincoli=vincoli_export,
                nn_segn=self.nn_segn,
                nn_val=[entry.get() for entry in self.nn_val],
                soluzione=self.last_soluzione,
                tabella=self.last_tabella,
                nomi_var=self.last_nomi_var,
                passi=self.last_passi,
                pivot=self.last_pivot,
                mostra_quadranti=self.var_mostra_quadranti.get(),
                c=self.last_c,
                A=self.last_A,
                b=self.last_b,
                img_grafico=self.last_img_grafico
            )
            messagebox.showinfo("Successo", "PDF salvato con successo.")
        except Exception as e:
            messagebox.showerror("Errore PDF", f"Errore durante l'esportazione del PDF: {e}")

if __name__ == "__main__":
    root = tk.Tk()
    gui = RisolutorePLGUI(root)
    root.mainloop()