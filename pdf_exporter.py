from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
import numpy as np
import utils  # Absolute import
import graph  # Absolute import

def export_to_pdf(
    filename, variables, obj_type, obj_coeffs, obj_signs, constraints, nn_signs, nn_values,
    integer_vars, solution, slack_values, dual_solution, alternative_solutions,
    original_constraint_count, problem_history, show_all_quadrants, c, A, steps_text
):
    """Export the solution and graph to PDF using ReportLab."""
    try:
        doc = SimpleDocTemplate(filename)
        styles = getSampleStyleSheet()
        story = []

        # Title
        title_style = ParagraphStyle(
            name='TitleStyle',
            parent=styles['Title'],
            fontSize=18,
            spaceAfter=12
        )
        story.append(Paragraph("<b>Soluzione del Problema di Programmazione Lineare</b>", title_style))
        story.append(Spacer(1, 0.2 * inch))

        # Problem Definition
        story.append(Paragraph("<b>Definizione del Problema:</b>", styles['Heading2']))
        story.append(Spacer(1, 0.1 * inch))

        # Objective Function
        obj_str = f"{obj_type} z = "
        for i, (sign, coeff, var) in enumerate(zip(obj_signs, obj_coeffs, variables)):
            if i > 0:
                obj_str += f" {sign} "
            obj_str += f"{coeff} {var}"
        story.append(Paragraph(f"Funzione Obiettivo: {obj_str}", styles['Normal']))
        story.append(Spacer(1, 0.1 * inch))

        # Constraints
        story.append(Paragraph("Vincoli:", styles['Normal']))
        for frame, signs, entries, ineq, rhs_entry in constraints:
            constraint_str = ""
            for i, (sign, entry, var) in enumerate(zip(signs, entries, variables)):
                if i > 0:
                    constraint_str += f" {sign} "
                constraint_str += f"{entry.get()} {var}"
            constraint_str += f" {ineq.get()} {rhs_entry.get()}"
            story.append(Paragraph(constraint_str, styles['Normal']))
        story.append(Spacer(1, 0.1 * inch))

        # Non-Negativity Constraints
        story.append(Paragraph("Vincoli di Non-Negatività:", styles['Normal']))
        for sign, val, var in zip(nn_signs, nn_values, variables):
            story.append(Paragraph(f"{var} {sign} {val}", styles['Normal']))
        story.append(Spacer(1, 0.2 * inch))

        # Solution
        story.append(Paragraph("<b>Soluzione Ottima:</b>", styles['Heading2']))
        story.append(Spacer(1, 0.1 * inch))

        if solution is not None:
            sol_str = ", ".join([f"{var} = {val:.4f}" for var, val in zip(variables, solution)])
            story.append(Paragraph(sol_str, styles['Normal']))
            story.append(Paragraph(f"Valore Ottimo: {np.dot(c, solution):.4f}", styles['Normal']))

            # Slack Variables
            if slack_values is not None:
                story.append(Paragraph("Variabili di Slack:", styles['Normal']))
                slack_str = ", ".join([f"s{i+1} = {val:.4f}" for i, val in enumerate(slack_values)])
                story.append(Paragraph(slack_str, styles['Normal']))

            # Dual Solution
            if dual_solution is not None:
                story.append(Paragraph("Soluzione Duale:", styles['Normal']))
                dual_str = ", ".join([f"y{i+1} = {val:.4f}" for i, val in enumerate(dual_solution)])
                story.append(Paragraph(dual_str, styles['Normal']))

            # Alternative Solutions
            if alternative_solutions:
                story.append(Paragraph("Soluzioni Ottime Alternative:", styles['Normal']))
                for alt_sol in alternative_solutions:
                    alt_sol_str = ", ".join([f"{var} = {val:.4f}" for var, val in zip(variables, alt_sol)])
                    story.append(Paragraph(alt_sol_str, styles['Normal']))

            # Simplex Tableau History
            if problem_history:
                story.append(Paragraph("<b>Passaggi dell'Algoritmo del Simplesso:</b>", styles['Heading2']))
                story.append(Spacer(1, 0.1 * inch))

                for title, tableau in problem_history:
                    story.append(Paragraph(f"<b>{title}</b>", styles['Heading3']))
                    story.append(Spacer(1, 0.05 * inch))
                    n_vars = len(variables)
                    n_slacks = tableau.shape[1] - n_vars - 1 - (tableau.shape[0] - 1)
                    data = (
                        ["Base"]
                        + [f"x<sub>{i+1}</sub>" for i in range(n_vars)]
                        + [f"s<sub>{i+1}</sub>" for i in range(n_slacks)]
                        + ["RHS"]
                    )
                    table_data = [data]
                    for i in range(tableau.shape[0]):
                        if i < tableau.shape[0] - 1:
                            # Find the basic variable for this row
                            basic_var_index = None
                            for j in range(n_vars + n_slacks):
                                col = tableau[:tableau.shape[0] - 1, j]
                                if np.allclose(col, np.eye(tableau.shape[0] - 1)[:, i]):
                                    basic_var_index = j
                                    break
                            if basic_var_index is not None:
                                if basic_var_index < n_vars:
                                    basic_var = f"x<sub>{basic_var_index + 1}</sub>"
                                else:
                                    basic_var = f"s<sub>{basic_var_index - n_vars + 1}</sub>"
                            else:
                                basic_var = ""
                            row = [basic_var] + [f"{x:.4f}" for x in tableau[i].tolist()]
                            table_data.append(row)
                        else:
                            row = ["z"] + [f"{x:.4f}" for x in tableau[i].tolist()]
                            table_data.append(row)

                    table = Table(table_data)
                    table.setStyle(TableStyle([
                        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                        ('GRID', (0, 0), (-1, -1), 1, colors.black),
                        ('LEFTPADDING', (0, 0), (-1, -1), 5),
                        ('RIGHTPADDING', (0, 0), (-1, -1), 5),
                    ]))
                    story.append(table)
                    story.append(Spacer(1, 0.2 * inch))
        else:
            story.append(Paragraph("Nessuna soluzione trovata.", styles['Normal']))

        doc.build(story)
        return True
    except Exception as e:
        print(f"Error exporting to PDF: {e}")
        return False