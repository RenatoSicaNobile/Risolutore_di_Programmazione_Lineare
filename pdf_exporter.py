from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image as RLImage, Table, TableStyle, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch, cm
from reportlab.lib import colors
import tempfile

def esporta_in_pdf(
    filename, variabili, tipo_fo, fo_coeff, fo_segn, vincoli, nn_segn, nn_val,
    soluzione, tabella, nomi_var, passi, pivot, mostra_quadranti, c, A, b, img_grafico
):
    try:
        doc = SimpleDocTemplate(filename, pagesize=letter)
        styles = getSampleStyleSheet()
        story = []

        stile_titolo = ParagraphStyle(
            name='StileTitolo',
            parent=styles['Title'],
            fontSize=18,
            spaceAfter=12
        )
        story.append(Paragraph("<b>Soluzione del Problema di Programmazione Lineare</b>", stile_titolo))
        story.append(Spacer(1, 0.2 * inch))

        story.append(Paragraph("<b>Definizione del Problema:</b>", styles['Heading2']))
        story.append(Spacer(1, 0.1 * inch))

        fo_str = f"{tipo_fo} z = "
        for i, (sign, coeff, var) in enumerate(zip(fo_segn, fo_coeff, variabili)):
            if i > 0:
                fo_str += f" {sign} "
            fo_str += f"{coeff} {var}"
        story.append(Paragraph(f"Funzione Obiettivo: {fo_str}", styles['Normal']))
        story.append(Spacer(1, 0.1 * inch))

        story.append(Paragraph("<b>Vincoli:</b>", styles['Heading3']))
        for vin in vincoli:
            vincolo_str = ""
            for i, (sign, val, var) in enumerate(zip(vin["segni"], vin["valori"], variabili)):
                if i > 0:
                    vincolo_str += f" {sign} "
                vincolo_str += f"{val} {var}"
            vincolo_str += f" {vin['ineq']} {vin['rhs']}"
            story.append(Paragraph(vincolo_str, styles['Normal']))
        story.append(Spacer(1, 0.1 * inch))

        story.append(Paragraph("<b>Vincoli di Non-Negatività:</b>", styles['Normal']))
        for sign, val, var in zip(nn_segn, nn_val, variabili):
            story.append(Paragraph(f"{var} {sign} {val}", styles['Normal']))
        story.append(Spacer(1, 0.2 * inch))

        if img_grafico is not None:
            story.append(Paragraph("<b>Grafico della Regione Ammissibile e Soluzione:</b>", styles['Heading2']))
            with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmpfile:
                img_grafico.save(tmpfile.name)
                story.append(RLImage(tmpfile.name, width=13*cm, height=13*cm))
            story.append(Spacer(1, 0.2 * inch))

        story.append(Paragraph("<b>Soluzione Ottima:</b>", styles['Heading2']))
        story.append(Spacer(1, 0.1 * inch))

        if soluzione is not None:
            sol_str = ", ".join([f"{var} = {val:.4f}" for var, val in zip(variabili, soluzione)])
            story.append(Paragraph(sol_str, styles['Normal']))
            if c is not None and soluzione is not None:
                story.append(Paragraph(f"Valore Ottimo: {sum(c[i]*soluzione[i] for i in range(len(c))):.4f}", styles['Normal']))
            story.append(Spacer(1, 0.1 * inch))
        else:
            story.append(Paragraph("Nessuna soluzione ottima trovata.", styles['Normal']))

        story.append(PageBreak())
        story.append(Paragraph("<b>Passaggi Matematici e Motivazioni:</b>", styles['Heading2']))
        story.append(Spacer(1, 0.1 * inch))

        if passi is not None:
            for idx, passo in enumerate(passi):
                html_passo = "<br/>".join(passo.split("\n"))
                story.append(Paragraph(f"<b>Passo {idx+1}:</b><br/>{html_passo}", styles['Normal']))
                story.append(Spacer(1, 0.1 * inch))

        if tabella is not None and nomi_var is not None:
            story.append(Paragraph("<b>Tableau Finale:</b>", styles['Heading3']))
            dati_tabella = [nomi_var] + [[f"{x:.4f}" if isinstance(x, float) else str(x) for x in riga] for riga in tabella]
            t = Table(dati_tabella, hAlign='LEFT')
            t.setStyle(TableStyle([('BACKGROUND', (0,0), (-1,0), colors.grey),
                                   ('TEXTCOLOR',(0,0),(-1,0),colors.whitesmoke),
                                   ('ALIGN',(0,0),(-1,-1),'CENTER'),
                                   ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                                   ('BOTTOMPADDING', (0,0), (-1,0), 12),
                                   ('BACKGROUND',(0,1),(-1,-1),colors.beige)]))
            story.append(t)
            story.append(Spacer(1, 0.2 * inch))

        if pivot:
            story.append(Paragraph("<b>Pivot scelti:</b>", styles['Heading3']))
            piv_str = "; ".join([str(p) for p in pivot])
            story.append(Paragraph(piv_str, styles['Normal']))

        story.append(Paragraph("<b>Motivazione:</b> L'algoritmo Simplex viene utilizzato per risolvere problemi di programmazione lineare trovando la soluzione ottima, passando da una base all'altra fino a che la funzione obiettivo non può più essere migliorata. I passaggi sopra illustrano tutte le iterazioni con i dettagli dei pivot e dei tableaux.", styles['Normal']))
        story.append(Spacer(1, 0.2 * inch))

        doc.build(story)
    except Exception as e:
        import traceback
        print(traceback.format_exc())
        raise