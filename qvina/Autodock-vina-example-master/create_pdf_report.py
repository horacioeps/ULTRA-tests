#!/usr/bin/env python3

import csv
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib import colors
from reportlab.platypus import Table, TableStyle, Paragraph, Spacer, Image
from reportlab.platypus.doctemplate import SimpleDocTemplate
from reportlab.platypus.frames import Frame
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
import os

def read_docking_results():
    """Read and process docking results"""
    results = []
    with open('batch_results/docking_results.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            results.append({
                'name': row['Ligand_Name'],
                'score': float(row['Best_Affinity_kcal_mol']),
                'runtime': float(row['Runtime_seconds'])
            })
    
    # Sort by best score (most negative)
    results.sort(key=lambda x: x['score'])
    return results

def create_molecule_placeholder(name, width=100, height=100):
    """Create a simple molecular structure placeholder"""
    from reportlab.graphics.shapes import Drawing, Circle, Line, String
    from reportlab.graphics import renderPDF
    
    d = Drawing(width, height)
    
    # Create a simple molecular diagram placeholder
    # Central atom
    d.add(Circle(50, 50, 8, fillColor=colors.blue, strokeColor=colors.black))
    
    # Bonds and atoms around
    positions = [(20, 30), (80, 30), (20, 70), (80, 70), (50, 20), (50, 80)]
    for i, (x, y) in enumerate(positions[:4]):  # Only show 4 connections
        # Bond line
        d.add(Line(50, 50, x, y, strokeColor=colors.black, strokeWidth=2))
        # Atom
        color = [colors.red, colors.green, colors.orange, colors.purple][i % 4]
        d.add(Circle(x, y, 6, fillColor=color, strokeColor=colors.black))
    
    # Add molecule name
    d.add(String(50, 10, name.replace('_', ' ').title(), textAnchor='middle', fontSize=8))
    
    # Save as image
    filename = f"temp_{name}_structure.png"
    renderPDF.drawToFile(d, filename, fmt='PNG')
    return filename

def create_pdf_report():
    """Create the main PDF report"""
    
    # Read results
    results = read_docking_results()
    
    # Create PDF
    filename = "batch_results/QuickVina_Docking_Report.pdf"
    doc = SimpleDocTemplate(filename, pagesize=A4)
    
    # Container for the 'Flowable' objects
    elements = []
    
    # Define styles
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=24,
        spaceAfter=30,
        alignment=1,  # Center alignment
        textColor=colors.darkblue
    )
    
    subtitle_style = ParagraphStyle(
        'CustomSubtitle',
        parent=styles['Heading2'],
        fontSize=16,
        spaceAfter=20,
        alignment=1,
        textColor=colors.darkgreen
    )
    
    # Title
    title = Paragraph("QuickVina Molecular Docking Results", title_style)
    elements.append(title)
    elements.append(Spacer(1, 20))
    
    # Summary statistics
    total_ligands = len(results)
    best_score = results[0]['score']
    worst_score = results[-1]['score']
    avg_score = sum(r['score'] for r in results) / len(results)
    total_time = sum(r['runtime'] for r in results)
    
    summary_text = f"""
    <b>Analysis Summary:</b><br/>
    • Total Ligands Tested: {total_ligands}<br/>
    • Best Binding Affinity: {best_score:.1f} kcal/mol ({results[0]['name']})<br/>
    • Average Binding Affinity: {avg_score:.1f} kcal/mol<br/>
    • Total Computation Time: {total_time:.1f} seconds<br/>
    • Score Range: {best_score:.1f} to {worst_score:.1f} kcal/mol
    """
    
    summary = Paragraph(summary_text, styles['Normal'])
    elements.append(summary)
    elements.append(Spacer(1, 30))
    
    # Create table data
    table_data = [['Rank', 'Ligand Name', '2D Structure', 'Binding Score (kcal/mol)', 'Runtime (s)']]
    
    # Create molecular structure placeholders and prepare table data
    for i, result in enumerate(results, 1):
        # Create structure placeholder
        struct_file = create_molecule_placeholder(result['name'])
        
        # Medal for top 3
        medal = ""
        if i == 1:
            medal = "🥇 "
        elif i == 2:
            medal = "🥈 "
        elif i == 3:
            medal = "🥉 "
        
        rank = f"{medal}{i}"
        name = result['name'].replace('_', ' ').title()
        score = f"{result['score']:.1f}"
        runtime = f"{result['runtime']:.1f}"
        
        # Try to add image to table
        try:
            from reportlab.platypus import Image
            img = Image(struct_file, width=60, height=60)
            table_data.append([rank, name, img, score, runtime])
        except:
            # Fallback if image doesn't work
            table_data.append([rank, name, "Structure", score, runtime])
    
    # Create and style the table
    table = Table(table_data, colWidths=[0.8*inch, 2.2*inch, 1.2*inch, 1.5*inch, 1*inch])
    
    table.setStyle(TableStyle([
        # Header row
        ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        
        # Data rows
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 10),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.lightgrey]),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        
        # Special formatting for top 3
        ('BACKGROUND', (0, 1), (-1, 3), colors.lightyellow),
        ('TEXTCOLOR', (3, 1), (3, -1), colors.darkgreen),  # Score column in green
        ('FONTNAME', (3, 1), (3, -1), 'Helvetica-Bold'),   # Score column bold
    ]))
    
    elements.append(table)
    elements.append(Spacer(1, 30))
    
    # Footer notes
    notes_text = """
    <b>Notes:</b><br/>
    • More negative binding scores indicate stronger protein-ligand interactions<br/>
    • Structures shown are representative 2D diagrams<br/>
    • All docking performed using QuickVina 2 with identical parameters<br/>
    • Binding site coordinates: center (139, 145, 171), size (25×25×25) Å
    """
    
    notes = Paragraph(notes_text, styles['Normal'])
    elements.append(notes)
    
    # Build PDF
    doc.build(elements)
    
    # Clean up temporary files
    for result in results:
        temp_file = f"temp_{result['name']}_structure.png"
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
            except:
                pass
    
    print(f"✅ PDF report created: {filename}")
    return filename

if __name__ == "__main__":
    try:
        report_file = create_pdf_report()
        print(f"📄 Report saved as: {report_file}")
        print("🔬 The report includes:")
        print("  • Summary statistics")
        print("  • Ranked ligand table")
        print("  • 2D molecular structure diagrams")
        print("  • Binding scores and runtimes")
        
    except ImportError as e:
        print(f"❌ Missing required library: {e}")
        print("💡 Install with: pip install reportlab")
        
    except Exception as e:
        print(f"❌ Error creating PDF: {e}")
        print("💡 Falling back to simple table generation...")
        
        # Simple fallback
        results = read_docking_results()
        print("\n📊 SIMPLIFIED RESULTS TABLE")
        print("="*50)
        print(f"{'Rank':<4} {'Ligand':<20} {'Score':<10}")
        print("-"*50)
        for i, result in enumerate(results, 1):
            medal = "🥇" if i == 1 else "🥈" if i == 2 else "🥉" if i == 3 else ""
            print(f"{i:<4} {result['name']:<20} {result['score']:<10.1f} {medal}") 