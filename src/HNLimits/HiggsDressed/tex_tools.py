import pylatex
from pylatex.utils import NoEscape
import os
import requests

def add_plot(doc, case):
    # caption with the citations
    citations = fr"Constraints on {case.latexflavor} as a function of the HNL mass $m_N$. Limits shown: "
    names = sorted(case.limits.plot_label)
    for name in names:
        inspires_key = list(case.limits.reference[case.limits.plot_label == name])[0].replace(' ', '')
        citations += fr'{name}~\cite{{{inspires_key}}}, '.replace('\\\\','')
    citations = citations[:-2]+'.'

    # figure
    with doc.create(pylatex.Figure(position='h!')) as latexfig:
        latexfig.add_image(case.path_to_plot, width=NoEscape(r'1\textwidth'))
        latexfig.add_caption(NoEscape(citations))

# Basic document
def create_latex_doc(cases, TEX_PATH = 'tex_files/'):

    if not os.path.exists(TEX_PATH):
        os.makedirs(TEX_PATH)

    doc = pylatex.Document(f'{TEX_PATH}/hd_plots', documentclass=NoEscape(r'revtex4-1'))

    for case in cases:
        add_plot(doc, case)

    with open(f'{TEX_PATH}/hd_plots.bib','w', encoding='utf-8') as f:
        added_ids = []
        for case in cases:
            for ref in case.limits.reference:
                if not ref in added_ids:
                    response = requests.get(f"https://inspirehep.net/api/literature?q=texkeys:{ref}&format=bibtex")
                    if response.status_code == 200:
                        f.write((response.content).decode("utf-8") )
                        added_ids.append(ref)
                    else:
                        print(f"Could not find Inspire entry for texkey={ref}. Skipping it.")


    doc.append(pylatex.Command('bibliographystyle', arguments=NoEscape(r'apsrev4-1')))
    doc.append(pylatex.Command('bibliography', arguments=NoEscape(r'hd_plots')))

    doc.generate_pdf(clean_tex=True)
    doc.generate_tex()