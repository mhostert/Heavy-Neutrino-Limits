import pylatex
from pylatex.utils import NoEscape
import os
import requests

def default_noentry(key):
    ''' Default entry for bibentry not found in Inspire '''
    key = key.replace(' ', '').lower()
    return '@misc{'+key+',\n\tnote\t= "{No entry found in Inspire for '+key+'}"\n}\n'

def add_plot(doc, case):
    ''' add a Figure environment with the caption as the citations '''
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


def create_latex_doc(cases, TEX_PATH = 'tex_files/'):
    ''' Create a latex document with the plots and the references '''
    
    if not os.path.exists(TEX_PATH):
        os.makedirs(TEX_PATH)

    doc = pylatex.Document(f'{TEX_PATH}/mixing_plots', documentclass=NoEscape(r'revtex4-1'))

    for case in cases:
        add_plot(doc, case)

    with open(f'{TEX_PATH}/mixing_plots.bib', 'w', encoding='utf-8') as f:
        added_ids = []
        for case in cases:
            for ref in case.limits.reference:
                if not (ref in added_ids):
                    # get the bib entry from inspire
                    response = requests.get(f"https://inspirehep.net/api/literature?q=texkeys:{ref}&format=bibtex")
                    
                    # check api found an entry and that it's not empty
                    if response.status_code == 200 and response.content:
                        f.write((response.content).decode("utf-8") )
                    else:
                        # write a note as a bib entry
                        ne = default_noentry(ref)
                        f.write(ne)
                        print(f"Could not find Inspire entry for texkey={ref}. Adding @misc entry instead.")
                    added_ids.append(ref)


    doc.append(pylatex.Command('bibliographystyle', arguments=NoEscape(r'apsrev4-1')))
    doc.append(pylatex.Command('bibliography', arguments=NoEscape(r'mixing_plots')))

    doc.generate_pdf(clean_tex=True)
    doc.generate_tex()