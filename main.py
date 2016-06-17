from HFRadar import HFRadar
from pylatex import Document, Package, NoEscape, Section, Itemize
from pylatex.base_classes import Options
import sys


class HFReportGenerator:
    def __init__(self, year, month):
        self.year = int(year)
        self.month = int(month)
        self.doc = None
        self.set_up_document()
        HFRadar(self.year, self.month, self.doc)
        self.doc.generate_pdf(extra_compiler_args='--xelatex')

    def write_introduction(self):
        with self.doc.create(Section('Introduction', numbering=False)):
            self.doc.append('The Fixed Stations monthly reports are generated automatically and are based on the'
                            ' information available in the thredds server for the fixed stations managed by SOCIB.')
            self.doc.append('The selection of variables is based upon the list of important variables from the SOCIB'
                            ' DataDiscovery Service.')
            self.doc.append('For each station, the report presents:')
            with self.doc.create(Itemize()) as itemize:
                itemize.add_item('A data summary describing the data availability, quality flags, mean and extremal'
                                 ' values over the period of interest')
                itemize.add_item('Plots of time series for the selected variables: time series with all the data points'
                                 ' with their associated quality flag and time series displaying only the good data')

    def write_title_page(self):
        self.doc.append(NoEscape(r'\thispagestyle{empty}\begin{titlepage}\centering'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\vspace*{3.5cm}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'{\Huge \thetitle}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\vspace{1cm}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\includegraphics[height=1.5cm]{logo_socib.eps}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\vspace{1cm}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'{\large\theauthor}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\vspace{1cm}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'{\large \thedate}'))
        self.doc.append(NoEscape(r'\end{titlepage}'))
        self.doc.append(NoEscape(r'\pagestyle{fancy}'))

    def write_header_preamble(self):
        self.doc.documentclass.options = Options('a4paper', 'svgnames', 'final')
        self.doc.packages.append(Package('hyperref'))
        self.doc.packages.append(Package('caption', options=['justification=centering']))
        self.doc.packages.append(Package('titling'))
        self.doc.packages.append(Package('sectsty'))
        self.doc.packages.append(Package('xcolor'))
        self.doc.packages.append(Package('booktabs'))
        self.doc.packages.append(Package('url,doi'))
        self.doc.packages.append(Package('lastpage'))
        self.doc.packages.append(Package('helvet'))
        self.doc.packages.append(Package('fancyhdr'))
        self.doc.packages.append(Package('datetime'))
        self.doc.packages.append(Package('graphicx'))

        self.doc.preamble.append(NoEscape(r'\definecolor{bluesocib}{rgb}{0.3411,0.7804,0.804}'
                                          r'\sectionfont{\color{bluesocib}}'))
        self.doc.preamble.append(NoEscape(r'\title{SOCIB Fixed Station Data Report}'))
        self.doc.preamble.append(NoEscape(r'\author{SOCIB Data Center}'))
        self.doc.preamble.append(NoEscape(r'\date{\today}'))
        self.doc.preamble.append(NoEscape(r'\ddmmyyyydate'))
        self.doc.preamble.append(NoEscape(r'\hypersetup{bookmarksopen=true,bookmarksnumbered=true,  pdffitwindow=false,'
                                          r' pdfstartview=FitH,pdftoolbar=true,pdfmenubar=true,pdfwindowui=true,'
                                          r'pdfauthor=\theauthor,pdftitle=\thetitle,pdfsubject=SOCIB Documentation,'
                                          r'bookmarksopenlevel=2,colorlinks=true,breaklinks=true,'
                                          r'linkcolor=black,anchorcolor=black,citecolor=black,filecolor=black,'
                                          r'menucolor=black,urlcolor=black}'))

        self.doc.preamble.append(NoEscape(r'\DeclareGraphicsExtensions{.jpg,.pdf,.png,.eps}'
                                          r'\renewcommand{\captionfont}{\it \small}'))
        self.doc.preamble.append(NoEscape(r'\setlength{\textwidth}{16cm}\setlength{\textheight}{23.5cm}'
                                          r'\setlength{\headheight}{50.3pt}\setlength{\footskip}{45pt}'
                                          r'\setlength{\hoffset}{-1.5cm}\setlength{\voffset}{-2.5cm}'
                                          r'\setlength{\unitlength}{1cm}\setlength{\parindent}{0pt}'
                                          r'\parskip 0.25cm\setcounter{tocdepth}{2}'
                                          r'\fancyhf{}\lhead{ \fancyplain{}{\hspace{-1.cm}'
                                          r'\includegraphics[height=1.25cm]{logo_socib.eps}} }'
                                          r'\cfoot{\hspace*{.9\textwidth}\sffamily \bfseries '
                                          r'\thepage\ / \pageref{LastPage}}\renewcommand{\headrulewidth}{0pt}  '
                                          r'\fancyfoot[L]{    '
                                          r'\parbox[b]{\dimexpr\linewidth\relax}'
                                          r'{\sffamily \bfseries \hspace*{.05\textwidth}'
                                          r'\textcolor{bluesocib}{www.socib.es} \hfill \\    '
                                          r'{\color{bluesocib}\rule{\dimexpr\linewidth\relax}{2pt}}\\}}'
                                          r'\fancyfootoffset{0.2\textwidth}'))

    def set_up_document(self):
        self.doc = Document('SOCIB_HFRadar_Report' + str(self.year) + str(self.month).zfill(2))
        self.write_header_preamble()
        self.write_title_page()
        self.doc.append(NoEscape(r'\pagebreak'))
        self.write_introduction()
        self.doc.append(NoEscape(r'\pagebreak'))
        self.doc.append(NoEscape(r'\tableofcontents'))
        self.doc.append(NoEscape(r'\pagebreak'))
        with self.doc.create(Section('test')):
            self.doc.append('Test')


def main():
    """
    Requires year and month as input parameters.
    year: YYYY
    month: MM

    Example: python main 2016 3
    """
    HFReportGenerator(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
