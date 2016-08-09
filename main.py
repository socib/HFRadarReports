from report_utils import get_years_and_months_ranges
import report_configuration as c
from HFRadar import HFRadar
from pylatex import Document, Package, NoEscape, Section, Itemize
from pylatex.base_classes import Options
import calendar
import sys


class HFReportGenerator:
    def __init__(self, year, month):
        self.year = int(year)
        self.month = int(month)
        self.month_str = calendar.month_name[self.month]
        self.doc = None
        self.set_up_document()
        HFRadar(self.year, self.month, self.doc)
        # self.doc.generate_tex()
        self.doc.generate_pdf(extra_compiler_args='--xelatex')

    def write_introduction(self):
        with self.doc.create(Section('Introduction', numbering=False)):
            self.doc.append('This monthly report aims to extract useful and meaningful information from the HF radar'
                            ' data using qualitative and quantitative data analysis methods. The document is'
                            ' automatically generated based on the information available in the THREDDS Data Server'
                            ' (TDS) Catalog for the HF Radar managed by SOCIB.')
            self.doc.append('Automatic data processing includes:')
            with self.doc.create(Itemize()) as itemize:
                itemize.add_item('Monthly means of the direction vectors, statistics (time series and data tables)')
                itemize.add_item(NoEscape(r'Comparisons of the horizontal current components derived from HF radar and'
                                          r' the pointwise subsurface currents from the current-meter (1.5 m) deployed'
                                          r' in the Ibiza Channel (at location 38$^{\circ}$49.46$^\prime$N and'
                                          r' 0$^{\circ}$47.02$^\prime$W), which allow us to evaluate the radar'
                                          r' performance and identify temporal periods or malfunctioning of the radar'
                                          r' (or the current-meter).'))
                itemize.add_item('Note that figures are using the oceanographic convention (currents pointing in the'
                                 ' direction the flow is toward). ')

    def write_title_page(self):
        self.doc.append(NoEscape(r'\thispagestyle{empty}\begin{titlepage}\centering'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\vspace*{3.5cm}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'{\Huge \thetitle}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\vspace{1cm}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\includegraphics[height=1.5cm]{' + c.settings.socib_logo_path + r'}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\vspace{1cm}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'{\large\theauthor}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'\vspace{1cm}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'{\LARGE ' + self.month_str + ' ' + str(self.year) + r' Monthly Report}'))
        self.doc.append(NoEscape(r''))
        self.doc.append(NoEscape(r'{\normalsize Document generated on \thedate}'))
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
        self.doc.preamble.append(NoEscape(r'\title{SOCIB HF Radar Data Report}'))
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
                                          r'\includegraphics[height=1.25cm]{' + c.settings.socib_logo_path + r'}} }'
                                          r'\cfoot{\hspace*{.9\textwidth}\sffamily \bfseries '
                                          r'\thepage\ / \pageref{LastPage}}\renewcommand{\headrulewidth}{0pt}  '
                                          r'\fancyfoot[L]{    '
                                          r'\parbox[b]{\dimexpr\linewidth\relax}'
                                          r'{\sffamily \bfseries \hspace*{.05\textwidth}'
                                          r'\textcolor{bluesocib}{www.socib.es} \hfill \\    '
                                          r'{\color{bluesocib}\rule{\dimexpr\linewidth\relax}{2pt}}\\}}'
                                          r'\fancyfootoffset{0.2\textwidth}'))

    def set_up_document(self):
        self.doc = Document(c.settings.document_output_directory + 'SOCIB_HFRadar_Report' + str(self.year) +
                            str(self.month).zfill(2))
        self.write_header_preamble()
        self.write_title_page()
        self.doc.append(NoEscape(r'\pagebreak'))
        self.write_introduction()
        self.doc.append(NoEscape(r'\pagebreak'))
        self.doc.append(NoEscape(r'\tableofcontents'))
        self.doc.append(NoEscape(r'\pagebreak'))


def main():
    """
    Requires year and month as input parameters.
    year: YYYY
    month: MM
    Example:
        python main 2016 3
    OR
        python main 2013 2016 8 7
        will process all monts in the timespan
    """
    if len(sys.argv) == 3:
        # single station
        HFReportGenerator(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 5:
        # process time span
        years, months = get_years_and_months_ranges(*map(int, sys.argv[1:5]))
        for i in range(0, len(years)):
            try:
                HFReportGenerator(years[i], months[i])
            except RuntimeError:
                continue


if __name__ == "__main__":
    main()
