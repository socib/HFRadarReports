# SOCIB HF Radar Report Generator

The SOCIB HF Radar Report Generator is a python command line tool to create automatically the monthly reports of the HF Radar managed by SOCIB. The monthly reports aim to extract useful and meaningful information from the HF radar data using qualitative and quantitative data analysis methods. The documents are automatically generated based on the information available in the THREDDS Data Server (TDS) Catalog for the HF Radar managed by SOCIB.

The automatic data processing includes:
- Monthly means of the direction vectors, statistics (time series and data tables) 
- Comparisons of the horizontal current components derived from HF radar and the pointwise subsurface currents from the currentmeter (1.5 m) deployed in the Ibiza Channel (at location 38 ◦ 49.46 N and 0 ◦ 47.02 W), which allow us to evaluate the radar performance and identify temporal periods or malfunctioning of the radar (or the current-meter).
- Note that figures are using the oceanographic convention (currents pointing in the direction the flow is toward).

## Dependencies
- numpy
- netCDF4
- matplotlib
- pylatex
- requests
- pandas
- pytz
- scipy
- oct2py
- windrose
- adjustText
- okean

Furthermore, mpl_toolkits.basemap (download here: https://sourceforge.net/projects/matplotlib/files/matplotlib-toolkits/), octave and latexmk is required.

The package has been developed and tested using Ubuntu 14 and is currently only supported by python 2.7.

## Install
Clone this repository using git. To install the dependencies, run pip install -r requirement.txt or install them manually.

Inside your octave session (access by typing octave into your terminal), the following packages also need to be installed with the following commands:

>pkg install -forge specfun

>pkg install -forge control

>pkg install -forge general

>pkg install -forge signal

Last but not least, to start processing, run the reportGeneratorHFRadar.py (see next section).

## Usage
There are two possibilities to use the tool:
- a) Specify the year and month of the reportGeneratorHFRadar.py to process the month (i.e. python reportGeneratorHFRadar.py 2016 7)
- b)Specify the timespan you want to process of the reportGeneratorHFRadar.py to process the timespan (i.e. "python reportGeneratorHFRadar.py 2015 2016 2 7" will process all monthly datasets between 02.2015 and 07.2016)
