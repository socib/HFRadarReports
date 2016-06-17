from report_utils import *
from netCDF4 import Dataset
import report_configuration as c
import numpy as np
import logging
from pylatex import Document, Package, NoEscape, Section, Itemize

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
handler = logging.FileHandler('HF_radar.log')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class HFRadar:
    def __init__(self, year, month, doc):
        self.doc = doc
        self.year = year
        self.month = month
        self.link = []
        self.root = []
        self.buoy_root = []
        logger.info('Generating High Resolution Basemap')
        self.high_res_basemap = get_basemap('h')
        self.variables = {}
        self.buoy_variables = {}
        self.get_root()
        self.get_buoy_root()
        self.get_variables()
        self.do_processing()

    def get_root(self):
        self.link = get_thredds_opendap_link('hf_radar', 'hf_radar_ibiza-scb_codarssproc001', 1, 'dep0001',
                                             'hf-radar-ibiza_scb-codarssproc001', self.year, self.month)
        try:
            self.root = Dataset(self.link)
        except RuntimeError:
            logger.error('File does not exist. ' + self.link, exc_info=True)

    def get_buoy_root(self):
        buoy_link = get_thredds_opendap_link('mooring/currentmeter', 'buoy_canaldeibiza-scb_dcs002', 1, 'dep0001',
                                             'buoy-canaldeibiza_scb-dcs002', self.year, self.month)
        try:
            self.buoy_root = Dataset(buoy_link)
        except RuntimeError:
            logger.error('File does not exist. ' + buoy_link, exc_info=True)

    def read_variable(self, var_name):
        try:
            self.variables[var_name] = self.root.variables[var_name]
        except RuntimeError:
            logger.error('No variable exists with that name: ' + var_name, exc_info=True)

    def read_buoy_variable(self, var_name):
        try:
            self.buoy_variables[var_name] = self.buoy_root.variables[var_name]
        except RuntimeError:
            logger.error('No variable exists with that name: ' + var_name, exc_info=True)

    def get_variables(self):
        for var_name in c.settings.variables_of_interest:
            self.read_variable(var_name)
        for var_name in c.settings.buoy_variables_of_interest:
            self.read_buoy_variable(var_name)

    def do_processing(self):
        with self.doc.create(Section('Monthly Mean Direction Plot')):
            self.monthly_mean()

    def monthly_mean(self):
        logger.info('Generating Monthly Mean Direction Plot')
        u_mean = get_temporal_mean_from_grid(get_data_array(self.variables["U"]))
        v_mean = get_temporal_mean_from_grid(get_data_array(self.variables["V"]))
        np_longrid, np_langrid = np.meshgrid(get_data_array(self.variables["LON"]), get_data_array(self.variables["LAT"]))
        hf_monthly_mean_direction_plot(self.doc, u_mean, v_mean, np_longrid, np_langrid, self.high_res_basemap)

