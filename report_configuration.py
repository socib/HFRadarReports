import matplotlib.dates as md


class _Settings(object):
    """Global Variables used in the report generation.

    The __init__ method uses the year (YYYY) and month (MM) as input. It is specified at the end of the file.

    Note:
        To use these settings, just import reportConfiguration as cfg and then just access it via cfg.settings.year etc.

    Examples:
        import reportConfiguration as cfg
        print cfg.settings.month_str

    """
    @property
    def xfmt(self):
        return self._XFMT

    @property
    def variables_of_interest(self):
        return self._VARIABLES_OF_INTEREST

    @property
    def buoy_variables_of_interest(self):
        return self._BUOY_VARIABLES_OF_INTEREST

    @property
    def form_path(self):
        return self._form_path

    @property
    def galf_path(self):
        return self._galf_path

    @property
    def totals_path(self):
        return self._totals_path

    @property
    def filesize_threshold(self):
        return [self._galf_filesize_threshold, self._form_filesize_threshold, self._total_filesize_threshold]

    @property
    def statistics_variable_names(self):
        return self._statistics_variable_names

    @property
    def threshold_parameters(self):
        return self._threshold_parameters

    @property
    def close_to_buoy_statistics_variable_names(self):
        return self._close_to_buoy_statistics_variable_names

    @property
    def buoy_paths(self):
        return self._buoy_paths

    @property
    def document_output_directory(self):
        return self._document_output_directory

    @property
    def socib_logo_path(self):
        return self._socib_logo_path

    @property
    def logging_path(self):
        return self._logging_path

    @property
    def octave_modified_t_tide_path(self):
        return self._octave_modified_t_tide_path

    def __init__(self):
        self._XFMT = md.DateFormatter('%Y-%m-%d')
        self._VARIABLES_OF_INTEREST = [
          "time",
          "LAT",
          "LON",
          "U",
          "V",
          "WSPE",
          "WSPE_DIR",
          "SSN",
          "RADV",
          "RABA_GALF",
          "RABA_FORM",
          "RABA_DIFF_GALF",
          "RABA_DIFF_FORM"
        ]
        self._BUOY_VARIABLES_OF_INTEREST = [
            "time",
            "LAT",
            "LON",
            "CUR_DIR",
            "CUR_SPE"
        ]
        self._form_path = '/home/radar/data_rt/radar_system_eivissa/SCB-CODARSSSITE002/raw_archive/'
        self._galf_path = '/home/radar/data_rt/radar_system_eivissa/SCB-CODARSSSITE001/raw_archive/'
        self._totals_path = '/home/radar/data_rt/radar_system_eivissa/SCB-CODARSSPROC001/raw_archive/'
        self._galf_filesize_threshold = 165.0
        self._form_filesize_threshold = 110.0
        self._total_filesize_threshold = 60
        self._statistics_variable_names = ['RADV', 'SSN', 'RABA_GALF', 'RABA_FORM', 'RABA_DIFF_GALF', 'RABA_DIFF_FORM']
        self._threshold_parameters = {'RADV': [500, 3000], 'SSN': [20, None], 'RABA_GALF': [-130, -90],
                                      'RABA_DIFF_GALF': [0, 20], 'RABA_FORM': [-110, -70], 'RABA_DIFF_FORM': [0, 20]}
        self._close_to_buoy_statistics_variable_names = ['U', 'V', 'WSPE', 'WSPE_DIR', 'U_QAL', 'V_QAL', 'COVARIANCE_QAL']
        self._buoy_paths = {
            'folder': ['mooring/currentmeter', 'mooring/currentmeter'],
            'sub_folder': ['buoy_canaldeibiza-scb_dcs002', 'buoy_canaldeibiza-scb_fsi002'],
            'station_name': ['buoy-canaldeibiza_scb-dcs002', 'buoy-canaldeibiza_scb-fsi002']
        }
        self._document_output_directory = 'Reports/'
        self._socib_logo_path = 'logo_socib.eps'
        self._logging_path = 'HFRadar.log'
        self._octave_modified_t_tide_path = 't_tide_octave'

settings = _Settings()
