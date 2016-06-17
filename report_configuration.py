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


settings = _Settings()
