from report_utils import *
from netCDF4 import Dataset
import report_configuration as c
import numpy as np
import logging
import pytz
import calendar
from datetime import timedelta
from pylatex import NoEscape, Section

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
        self.month_str = calendar.month_name[self.month]
        self.link = []
        self.root = None
        self.buoy_root = None
        self.converted_time = []
        self.x_limits = []
        self.closest_lat_idx = []
        self.closest_lon_idx = []
        self.lat = []
        self.lon = []
        self.time = []
        self.ibiz_avail = []
        self.lat_lon_percent = []
        self.high_res_basemap = get_basemap('h')
        self.variables = {}
        self.buoy_variables = {}
        self.get_root()
        self.get_variables()
        self.convert_time()
        self.do_processing()

    def get_root(self):
        self.link = get_thredds_opendap_link('hf_radar', 'hf_radar_ibiza-scb_codarssproc001', 1, 'dep0001',
                                             'hf-radar-ibiza_scb-codarssproc001', self.year, self.month)
        try:
            self.root = Dataset(self.link)
        except RuntimeError:
            logger.error('File does not exist. ' + self.link, exc_info=True)
        buoy_link = get_thredds_opendap_link('mooring/currentmeter', 'buoy_canaldeibiza-scb_dcs002', 1, 'dep0001',
                                             'buoy-canaldeibiza_scb-dcs002', self.year, self.month)
        try:
            self.buoy_root = Dataset(buoy_link)
        except RuntimeError:
            logger.error('File does not exist. ' + buoy_link, exc_info=True)

    def read_variable(self, var_name):
        try:
            self.variables[var_name] = self.root.variables[var_name]
            if var_name == 'time':
                self.time = get_data_array(self.variables[var_name])
            elif var_name == 'LAT':
                self.lat = get_data_array(self.variables[var_name])
            elif var_name == 'LON':
                self.lon = get_data_array(self.variables[var_name])
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

    def convert_time(self):
        dates = [datetime.fromtimestamp(ts, tz=pytz.utc) for ts in self.time]
        self.converted_time = md.date2num(dates)
        self.x_limits = get_x_limits(self.time[0])

    def write_section(self, title, func):
        with self.doc.create(Section(title)):
            func()
            self.doc.append(NoEscape(r'\clearpage'))

    def do_processing(self):
        sections = [['Monthly Mean Direction Plot', self.monthly_mean],
                    ['Time Series at Closest Buoy Ibiza Grid Point', self.timeseries_at_buoy_ibiza],
                    ['Data Tables at Closest Buoy Ibiza Grid Point', self.tables_at_buoy_ibiza],
                    ['Comparison Graphs', self.comparison_radar_buoy],
                    ['Temporal Availability', self.temporal_availability],
                    ['Spatial Availability', self.spatial_availability],
                    ['Temporal and Spatial Availability', self.spatial_and_temporal_availability],
                    ['Percent Filesize above Threshold', self.filesize_threshold],
                    ['Statistics from QC Variables', self.compute_statistics],
                    ['Threshold Graphs', self.create_threshold_graphs],
                    ['Histogram Radial Files per 10 Days.', self.create_histogram]]
        for section in sections:
            self.write_section(section[0], section[1])

    def monthly_mean(self):
        logger.info('Generating Monthly Mean Direction Plot')
        self.doc.append('This plot represents the averaged monthly mean directions of the U and V components.')
        u_mean = get_temporal_mean_from_grid(get_data_array(self.variables["U"]))
        v_mean = get_temporal_mean_from_grid(get_data_array(self.variables["V"]))
        np_longrid, np_langrid = np.meshgrid(self.lon, self.lat)
        hf_monthly_mean_direction_plot(self.doc, u_mean, v_mean, np_longrid, np_langrid, self.high_res_basemap)

    def tables_at_buoy_ibiza(self):
        variables_of_interest = c.settings.close_to_buoy_statistics_variable_names
        self.doc.append('Summarising data tables for the variables at the closest grid point respective to the'
                        ' Ibiza Buoy.')
        for v in variables_of_interest:
            variable = self.root.variables[v]
            variable_data, variable_dim = self.get_closest_grid_data(variable)
            qc_v_name = get_qc_variable_name(variable)
            if qc_v_name is not None:
                qc_variable = self.root.variables[qc_v_name]
                qc_variable_data, qc_variable_dim = self.get_closest_grid_data(qc_variable)
                bad_idx = qc_variable_data != 1
                good_data = variable_data
                good_data[bad_idx] = np.nan
                if v != 'WSPE_DIR':
                    cur_mean, cur_mean_time, cur_min, cur_min_time, cur_max, cur_max_time = \
                        get_mean_min_max(good_data, self.converted_time)
                    sum_good, sum_prob_good, sum_prob_bad, sum_bad, sum_spike, sum_nan, percent_rest, \
                        percent_good, percent_nan = get_sums_and_percents(qc_variable_data)
                    write_table(self.doc, '|cc|c|cccc|c|c|ccc|', variable, self.month_str, self.year, cur_mean=cur_mean,
                                cur_max=cur_max, cur_min=cur_min, cur_max_time=cur_max_time, cur_min_time=cur_min_time,
                                sum_good=sum_good, sum_prob_good=sum_prob_good, sum_prob_bad=sum_prob_bad,
                                sum_bad=sum_bad, sum_spike=sum_spike, sum_nan=sum_nan, percent_rest=percent_rest,
                                percent_good=percent_good, percent_nan=percent_nan)
                else:
                    sum_good, sum_prob_good, sum_prob_bad, sum_bad, sum_spike, sum_nan, percent_rest, \
                        percent_good, percent_nan = get_sums_and_percents(qc_variable_data)
                    write_table(self.doc, '|cc|c|cccc|c|', variable, self.month_str, self.year, sum_good=sum_good,
                                sum_prob_good=sum_prob_good, sum_prob_bad=sum_prob_bad, sum_bad=sum_bad,
                                sum_spike=sum_spike, sum_nan=sum_nan, percent_rest=percent_rest,
                                percent_good=percent_good, percent_nan=percent_nan)
            else:
                if v != 'WSPE_DIR':
                    cur_mean, cur_mean_time, cur_min, cur_min_time, cur_max, cur_max_time = \
                        get_mean_min_max(variable_data, self.converted_time)
                    write_table(self.doc, '|cc|ccc|', variable, self.month_str, self.year, cur_mean=cur_mean,
                                cur_min=cur_min, cur_min_time=cur_min_time, cur_max=cur_max, cur_max_time=cur_max_time)

    def get_closest_grid_data(self, variable):
        if get_dimension(variable) == 1:
            return get_data_array(variable), 1
        elif get_dimension(variable) == 3:
            return get_data_array(variable)[:, self.closest_lat_idx, self.closest_lon_idx], 3

    def timeseries_at_buoy_ibiza(self):
        # Note: they are sharing the amplifiers here currently
        buoy_lat, buoy_lon = get_data_array(self.buoy_variables["LAT"]), get_data_array(self.buoy_variables["LON"])
        self.closest_lat_idx, self.closest_lon_idx = get_idx_closest_grid_point(self.lat, self.lon, buoy_lat, buoy_lon)
        self.doc.append('Closest grid point LAT: %.6f' % self.lat[self.closest_lat_idx] + ' degree\n')
        self.doc.append('Closest grid point LON: %.6f' % self.lon[self.closest_lon_idx] + ' degree\n\n')
        variables_of_interest = ['U', 'V', 'WSPE', 'WSPE_DIR', 'U_QAL', 'V_QAL', 'COVARIANCE_QAL']
        self.doc.append('Represents the variables ' + str(variables_of_interest) + ' at the closest grid point'
                                                                                   ' respective to the Ibiza Buoy and,'
                                                                                   ' if available, the corresponding'
                                                                                   ' buoy data.')
        for variable_name in variables_of_interest:
            cur_variable = self.root.variables[variable_name]
            cur_qc_variable_name = get_qc_variable_name(cur_variable)
            if cur_qc_variable_name is not None:
                cur_qc_variable = self.root.variables[cur_qc_variable_name]
                cur_qc_data = get_data_array(cur_qc_variable)[:, self.closest_lat_idx, self.closest_lon_idx]
            else:
                cur_qc_data = None
            cur_data = get_data_array(cur_variable)[:, self.closest_lat_idx, self.closest_lon_idx]
            if variable_name == 'WSPE_DIR':
                amplifier_variable = get_data_array(
                    self.variables['WSPE'])[:, self.closest_lat_idx, self.closest_lon_idx]
                #plot_quiver_direction(self.doc, self.converted_time, cur_data, 'quiver wspe dir', amplifier_variable,
                #                      input_month_title=self.month_str + ' ' + str(self.year))
                # Let us bypass the buoy graph here for visual comparison)
                buoy_dir = self.buoy_root.variables["CUR_DIR"]
                buoy_dir_data = get_data_array(buoy_dir)
                buoy_spe = self.buoy_root.variables["CUR_SPE"]
                buoy_spe_data = get_data_array(buoy_spe)
                plot_quiver_direction_overlapping(self.doc, self.converted_time, cur_data,
                                                  'Comparison of directions from HF Radar closest grid point (bottom)'
                                                  ' and Ibiza Buoy (top)', buoy_dir_data,
                                                  lower_amplifier=amplifier_variable,
                                                  input_month_title=self.month_str + ' ' + str(self.year),
                                                  upper_amplifier=amplifier_variable)
                #plot_quiver_direction(self.doc, self.converted_time, buoy_dir_data, 'quiver cur dir buoy',
                #                      buoy_spe_data, input_month_title=self.month_str + ' ' + str(self.year))
                if cur_qc_data is not None:
                    # laziness: no checks performed that variable exists
                    cur_data_good_idx = cur_qc_data == 1
                    amplifier_qc_variable_name = get_qc_variable_name(self.variables['WSPE'])
                    amplifier_qc_data = get_data_array(self.root.variables[amplifier_qc_variable_name])[:, self.closest_lat_idx, self.closest_lon_idx]
                    amplifier_qc_good_idx = amplifier_qc_data == 1
                    combined_data_good_idx = np.logical_or(cur_data_good_idx, amplifier_qc_good_idx)
                    combined_data = cur_data[combined_data_good_idx]
                    #plot_quiver_direction(self.doc, self.converted_time[combined_data_good_idx], combined_data,
                    #                      'quiver wspe dir good data only', amplifier_variable[combined_data_good_idx],
                    #                      input_month_title=self.month_str + ' ' + str(self.year))
                buoy_qc_variable_name = get_qc_variable_name(buoy_dir)
                buoy_amplifier_qc_variable_name = get_qc_variable_name(buoy_spe)
                buoy_qc_variable_data = get_data_array(self.buoy_root.variables[buoy_qc_variable_name])
                buoy_qc_amplifier_data = get_data_array(self.buoy_root.variables[buoy_amplifier_qc_variable_name])
                buoy_dir_good_idx = buoy_qc_variable_data == 1
                buoy_amplifier_good_idx = buoy_qc_amplifier_data == 1
                buoy_combined_good_idx = np.logical_or(buoy_dir_good_idx, buoy_amplifier_good_idx)
                buoy_combined_data = buoy_dir_data[buoy_combined_good_idx]
                buoy_combined_amplifier_data = buoy_spe_data[buoy_combined_good_idx]
                hf_buoy_combined_idx = np.logical_and(combined_data_good_idx, buoy_combined_good_idx)
                buoy_combined_data = buoy_dir_data[hf_buoy_combined_idx]
                hf_combined_data = cur_data[hf_buoy_combined_idx]
                # plot_quiver_direction(self.doc, self.converted_time[buoy_combined_good_idx], buoy_combined_data,
                #                       'quiver cur dir buoy good data only', buoy_combined_amplifier_data,
                #                       input_month_title=self.month_str + ' ' + str(self.year))
                plot_quiver_direction_overlapping(self.doc, self.converted_time[hf_buoy_combined_idx], hf_combined_data,
                                                  'Comparison of directions from HF Radar closest grid point (bottom)'
                                                  ' and Ibiza Buoy (top) GOOD DATA ONLY', buoy_combined_data,
                                                  lower_amplifier=amplifier_variable[hf_buoy_combined_idx],
                                                  input_month_title=self.month_str + ' ' + str(self.year),
                                                  upper_amplifier=amplifier_variable[hf_buoy_combined_idx])
            else:
                plot_1d(self.doc, self.converted_time, cur_data, cur_variable.units, variable_name, cur_qc_data,
                        input_month_title=self.month_str + ' ' + str(self.year))
                if cur_qc_data is not None:
                    good_idx = cur_qc_data == 1
                    good_data = cur_data
                    good_data[np.logical_not(good_idx)] = np.nan
                    plot_1d(self.doc, self.converted_time, good_data, cur_variable.units,
                            variable_name + ' good data only', input_month_title=self.month_str + ' ' + str(self.year))

    def comparison_radar_buoy(self):
        self.doc.append('The following figures are showing the speed and direction observed by the SOCIB HF Radar of'
                        ' the closest grid point with respect to the position of the SOCIB Ibiza Buoy.\n')
        self.doc.append('On a rotating basis, an overlapping graph and the differences between these two datasets are'
                        ' shown.\n')
        self.doc.append('The compared datasets are WSPE --> CUR_SPE and WSPE_DIR --> CUR_DIR where the WSPE variables'
                        ' are measured by the HF radar and CUR variables the respective buoy variables.')
        temp_dates = [datetime.fromtimestamp(ts, tz=pytz.utc) for ts in get_data_array(self.buoy_variables["time"])]
        buoy_converted_time = md.date2num(temp_dates)
        compare_variables_names = [["WSPE", "CUR_SPE"],
                                   ["WSPE_DIR", "CUR_DIR"]]
        for comparison_var_names in compare_variables_names:
            hf_variable = self.variables[comparison_var_names[0]]
            hf_units = hf_variable.units
            logger.debug(hf_units)
            hf_data = get_data_array(hf_variable)[:, self.closest_lat_idx, self.closest_lon_idx]
            buoy_variable = self.buoy_variables[comparison_var_names[1]]
            buoy_units = buoy_variable.units
            logger.debug(buoy_units)
            inverse_conversion_factor = 1
            if hf_units == 'm s-1' and buoy_units == 'cm s-1':
                inverse_conversion_factor = 100
            elif hf_units == 'cm s-1' and buoy_units == 'm s-1':
                inverse_conversion_factor = 1.0/100
            buoy_data = get_data_array(buoy_variable)/inverse_conversion_factor
            plot_overlapping_1d_graphs(self.doc, self.converted_time, hf_data, buoy_converted_time, buoy_data,
                                       title_str='Buoy {0} minus HF {1}'.format(comparison_var_names[0],
                                                                                comparison_var_names[1]),
                                       y_label=hf_units)
            diff = get_differences(hf_data, buoy_data)
            plot_1d(self.doc, self.converted_time, diff, hf_units,
                    'Buoy {0} minus HF {1}'.format(comparison_var_names[0], comparison_var_names[1]),
                    input_month_title=str(self.year) + ' ' + str(self.month))
            if hf_units == 'degree' or buoy_units == 'degree':
                amplifier = get_data_array(self.variables["WSPE"][:, self.closest_lat_idx, self.closest_lon_idx])
                amplifier[np.isnan(amplifier)] = 0
                plot_quiver_direction(self.doc, self.converted_time, diff,
                                      'Buoy {0} minus HF {1}'.format(comparison_var_names[0], comparison_var_names[1]),
                                      amplifier, input_month_title=str(self.year) + ' ' + str(self.month))

    def temporal_availability(self):
        self.doc.append('This graph shows the temporal availability of both radial sites managed by SOCIB. A continues'
                        ' line indicates full data availability whilst an interrupted line indicates missing files for'
                        ' these dates.')
        form_path = c.settings.form_path + str(self.year) + '/' + str(self.month).zfill(2) + '/'
        galf_path = c.settings.galf_path + str(self.year) + '/' + str(self.month).zfill(2) + '/'
        self.ibiz_avail = plot_hf_temporal_availability(self.doc, form_path, galf_path, self.time, self.year,
                                                        self.month)

    def spatial_availability(self):
        self.doc.append('These maps show the total coverage of available data at each gridpoint.')
        time_ibiz_avail = self.time[np.logical_not(self.ibiz_avail)]
        wspe_ibiz_avail = get_data_array(self.variables["WSPE"])[np.logical_not(self.ibiz_avail), :, :]
        lat_lon_counts = np.meshgrid(np.zeros((1, len(self.lon))), np.zeros((1, len(self.lat))))[0]
        for i in range(len(time_ibiz_avail)):
            cur_wspe = wspe_ibiz_avail[i, :, :]
            cur_data_avail_idx = ~np.isnan(cur_wspe)
            lat_lon_counts[cur_data_avail_idx] += 1
        self.lat_lon_percent = lat_lon_counts/float(len(time_ibiz_avail))*100

        lon2, lat2 = np.meshgrid(self.lon, self.lat)
        masked_array = np.ma.array(self.lat_lon_percent, mask=self.lat_lon_percent == 0.0)

        temp_basemap = get_basemap('h')
        temp_basemap1 = get_basemap('h')

        plot_spatial_availability(self.doc, temp_basemap, lon2, lat2, masked_array, temp_basemap1)

    def spatial_and_temporal_availability(self):
        self.doc.append('This figure shows the temporal and spatial availability of all gridpoints that contain at once'
                        ' during the regarded month data.')
        spatial_avail = np.reshape(self.lat_lon_percent, (1, len(self.lon)*len(self.lat)))[0]
        sorted_spatial_avail = np.sort(spatial_avail)
        non_zero_sorted_spatial_avail = sorted_spatial_avail[sorted_spatial_avail != 0.]
        non_zero_sorted_spatial_avail = non_zero_sorted_spatial_avail[::-1]
        sum_is_nan_points = sum(self.ibiz_avail)
        non_zero_sorted_spatial_avail = np.append(non_zero_sorted_spatial_avail, np.zeros((1, sum_is_nan_points)))
        temp_avail = np.array(map(float, np.arange(0, len(non_zero_sorted_spatial_avail))))
        temporal_avail_percent = temp_avail/len(temp_avail)*100
        plot_temporal_and_spatial_availability(self.doc, temporal_avail_percent, non_zero_sorted_spatial_avail,
                                               temp_avail)

    def filesize_threshold(self):
        self.doc.append('Represents the percent availability of files within the regarded month that exceed the defined'
                        ' thresholds for file sizes. Missing files are not regarded.')
        form_path = c.settings.form_path + str(self.year) + '/' + str(self.month).zfill(2) + '/'
        galf_path = c.settings.galf_path + str(self.year) + '/' + str(self.month).zfill(2) + '/'
        form_elements = sorted(os.listdir(form_path))
        galf_elements = sorted(os.listdir(galf_path))
        thresholds = c.settings.filesize_threshold
        galf_good, form_good, total_good = get_hf_radial_sites_file_sizes(galf_path, galf_elements, form_path,
                                                                          form_elements, thresholds[0], thresholds[1],
                                                                          thresholds[2])
        combined_results = [[thresholds[0], galf_good, 'Galfi'],
                            [thresholds[1], form_good, 'Formentera'],
                            [thresholds[2], total_good, 'Totals']]
        for result in combined_results:
            self.doc.append('{0} percent files above threshold ('.format(result[2]) + str(
                int(result[0])) + ')Kb: %.2f%%' % result[1])
            self.doc.append(NoEscape(r'\\'))

    def compute_statistics(self):
        self.doc.append('The following tables show the mean, standard deviation, minimum, maximum and percentage of'
                        ' good data with respect to their associated QC variable.')
        for v in c.settings.statistics_variable_names:
            v = self.root.variables[v]
            v_data = get_data_array(v)
            v_qc_name = get_qc_variable_name(v)
            if v_qc_name is not None:
                v_qc_data = get_data_array(self.root.variables[v_qc_name])
            else:
                v_qc_data = None
            v_mean, v_std, v_min, v_max, v_good_percent = get_hf_simple_statistics(v_data, v_qc_data)
            write_table(self.doc, '|c|c|c|c|c|', v, self.month_str, self.year, v_mean=v_mean, v_std=v_std, v_min=v_min,
                        v_max=v_max, v_good_percent=v_good_percent)

    def create_threshold_graphs(self):
        self.doc.append('The following figures show the time series and their acceptable values for the defined quality'
                        ' control parameters.')
        for name, thresholds in c.settings.threshold_parameters.iteritems():
            data = get_data_array(self.root.variables[name])
            plot_threshold_graph(self.doc, self.converted_time, data, thresholds[0], thresholds[1],
                                 self.root.variables[name])

    def create_histogram(self):
        self.doc.append('This bar chart shows the number of available radial files per 10 days.')
        cur_end_of_month = calendar.monthrange(self.year, self.month)[1]
        file_prefix = 'RDLi_'
        file_suffix = '.ruv'
        stations_subfolders = {'GALF': c.settings.galf_path, 'FORM': c.settings.form_path}
        start_date = datetime(self.year, self.month, cur_end_of_month)
        stations_bins = {'GALF': [], 'FORM': []}
        stations_end_time = []
        for station in stations_subfolders:
            temp_station_base = stations_subfolders[station]
            for i in range(10, 360, 10):
            #for i in range(10, 40, 10):
                cur_date = start_date - timedelta(days=i)
                stations_end_time.append(cur_date + timedelta(days=10))
                current_bin_file_counter = 0
                for k in range(1, 11, 1):
                    cur_temp_date = cur_date + timedelta(days=k)
                    cur_temp_year = cur_temp_date.year
                    cur_temp_month = cur_temp_date.month
                    cur_temp_day = cur_temp_date.day
                    temp_station_folder = temp_station_base + '/' + str(cur_temp_year) + '/' + str(
                        cur_temp_month).zfill(2) + '/' + file_prefix + station + '_'
                    for j in range(0, 24, 1):
                        cur_hour = j
                        file_date_identifier = str(cur_temp_year) + '_' + str(cur_temp_month).zfill(2) + '_' + str(
                            cur_temp_day).zfill(2) + '_' + str(cur_hour).zfill(2) + '00'
                        temp_file_path = temp_station_folder + file_date_identifier + file_suffix
                        if not os.path.isfile(temp_file_path):
                            pass
                        else:
                            current_bin_file_counter += 1
                stations_bins[station].append(current_bin_file_counter)
        plot_histogram_radial_files(self.doc, stations_end_time, stations_bins)

    def create_power_spectrum(self):
        pass