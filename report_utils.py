import report_configuration as c
import numpy as np
import requests
import time
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib as mp
import logging
import os
import pytz
from textwrap import wrap
from datetime import datetime
from pylatex import Figure, NoEscape, Table, Tabular, MultiColumn, MultiRow
from mpl_toolkits.basemap import Basemap
from scipy import signal, fft, interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
# noinspection PyUnresolvedReferences
from ttide.t_tide import t_tide
from matplotlib.patches import Ellipse
import calendar
from oct2py import octave
from matplotlib import ticker
from windrose import WindroseAxes
from adjustText import adjust_text
from okean import gshhs

# Command does not really work currently, but let's keep it if we use a similar approach later
# mp.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
handler = logging.FileHandler('Utils.log')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s p%(process)s {%(pathname)s:%(lineno)d} - %(name)s -'
                              ' %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def get_years_and_months_ranges(start_year, end_year, start_month, end_month):
    year_difference = end_year - start_year
    month_difference = end_month - start_month
    counter_difference = year_difference * 12 + month_difference
    cur_year = start_year
    out_years = []
    out_months = []
    for x in range(0, counter_difference + 1):
        cur_month = (start_month + x - 1) % 12 + 1
        if ((start_month + x) % 12) == 1:
            cur_year += 1
        out_years.append(cur_year)
        out_months.append(cur_month)
    return out_years, out_months


def calc_distance(lat1, lon1, lat2, lon2):
    # Using haversine formula to calculate the spherical distance
    radius = 6370.
    dlat = np.deg2rad(lat2) - np.deg2rad(lat1)
    dlon = np.deg2rad(lon2) - np.deg2rad(lon1)
    a = np.power(np.sin(dlat/2.), 2) + np.cos(np.deg2rad(lat1)) * np.cos(np.deg2rad(lat2)) * np.power(np.sin(dlon/2.), 2)
    dist = 2. * radius * np.arcsin(np.sqrt(a))
    return dist


def t_tide_harmonic_analysis(doc, u_data, v_data, cur_time, year, month, latitude_array, longitude_array, qc_u, qc_v,
                             wdir, wspe, np_longrid, np_latgrid):
    general_shape = u_data.shape
    desired_constituents = ['K1  ', 'M2  ', 'S2  ']
    mp.rcParams['ytick.labelsize'] = 10
    mp.rcParams['xtick.labelsize'] = 10
    mp.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(3, 1, figsize=(7, 11), sharey=True)

    axes[0].set_title(desired_constituents[0])
    axes[1].set_title(desired_constituents[1])
    axes[2].set_title(desired_constituents[2])



    # basemaps = [Basemap(projection='cyl', llcrnrlat=38.35, urcrnrlat=39.05, llcrnrlon=0.55, urcrnrlon=1.45,
    #                     lat_ts=35., resolution='h', ax=axes[0])]
    x, y = gshhs.get_coastline(xlim=[0.55, 1.45], ylim=[38.35, 39.05], res='f')
    for a in axes:
        a.plot(x, y, 'k')
        a.set_xlim([0.55, 1.45])
        a.set_ylim([38.35, 39.05])
        a.set(adjustable='box-forced', aspect='equal')
        a.set_ylabel('$^{\circ}$N', rotation=0, horizontalalignment='right')
        a.set_xlabel('$^{\circ}$E')
    # longrid, latgrid = basemaps[0](np_longrid, np_latgrid)
    #
    # for cur_map in basemaps:
    #     cur_map.drawcoastlines(linewidth=.25, zorder=4)
    #     cur_map.drawmapboundary(fill_color='white')
    #     cur_map.drawparallels(np.arange(30., 40., 1.), labels=[True, False, True, False], zorder=2)
    #     cur_map.drawmeridians(np.arange(-10., 2., 1.), labels=[False, True, False, True], zorder=2)
    for m in range(0, general_shape[1]):
        cur_latitude = latitude_array[m]
        for n in range(0, general_shape[2]):
            cur_longitude = longitude_array[n]
            cur_u = u_data[:, m, n]
            cur_v = v_data[:, m, n]
            cur_qc_u = qc_u[:, m, n]
            cur_qc_u_idx = cur_qc_u == 1
            cur_u[~cur_qc_u_idx] = np.nan
            cur_qc_v = qc_v[:, m, n]
            cur_qc_v_idx = cur_qc_v == 1
            cur_v[~cur_qc_v_idx] = np.nan
            cur_u, cur_v = compute_u_v_components(wdir[:, m, n], wspe[:, m, n])
            u_percent_nan = len(np.where(np.isnan(cur_u))[0])/float(len(cur_u)) * 100
            v_percent_nan = len(np.where(np.isnan(cur_v))[0])/float(len(cur_v)) * 100
            if (u_percent_nan > 99) or (v_percent_nan > 99):
                logger.info('>60% NaNs found (index {0}, {1})'
                            ' - U NaNs: {2:.2g}%, V NaNs: {3:.2f}%'.format(m, n, u_percent_nan, v_percent_nan))
                continue
            # cur_u = cur_u - np.nanmean(cur_u)
            # cur_v = cur_v - np.nanmean(cur_v)
            cur_u_interpolated = linear_fill(cur_time, cur_u)
            cur_v_interpolated = linear_fill(cur_time, cur_v)
            complex_u_v = cur_u_interpolated + 1j * cur_v_interpolated
            # list of constituents used, frequency of tidal constituents (cycles/hr), tidal constituents with confidence
            #  intervals
            logger.info('Processing {0}, {1} ({2}, {3})'.format(cur_latitude, cur_longitude, m, n))
            try:
                octave.addpath('/home/akrietemeyer/Documents/MATLAB/t_tide_octave/t_tide_v1.3beta')
                # noinspection PyUnusedLocal
                [const_names, const_freqs,
                 tide_const, prediction] = octave.t_tide(complex_u_v, 'interval', 1,  'start_time',
                                                         get_md_datenum(cur_time)[0], 'latitude', cur_latitude,
                                                         'output', 'none')
                # In case the python t_tide implementation should be used:
                # const_names, const_freqs, tide_const, prediction = t_tide(complex_u_v, dt=1,
                #                                                           stime=get_md_datenum(cur_time)[0],
                #                                                           lat=cur_latitude, output=False)
            except TypeError:
                logger.warning('Type error at tidal analysis.')
                continue
            cur_counter = 0
            for constituent in desired_constituents:
                cur_axis = axes[cur_counter]
                # cur_axis.set_title(constituent)
                idx = np.where(const_names == constituent)[0]
                ell_params = tide_const[idx][0]
                cur_ellipse = Ellipse(xy=(np_longrid[m, n], np_latgrid[m, n]), width=ell_params[0], height=ell_params[2],
                                      angle=ell_params[4])
                cur_axis.add_artist(cur_ellipse)
                cur_ellipse.set_clip_box(cur_axis.bbox)
                cur_ellipse.set_facecolor('white')
                cur_counter += 1
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.25)
    with doc.create(Figure(position='!htbp')) as plot:
        plot.add_plot()
        plot.add_caption('Tidal Ellipses for the main tidal constituents (K1, M2 and S2) in {0} {1}.'.format(
            calendar.month_name[month], str(year)))
    plt.clf()
    plt.close('all')
    mp.rcParams.update(mp.rcParamsDefault)


def clean_direction_and_speed(direction, var):
    """
    Remove nan and var=0 values in the two arrays
    if a var (wind speed) is nan or equal to 0, this data is
    removed from var array but also from dir array
    :param var:
    :param direction:
    """
    dirmask = np.isfinite(direction)
    varmask = (var != 0 & np.isfinite(var))
    ind = dirmask*varmask
    return direction[ind], var[ind]


def wind_rose_hist(direction, var, nsector, normed=False):
    if len(var) != len(direction):
        raise(ValueError("var and direction must have same length"))

    angle = 360. / nsector

    dir_bins = np.arange(-angle / 2, 360. + angle, angle, dtype=np.float)
    dir_edges = dir_bins.tolist()
    dir_edges.pop(-1)
    dir_edges[0] = dir_edges.pop(-1)
    dir_bins[0] = 0.

    bins = np.arange(10, 70, 10)

    var_bins = np.asarray(bins).tolist()
    var_bins.append(np.inf)

    table = np.lib.twodim_base.histogram2d(x=var, y=direction, bins=[var_bins, dir_bins], normed=False)[0]
    # add the last value to the first to have the table of North winds
    table[:, 0] = table[:, 0] + table[:, -1]
    # and remove the last col
    table = table[:, :-1]
    if normed:
        table = table * 100 / table.sum()
    wd_freq = np.sum(table, axis=0)
    return wd_freq


def plot_wind_rose(doc, speed, direction, cur_title, month_str, year, cur_y_lim=None):
    clean_wdir, clean_wspe = clean_direction_and_speed(direction, speed)
    if np.all(clean_wspe <= 10):
        logger.warning('Wind rose creation skipped. All speed values are below 10 cm/s.')
        return None
    with doc.create(Figure(position='htbp')) as plot:
        ax = WindroseAxes.from_ax()
        ax.bar(clean_wdir, clean_wspe, nsector=32, normed=True, opening=0.8, edgecolor='white',
               bins=np.arange(10, 70, 10))
        ax.set_legend()
        ax.set_title("\n".join(wrap(cur_title, 50)), y=1.05)
        table = ax._info['table']
        wd_freq = np.sum(table, axis=0)
        # In case a defined limit for the distributions is wanted to be set
        if cur_y_lim is not None:
            ax.set_ylim(cur_y_lim)
            ax.yaxis.set_major_locator(ticker.AutoLocator())
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        plot.add_plot(width=NoEscape(r'0.6\textwidth'))
        plot.add_caption(cur_title + ' in ' + month_str + ' ' + str(year))
        plt.clf()
        plt.close('all')
    return wd_freq


def plot_wind_bars_distribution(doc, wd_freq, cur_title, month_str, year, cur_y_lim=None):
    # care, I assume 32 bins here!
    if wd_freq is None:
        logger.info('Wind bar distribution skipped.')
        return
    with doc.create(Figure(position='htbp')) as plot:
        plt.bar(np.arange(32), wd_freq, align='center', color='gray')
        xlabels = ('N', '', '', '', 'N-E', '', '', '', 'E', '', '', '', 'S-E', '', '', '',
                   'S', '', '', '', 'S-W', '', '', '', 'W', '', '', '', 'N-W', '', '', '')
        xticks = np.arange(32)
        plt.gca().set_xticks(xticks)
        plt.ylabel('%')
        plt.xlabel('Direction')
        plt.draw()
        plt.gca().set_xticklabels(xlabels)
        plt.title(cur_title)
        if cur_y_lim is not None:
            plt.ylim(cur_y_lim)
        plt.draw()
        plot.add_plot(width=NoEscape(r'0.6\textwidth'))
        plot.add_caption(cur_title + ' in ' + month_str + ' ' + str(year) + '.')
        plt.clf()
        plt.close('all')


def average_spatially(data):
    out_data = np.nanmean(data[0::, :, :], axis=(1, 2))
    return out_data


def filter_components(doc, data, cur_time, y_label):
    y_label = transform_y_label(y_label)
    conv_time = get_md_datenum(cur_time)
    x_limits = [conv_time[0], conv_time[-1]]
    desired_filters = np.asarray([33, 24, 12, 19])          # in hours
    desired_filters = 1.0/desired_filters              # in Hz
    variance_storage = [np.nanvar(data)]
    variance_strings = ['unfiltered Spatially Averaged', 'low-pass 33h', 'low-pass 24h', 'low-pass 12h', 'low-pass 19h']

    for band in desired_filters:
        filtered_data, _ = apply_low_pass_filter(band, data, cur_time)
        variance_storage.append(np.nanvar(filtered_data)*100.)
    for i in range(0, 4):
        doc.append(NoEscape('Variance of {0}: {1:.5f}{2}.\n'.format(variance_strings[i], variance_storage[i],
                                                                    '$(cm\, s^{-1})^2$')))

    for band in desired_filters:
        filtered_data, interpolated_data = apply_low_pass_filter(band, data, cur_time)
        nan_idx = np.isnan(data)
        # nan insertion
        filtered_data[nan_idx] = np.nan
        interpolated_data[nan_idx] = np.nan
        with doc.create(Figure(position='htbp')) as plot:
            plt.figure()
            plt.gcf().set_size_inches(11, 3)
            plt.plot(conv_time, filtered_data, '-k')
            data_axis = plt.gca()
            data_axis.xaxis.set_major_formatter(c.settings.xfmt)
            data_axis.set_ylabel(y_label, rotation=0, horizontalalignment='right')
            data_axis.grid(b=False, which='major', color='k', linestyle='--', linewidth=0.25)
            plt.plot(conv_time, interpolated_data, '-b')
            plt.title('Low-Pass {0:.3f}Hz / {1}h and input data.'.format(band, 1./band))
            data_axis.set_xlim(x_limits)
            data_axis.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.gcf().autofmt_xdate()
            data_axis.set_xticks(np.arange(x_limits[0], x_limits[1]+1, 5.0))
            data_axis.locator_params(axis='y', nbins=6)
            plt.tight_layout()
            plot.add_plot(width=NoEscape(r'1\textwidth'))
            plot.add_caption('The blue line is the WSPE time serie. The black line indicates the 4th order'
                             ' low-pass Butterworth filter ({0})h. Gaps show NaN values for these'
                             ' dates.'.format(1./band))
            plt.clf()
            plt.close('all')


def apply_low_pass_filter(cutoff_frequency, data, cur_time):
    if np.any(np.isnan(data)):
        logger.info('Interpolating data. Check to use the data as output.')
        data = linear_fill(cur_time, data)
    order = 4
    b, a = signal.butter(order, cutoff_frequency, 'low', output='ba')
    filtered = signal.filtfilt(b, a, data)
    return filtered, data


def compare_u_v_components(doc, hf_dir, hf_spe, buoy_dir, buoy_spe, hf_u, hf_v, hf_converted_time,
                           same_y_limits=False):
    """
    Computes the u v components from two speed and direction variables.
    :param same_y_limits:
    :param doc:
    :param hf_dir:
    :param hf_spe:
    :param buoy_dir:
    :param buoy_spe:
    :param hf_u:
    :param hf_v:
    :param hf_converted_time:
    :return:
    """
    x_limits = [hf_converted_time[0], hf_converted_time[-1]]
    hf_computed_u, hf_computed_v = compute_u_v_components(hf_dir, hf_spe)
    buoy_computed_u, buoy_computed_v = compute_u_v_components(buoy_dir, buoy_spe)
    u_y_lim = [np.nanmin([hf_u, buoy_computed_u]), np.nanmax([hf_u, buoy_computed_u])]
    v_y_lim = [np.nanmin([hf_v, buoy_computed_v]), np.nanmax([hf_v, buoy_computed_v])]
    combined_y_lim = [np.nanmin([u_y_lim[0], v_y_lim[0]]), np.nanmax([u_y_lim[1], v_y_lim[1]])]
    with doc.create(Figure(position='htbp')) as plot:
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(hf_converted_time, hf_computed_u, '-r', label='hf computed')
        axarr[0].plot(hf_converted_time, hf_u, '--k', label='hf read-in')
        axarr[0].plot(hf_converted_time, buoy_computed_u, '--b', label='buoy computed')
        axarr[0].set_title('U component')
        axarr[0].set_xlim(x_limits)
        axarr[0].set_ylim(combined_y_lim)
        axarr[0].set_xticks(np.arange(x_limits[0], x_limits[1]+1, 5.0))
        axarr[0].set_ylabel(r'$ms^{-1}$', rotation=0, horizontalalignment='right')
        axarr[1].plot(hf_converted_time, hf_computed_v, '-r', label='hf computed')
        axarr[1].plot(hf_converted_time, hf_v, '--k', label='hf read-in')
        axarr[1].plot(hf_converted_time, buoy_computed_v, '--b', label='buoy computed')
        axarr[1].xaxis.set_major_formatter(c.settings.xfmt)
        axarr[1].set_ylabel(r'$ms^{-1}$', rotation=0, horizontalalignment='right')
        axarr[1].set_title('V component')
        axarr[1].set_ylim(combined_y_lim)
        if same_y_limits:
            for cur_axis in axarr:
                cur_limits = cur_axis.get_ylim()
                new_limit = np.max(np.abs(cur_limits)) + 0.02
                cur_axis.set_ylim([-new_limit, new_limit])
                # cur_axis.locator_params(axis='y', nbins=6)
        axarr[0].yaxis.set_major_locator(mp.ticker.MaxNLocator(nbins=6, symmetric=True, trim=False))
        # axarr[0].yaxis.set_major_formatter(mp.ticker.ScalarFormatter())
        axarr[1].yaxis.set_major_locator(mp.ticker.MaxNLocator(nbins=6, symmetric=True, trim=False))
        # axarr[1].yaxis.set_major_formatter(mp.ticker.ScalarFormatter())
        f.autofmt_xdate()
        f.suptitle('continuous red line = HFR computed U and V' + "\n" + 'discontinuous black line = HFR output U'
                                  ' and V' + "\n" + 'discontinuous blue line = BUOY output U and V')
        f.set_size_inches(11, 5)
        plt.tight_layout()
        plt.subplots_adjust(top=0.8)
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption('U and V comparisons. The red lines indicates the computed HFR U and V components. The'
                         ' discontinous black lines are the HFR output U and V components. The discontinous blue lines'
                         ' are the Buoy output U and V components. The red and black lines should overlap. This means'
                         ' that the U and V calculation provides the same result.')
        plt.clf()
        plt.close('all')


def compute_u_v_components(direction, spe=None):
    if spe is None:
        spe = np.ones((1, len(direction)))[0]
    u = spe * np.sin(np.deg2rad(direction))
    v = spe * np.cos(np.deg2rad(direction))
    return u, v


def get_title_name(variable):
    try:
        title_name = variable.long_name + ' (' + variable.name + ')'
    except AttributeError:
        logger.info(variable.name + ' has no long name.')
        try:
            title_name = variable.standard_name + ' (' + variable.name + ')'
        except AttributeError:
            logger.info(variable.name + ' has no standard name.')
            title_name = variable.name
    return title_name


def transform_to_full_time(time1, time2):
    if len(time1) > len(time2):
        return time1
    else:
        return time2


def transform_to_full_data(data1, data2, idx):
    if len(data1) > len(data2):
        data1_filled = np.empty((1, len(data1)))[0]
        data1_filled.fill(np.nan)
        data2_filled = np.empty((1, len(data1)))[0]
        data2_filled.fill(np.nan)
        data1_filled = data1
        data2_filled[idx] = data2
    else:
        data1_filled = np.empty((1, len(data2)))[0]
        data1_filled.fill(np.nan)
        data2_filled = np.empty((1, len(data2)))[0]
        data2_filled.fill(np.nan)
        data1_filled[idx] = data1
        data2_filled = data2
    return data1_filled, data2_filled


def get_same_idx(arr1, arr2):
    if len(arr1) > len(arr2):
        return np.in1d(arr1, arr2)
    else:
        return np.in1d(arr2, arr1)


def get_md_datenum(obs_time):
    dates = [datetime.fromtimestamp(ts, tz=pytz.utc) for ts in obs_time]
    return md.date2num(dates)


def plot_tidal_ellipses():
    """
    Here we prepare the tidal ellipses representation...
    :return:
    """
    pass


def get_tidal_ellipse_coordinates(smaj, smin, inc):
    """
    Might be replaced with plot_ellipse (https://casper.berkeley.edu/astrobaki/index.php/Plotting_Ellipses_in_Python)
    Or http://matplotlib.org/api/patches_api.html#matplotlib.patches.Ellipse:
        Ellipse(xy, width, height, angle=0.0, **kwargs)
        xy: center of ellipse
        width: total length (diameter) of horizontal axis
        height: total length (diameter) of vertical axis
        angle: rotation in degrees (anti-clockwise)
    :param smaj: semi major axis width
    :param smin: semi minor axis width
    :param inc: angle in degrees (anti-clockwise)
    :return: arrays with x and y coordinates of the ellipse (density is 100 points)
    """
    th = np.arange(0, 2.*np.pi+0.0001, 2.*np.pi/100)
    inc = inc / 180 * np.pi
    x = smaj * np.cos(th)
    y = smin * np.sin(th)

    (th, r) = cart2pol(x, y)
    (x, y) = pol2cart(th + inc, r)
    return x, y


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y


def spec_rot(u, v):
    """
    Taken from https://github.com/pyoceans/python-oceans/blob/master/oceans/ff_tools/ocfis.py
    Compute the rotary spectra from u,v velocity components
    Parameters
    ----------
    u : array_like
    zonal wind velocity [m s :sup:`-1`]
    v : array_like
    meridional wind velocity [m s :sup:`-1`]
    Returns
    -------
    cw : array_like
    Clockwise spectrum [TODO]
    ccw : array_like
    Counter-clockwise spectrum [TODO]
    puv : array_like
    Cross spectra [TODO]
    quv : array_like
    Quadrature spectra [ TODO]
    Notes
    -----
    The spectral energy at some frequency can be decomposed into two circularly
    polarized constituents, one rotating clockwise and other anti-clockwise.
    Examples
    --------
    TODO: puv, quv, cw, ccw = spec_rot(u, v)
    References
    ----------
    .. [1] J. Gonella Deep Sea Res., 833-846, 1972.
    """
    # Individual components Fourier series.
    fu, fv = list(map(np.fft.fft, (u, v)))
    # Auto-spectra of the scalar components.
    pu = fu * np.conj(fu)
    pv = fv * np.conj(fv)
    # Cross spectra.
    puv = fu.real * fv.real + fu.imag * fv.imag
    # Quadrature spectra.
    quv = -fu.real * fv.imag + fv.real * fu.imag
    # Rotatory components
    # TODO: Check the division, 4 or 8?
    cw = (pu + pv - 2 * quv) / 8.
    ccw = (pu + pv + 2 * quv) / 8.
    n = len(u)
    f = np.arange(0, n) / n
    return puv, quv, cw, ccw, f


def get_inertial_band_hours(latitude):
    return 12. / np.sin(np.deg2rad(latitude))


def selection_main_frequencies(latitude=None):
    freqs = dict()
    freqs['NO1'] = 0.0402686
    freqs['K1'] = 0.0417807
    freqs['J1'] = 0.0432929
    freqs['OO1'] = 0.0448308
    freqs['UPS1'] = 0.0463430
    freqs['EPS2'] = 0.0761773
    freqs['MU2'] = 0.0776895
    freqs['N2'] = 0.0789992
    freqs['M2'] = 0.0805114
    freqs['L2'] = 0.0820236
    freqs['S2'] = 0.0833333
    freqs['ETA2'] = 0.0850736
    freqs['MO3'] = 0.1192421
    freqs['M3'] = 0.1207671
    freqs['MK3'] = 0.1222921
    freqs['SK3'] = 0.1251141
    freqs['MN4'] = 0.1595106
    freqs['M4'] = 0.1610228
    freqs['SN4'] = 0.1623326
    freqs['MS4'] = 0.1638447
    if latitude is not None:
        freqs['inertial'] = 1./get_inertial_band_hours(latitude)
    return freqs


def check_freq_in_range(cur_freq, latitude=None):
    defined_freqs = selection_main_frequencies(latitude)
    cur_identifier = ''
    if defined_freqs['NO1'] < cur_freq < defined_freqs['J1']:
        cur_identifier = 'K1'
    if defined_freqs['K1'] < cur_freq < defined_freqs['OO1']:
        cur_identifier = 'J1'
    if defined_freqs['OO1'] < cur_freq < 0.0468:
        cur_identifier = 'UPS1'
    if 0.075 < cur_freq < defined_freqs['MU2']:
        cur_identifier = 'EPS2'
    if defined_freqs['EPS2'] < cur_freq < defined_freqs['J1']:
        cur_identifier = 'MU2'
    if defined_freqs['MU2'] < cur_freq < defined_freqs['M2']:
        cur_identifier = 'N2'
    if defined_freqs['N2'] < cur_freq < defined_freqs['L2']:
        cur_identifier = 'M2'
    if defined_freqs['M2'] < cur_freq < defined_freqs['S2']:
        cur_identifier = 'L2'
    if defined_freqs['L2'] < cur_freq < defined_freqs['ETA2']:
        cur_identifier = 'S2'
    if defined_freqs['S2'] < cur_freq < 0.086:
        cur_identifier = 'ETA2'
    if 0.118 < cur_freq < defined_freqs['M3']:
        cur_identifier = 'MO3'
    if 'inertial' in defined_freqs:
        if (defined_freqs['inertial'] - 0.002) < cur_freq < (defined_freqs['inertial'] + 0.002):
            cur_identifier = 'inertial'
    return cur_identifier


def get_important_peaks_in_freq_range(freq, power, latitude=None):
    # propably also possible as a one-liner but sorry i am stupid so KISS... keep it simple and stupid
    identifier_storage = []
    main_frequencies_points = []
    remainder_frequencies = []
    remainder_power = []
    # MAX POWER!
    max_power = np.nanmax(power)
    for i in range(0, len(freq)):
        cur_freq = freq[i]
        identifier = check_freq_in_range(cur_freq, latitude)
        if len(identifier_storage) > 0:
            if (identifier != '') & (identifier_storage[-1] != identifier):
                logger.info('freq found: ' + identifier)
                identifier_storage.append(identifier)
                main_frequencies_points.append([cur_freq, power[i]])
            else:
                remainder_frequencies.append(cur_freq)
                remainder_power.append(power[i])
        else:
            if identifier != '':
                logger.info('freq found: ' + identifier)
                identifier_storage.append(identifier)
                main_frequencies_points.append([cur_freq, power[i]])
            else:
                remainder_frequencies.append(cur_freq)
                remainder_power.append(power[i])
    return identifier_storage, main_frequencies_points, np.asarray(remainder_frequencies), np.asarray(remainder_power)


def _datacheck_peakdetect(x_axis, y_axis):
    if x_axis is None:
        x_axis = range(len(y_axis))
    if len(y_axis) != len(x_axis):
        raise ValueError(
                "Input vectors y_axis and x_axis must have same length")
    # needs to be a numpy array
    y_axis = np.array(y_axis)
    x_axis = np.array(x_axis)
    return x_axis, y_axis


def peakdet(v, delta, x=None):
    maxtab = []
    mintab = []

    if x is None:
        x = np.arange(len(v))

    v = np.asarray(v)

    if len(v) != len(x):
        logger.warning('Input vectors v and x must have same length')

    if not np.isscalar(delta):
        logger.warning('Input argument delta must be a scalar')

    if delta <= 0:
        logger.warning('Input argument delta must be positive')

    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN

    lookformax = True

    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)


def peakdetect(y_axis, x_axis=None, lookahead=200, delta=0):
    max_peaks = []
    min_peaks = []
    dump = []   # Used to pop the first hit which almost always is false

    # check input data
    x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis)
    # store data length for later use
    length = len(y_axis)

    # perform some checks
    if lookahead < 1:
        raise ValueError("Lookahead must be '1' or above in value")
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError("delta must be a positive number")

    # maxima and minima candidates are temporarily stored in
    # mx and mn respectively
    mn, mx = np.Inf, -np.Inf

    # Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x

        # look for max
        if y < mx-delta and mx != np.Inf:
            # Maxima peak candidate found
            # look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                max_peaks.append([mxpos, mx])
                dump.append(True)
                # set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf
                if index+lookahead >= length:
                    # end is within lookahead no more peaks can be found
                    break
                continue
            # else:  #slows shit down this does
            #    mx = ahead
            #    mxpos = x_axis[np.where(y_axis[index:index+lookahead]==mx)]

        # look for min
        if y > mn+delta and mn != -np.Inf:
            # Minima peak candidate found
            # look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].min() > mn:
                min_peaks.append([mnpos, mn])
                dump.append(False)
                # set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf
                if index+lookahead >= length:
                    # end is within lookahead no more peaks can be found
                    break
            # else:  #slows shit down this does
            #    mn = ahead
            #    mnpos = x_axis[np.where(y_axis[index:index+lookahead]==mn)]

    # Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            max_peaks.pop(0)
        else:
            min_peaks.pop(0)
        del dump
    except IndexError:
        # no peaks were found, should the function return empty lists?
        pass

    return [max_peaks, min_peaks]


def frequency_plotter(doc, freq, power, significant_freq_names, significant_freq_points, cur_title, month_str, year,
                      other_freqs_x=None, other_freqs_y=None, zoom=False):
    cur_x_lim = [0.025, 0.1]
    if other_freqs_x is not None:
        idx = (other_freqs_x >= 0.04) & (other_freqs_x <= 0.07)
        other_freqs_x = other_freqs_x[idx]
        other_freqs_y = other_freqs_y[idx]
    if zoom:
        cur_title += ' section of interest'
    with doc.create(Figure(position='htbp')) as plot:
        f = plt.figure()
        plt.plot(freq, power, 'k-', lw=0.5)
        ax = plt.gca()
        line_holder = []
        y_lims = ax.get_ylim()
        if zoom:
            texts = []
            for i in range(0, len(significant_freq_points)):
                act_freq = significant_freq_points[i][0]
                act_power = significant_freq_points[i][1]
                plt.plot(act_freq, act_power, 'ro')
                texts.append(plt.text(act_freq, act_power, significant_freq_names[i]))
                cs, = plt.plot([act_freq, act_freq], [y_lims[0], y_lims[1]], '--')
                line_holder.append(cs)
                # cur_line = significant_freq_points[i]
                # cs, = plt.plot(cur_line[0], cur_line[1], '--')
                # line_holder.append(cs)
            for i in range(0, len(other_freqs_x)):
                plt.plot(other_freqs_x[i], other_freqs_y[i], 'ro')
                texts.append(plt.text(other_freqs_x[i], other_freqs_y[i], '{:0.2f}'.format(1./other_freqs_x[i])))
        plt.legend(line_holder, significant_freq_names, numpoints=1, fontsize=10)
        plt.title(cur_title)
        ax.set_ylabel('${|dft|}^{2}$', rotation=0, horizontalalignment='right')
        ax.set_xlabel('cycles/hour')
        f.set_size_inches(11, 3)
        if zoom:
            y_min = 0
            selection_idx = np.logical_and(freq >= 0.004, freq <= 0.1)
            y_max = np.nanmax(power[selection_idx])
            plt.ylim([y_min, y_max])
            plt.xlim(cur_x_lim)
            adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), force_points=0.9,
                        expand_points=(1.2, 1.3))
        plt.tight_layout()
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption(NoEscape(cur_title + ' in ' + month_str + ' ' + str(year)))
        plt.clf()
        plt.close('all')
    if not zoom:
        frequency_plotter(doc, freq, power, significant_freq_names, significant_freq_points, cur_title, month_str, year,
                          other_freqs_x=other_freqs_x, other_freqs_y=other_freqs_y, zoom=True)


def compute_dft_spectrum(obs_time, component_data):
    filled_component_data = linear_fill(obs_time, component_data)

    x = filled_component_data
    fs = 1./3600.
    t = np.arange(0, 0.25, fs/2.)
    n = len(t)
    y = fft(x, n)
    power = np.power(abs(y[0:(n/2)]), 2)
    nyquist = 1/2.
    freq = np.arange(0, n/2)/(n/2.)*nyquist

    # max_peaks, min_peaks = peakdetect(power, lookahead=3)
    max_peaks, min_peaks = peakdet(power, 5)
    x_p, y_p = zip(*max_peaks)
    peak_freqs = freq[list(map(int, x_p))]
    peak_power = power[list(map(int, x_p))]
    return freq, power, peak_freqs, peak_power


def plot_energy_spectrum(doc, obs_time, u, v, month_str, year, latitude=None):
    cur_title = r'U component ${|dft|}^{2}$'
    freq, power, peak_freqs, peak_power = compute_dft_spectrum(obs_time, u)
    significant_freq_names, significant_freq_points, remainder_freqs, remainder_power = \
        get_important_peaks_in_freq_range(peak_freqs, peak_power, latitude)
    frequency_plotter(doc, freq, power, significant_freq_names, significant_freq_points, cur_title, month_str, year,
                      other_freqs_x=remainder_freqs, other_freqs_y=remainder_power, zoom=True)

    found_freqs = []
    for stored_freq in significant_freq_points:
        found_freqs.append(stored_freq[0])
    write_frequencies_table(doc, 'U component', significant_freq_names, found_freqs)

    cur_title = r'V component ${|dft|}^{2}$'
    freq, power, peak_freqs, peak_power = compute_dft_spectrum(obs_time, v)
    significant_freq_names, significant_freq_points, remainder_freqs, remainder_power = \
        get_important_peaks_in_freq_range(peak_freqs, peak_power, latitude)
    frequency_plotter(doc, freq, power, significant_freq_names, significant_freq_points, cur_title, month_str, year,
                      other_freqs_x=remainder_freqs, other_freqs_y=remainder_power, zoom=True)
    found_freqs = []
    for stored_freq in significant_freq_points:
        found_freqs.append(stored_freq[0])
    write_frequencies_table(doc, 'V component', significant_freq_names, found_freqs)


def write_frequencies_table(doc, cur_title, identifiers, frequencies):
    with doc.create(Table(position='htb')) as t:
        doc.append(NoEscape(r'\begin{center}'))
        with doc.create(Tabular('|c|c|c|', pos='htb')) as table:
            table.add_hline()
            table.add_row(((MultiColumn(3, align='|c|', data=cur_title)),))
            table.add_row(('Identifier', 'cycles/hour', 'hour'))
            table.add_hline()
            for i in range(0, len(identifiers)):
                table.add_row((identifiers[i], np.round(frequencies[i], 4), np.round(1./frequencies[i], 2)))
                table.add_hline()
        t.add_caption('Identified frequencies for ' + cur_title)
        doc.append(NoEscape(r'\end{center}'))


def linear_fill(x, y, kind=None):
    """
    Fill the gap in a time serie using linear interpolation.
    Taken from https://github.com/ctroupin/SOCIB_plots/blob/master/HFradar/energy_spectrum_radar.ipynb all credit goes
    to Charles Troupin, Socib.
    :param x: 1-D time array
    :param y: 1-D data array
    :param kind: str
    :return: filled gaps (NaN values) using one-dimensional linear interpolation (numpy.interp)
    """
    if kind is None:
        kind = 'linear'
    good_values = np.where(~np.isnan(y))
    missing_values = np.where(np.isnan(y))
    y_interp = np.copy(y)
    if kind == 'linear':
        y_interp[missing_values] = np.interp(x[missing_values], x[good_values], y[good_values])
    # elif kind == 'spline':
    #     f = interpolate.interp1d(x, y)
    #     y_interp[missing_values] = np.interp(x[missing_values], x[good_values], y[good_values])
    return y_interp


def get_good_data_only():
    # TODO implement
    pass


def plot_quiver_direction_overlapping(doc, cur_time, lower_direction, plot_title, upper_direction, upper_amplifier=None,
                                      s_name=None, x_limits=None, input_month_title=None, lower_amplifier=None,
                                      shared_qc_idx_upper=None, shared_qc_idx_lower=None):
    # TODO: jajaja i c, very dirty. but had to be done quick, so remove duplicate code bases when there is some time
    if x_limits is None:
        x_limits = [np.nanmin(cur_time) - 1, np.nanmax(cur_time) + 1]
    else:
        x_limits = [x_limits[0] - 1, x_limits[1] + 1]
    plt.rcParams.update({'font.size': 13})
    try:
        month_title = c.settings.month_str + ' ' + str(c.settings.year)
    except AttributeError:
        if input_month_title is None:
            logger.warning('Trying to access non-set setting variable.', exc_info=True)
            month_title = ''
        else:
            month_title = input_month_title
    if s_name is None:
        my_title = 'Evolution of ' + plot_title + ' in ' + month_title + '.'
    else:
        my_title = 'Evolution of ' + plot_title + ' at ' + s_name + ' in ' + month_title + '.'
    if np.all(np.isnan(lower_direction)) or np.all(np.isnan(upper_direction)):
        print 'Only NaNs found at ' + plot_title
        doc.append('Only NaNs encountered at one of the plots from ' + plot_title)
        return
    if upper_amplifier is not None:
        upper_amplifier = normalize_data(upper_amplifier, 'mean')
    if lower_amplifier is not None:
        lower_amplifier = normalize_data(lower_amplifier, 'mean')

    u, v = compute_u_v_components(lower_direction, lower_amplifier)
    upper_u, upper_v = compute_u_v_components(upper_direction, upper_amplifier)
    # u_raw = np.cos(np.deg2rad(lower_direction))
    # v_raw = np.sin(np.deg2rad(lower_direction))
    # u_raw_lower = np.cos(np.deg2rad(upper_direction))
    # v_raw_lower = np.sin(np.deg2rad(upper_direction))
    if (upper_amplifier is not None) or (lower_amplifier is not None):
        plot_title += ' normalized and scaled with speed'
    if lower_amplifier is None:
        vector_length = np.sqrt((u * u) + (v * v))
        u = u / vector_length
        v = v / vector_length
    # else:
    #     upper_amplifier = normalize_data(upper_amplifier, 'mean')
    #     u = -upper_amplifier * u_raw
    #     v = -upper_amplifier * v_raw
    if upper_amplifier is None:
        upper_vector_length = np.sqrt((upper_u * upper_u) + (upper_v * upper_v))
        upper_u = upper_u / upper_vector_length
        upper_v = upper_u / upper_vector_length
    # else:
    #     lower_amplifier = normalize_data(lower_amplifier, 'mean')
    #     lower_u = -lower_amplifier * u_raw_lower
    #     lower_v = -lower_amplifier * v_raw_lower

    with doc.create(Figure(position='htbp')) as plot:
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].set_title("\n".join(wrap(plot_title, 70)))
        try:
            if shared_qc_idx_upper is not None:
                upper_u[~shared_qc_idx_upper] = np.nan
                upper_v[~shared_qc_idx_upper] = np.nan
                axarr[0].quiver(cur_time, 0, upper_u, upper_v, color='r', alpha=0.5,
                                width=0.004, units='width', scale=1, scale_units='y', headlength=0, headwidth=1)
            else:
                axarr[0].quiver(cur_time, 0, upper_u, upper_v, color='r', alpha=0.5, width=0.004, units='width',
                                scale=1, scale_units='y', headlength=0, headwidth=1)
            axarr[0].set_xlim(x_limits)
            axarr[0].set_ylim(-1, 1)
            axarr[0].xaxis.set_major_formatter(c.settings.xfmt)
            axarr[0].set_xticks(np.arange(x_limits[0]+1, x_limits[1], 5.0))
            if shared_qc_idx_lower is not None:
                u[~shared_qc_idx_lower] = np.nan
                v[~shared_qc_idx_lower] = np.nan
                axarr[1].quiver(cur_time, 0, u, v, color='r', alpha=0.5, width=0.004,
                                units='width', scale=1, scale_units='y', headlength=0, headwidth=1)
            else:
                axarr[1].quiver(cur_time, 0, u, v, color='r', alpha=0.5, width=0.004, units='width', scale=1,
                                scale_units='y', headlength=0, headwidth=1)
            axarr[1].set_xlim(x_limits)
            axarr[1].set_ylim(-1, 1)
            #axarr[1].set_xticks([])
            f.autofmt_xdate()
            f.set_size_inches(11, 5)
            plt.tight_layout()
        except RuntimeWarning:
            logger.debug('Problem at plotting quiver.', exc_info=True)
        plt.subplots_adjust(top=0.8)
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption(my_title)
        plt.clf()
        plt.close('all')


def write_header_tabular(table, variable, spec_len):
    table.add_hline()
    table.add_row(((MultiColumn(spec_len, align='|c|',
                                data=get_standard_name(variable) + ' (' + variable.name + ')')),))


def write_data_only_tabular(table, variable, cur_mean, cur_min, cur_max, cur_min_time, cur_max_time):
    cur_units = variable.units
    table.add_row((MultiColumn(2, align='|c|', data=''), MultiColumn(3, align='c|', data='Data Statistics')))
    table.add_row((MultiColumn(1, align='|c', data=''), '', 'Mean (' + cur_units + ')',
                   'Min (' + cur_units + ')', 'Max (' + cur_units + ')'))
    table.add_hline()
    table.add_row((MultiRow(2, data=''), 'Value', "%.2f" % cur_mean, "%.2f" % cur_min, "%.2f" % cur_max))
    table.add_row(('', 'Time', '-', cur_min_time, cur_max_time))


def write_qc_only_tabular(table, sum_good, sum_prob_good, sum_prob_bad, sum_bad, sum_spike, sum_nan, percent_good,
                          percent_rest, percent_nan):
    table.add_row((MultiColumn(2, align='|c|', data=''), MultiColumn(6, align='c|', data='QC Flags')))
    table.add_row((MultiColumn(1, align='|c', data=''), '', 1, 2, 3, 4, 6, 9))
    table.add_hline()
    table.add_row((MultiRow(2, data=''), 'N', sum_good, sum_prob_good, sum_prob_bad, sum_bad, sum_spike, sum_nan))
    table.add_row(('', '%', percent_good, MultiColumn(4, align='c|', data=percent_rest), percent_nan))


def write_data_and_qc_tabular(table, variable, sum_good, sum_prob_good, sum_prob_bad, sum_bad, sum_spike, sum_nan,
                              cur_mean, cur_min, cur_max, percent_good, percent_rest, percent_nan, cur_min_time,
                              cur_max_time):
    cur_units = variable.units
    table.add_row((MultiColumn(2, align='|c|', data=''), MultiColumn(6, align='c|', data='QC Flags'),
                   MultiColumn(4, align='c|', data='Data Statistics')))
    table.add_row((
        MultiColumn(1, align='|c', data=''), '', 1, 2, 3, 4, 6, 9, MultiColumn(1, align='c', data=''),
        'Mean (' + cur_units + ')', 'Min (' + cur_units + ')',
        'Max (' + cur_units + ')'))
    table.add_hline()
    table.add_row((MultiRow(2, data=''), 'N', sum_good, sum_prob_good, sum_prob_bad, sum_bad, sum_spike,
                   sum_nan, 'Value', "%.2f" % cur_mean, "%.2f" % cur_min, "%.2f" % cur_max))
    table.add_row(('', '%', percent_good, MultiColumn(4, align='c|', data=percent_rest), percent_nan,
                   'Time', '-', cur_min_time, cur_max_time))


def write_specific_statistics_tabular(table, v, v_mean, v_std, v_min, v_max, v_good_percent):
    cur_units = v.units
    table.add_row(('Mean (%s)' % cur_units, 'Std (%s)' % cur_units, 'Min (%s)' % cur_units, 'Max (%s)' % cur_units,
                   '% Good'))
    table.add_hline()
    table.add_row(('%.2f' % v_mean, '%.2f' % v_std, '%.2f' % v_min, '%.2f' % v_max,
                   '%.2f' % v_good_percent))


def write_table(doc, tab_spec, variable, month_str, year, **kwargs):
    with doc.create(Table(position='htb')) as t:
        doc.append(NoEscape(r'\begin{center}'))
        with doc.create(Tabular(tab_spec, pos='htb')) as table:
            write_header_tabular(table, variable, tab_spec.count('c'))
            if tab_spec == '|cc|ccc|':
                write_data_only_tabular(table, variable, kwargs['cur_mean'], kwargs['cur_min'], kwargs['cur_max'],
                                        kwargs['cur_min_time'], kwargs['cur_max_time'])
            elif tab_spec == '|cc|c|cccc|c|':
                write_qc_only_tabular(table, kwargs['sum_good'], kwargs['sum_prob_good'],
                                      kwargs['sum_prob_bad'], kwargs['sum_bad'], kwargs['sum_spike'], kwargs['sum_nan'],
                                      kwargs['percent_good'], kwargs['percent_rest'], kwargs['percent_nan'])
            elif tab_spec == '|cc|c|cccc|c|c|ccc|':
                write_data_and_qc_tabular(table, variable, kwargs['sum_good'], kwargs['sum_prob_good'],
                                          kwargs['sum_prob_bad'], kwargs['sum_bad'], kwargs['sum_spike'],
                                          kwargs['sum_nan'], kwargs['cur_mean'], kwargs['cur_min'], kwargs['cur_max'],
                                          kwargs['percent_good'], kwargs['percent_rest'], kwargs['percent_nan'],
                                          kwargs['cur_min_time'], kwargs['cur_max_time'])
            elif tab_spec == '|c|c|c|c|c|':
                write_specific_statistics_tabular(table, variable, kwargs['v_mean'], kwargs['v_std'], kwargs['v_min'],
                                                  kwargs['v_max'], kwargs['v_good_percent'])
            else:
                logger.warning('Undefined tabular specification.')
            table.add_hline()
        t.add_caption('Data Summary for ' + get_standard_name(variable) + ' in ' + month_str + ' ' + str(year) + '.')
        doc.append(NoEscape(r'\end{center}'))


def get_time_str_from_converted_time_idx(idx, converted_time):
    return md.num2date(converted_time[idx]).strftime("%Y-%m-%d %H:%M:%S")


def get_mean_min_max(data, time_ref=None):
    if time_ref is None:
        # Makes only sense for timeless data.
        pass
    else:
        cur_mean = np.nanmean(data)
        mean_time_str = ''
        cur_min = np.nanmin(data)
        cur_min_idx = np.where(data == cur_min)
        min_time_str = get_time_str_from_converted_time_idx(cur_min_idx[0][0], time_ref)
        cur_max = np.nanmax(data)
        cur_max_idx = np.where(data == cur_max)
        max_time_str = get_time_str_from_converted_time_idx(cur_max_idx[0][0], time_ref)
        return cur_mean, mean_time_str, cur_min, min_time_str, cur_max, max_time_str


def plot_histogram_radial_files(doc, stations_end_time, stations_bins):

    with doc.create(Figure(position='htbp')) as plot:
        plt.figure(figsize=(11,5))
        galf_bar = plt.bar(stations_end_time[0:len(stations_end_time)/2], stations_bins['GALF'], color='cornflowerblue',
                           align='center', width=8)
        form_bar = plt.bar(stations_end_time[0:len(stations_end_time)/2], stations_bins['FORM'], color='lightcoral',
                           align='center', bottom=stations_bins['GALF'], width=8)
        ax = plt.gca()
        ax.set_xticks(stations_end_time[0:len(stations_end_time)/2])
        x_labels = []
        for dt in stations_end_time[0:len(stations_end_time)/2]:
            x_labels.append(str(dt.day) + '.' + str(dt.month) + '.' + str(dt.year))
        ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='right', rotation_mode='anchor')
        ax.tick_params(direction='out', axis='y', pad=15)
        ax.set_ylim([0, 600])
        plt.title('Number of radial files per 10 days.')
        plt.xticks(stations_end_time[0:len(stations_end_time)/2])
        plt.legend((galf_bar, form_bar), ('GALF', 'FORM'))
        plt.tight_layout()
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption('Number of radial files')
        plt.clf()
        plt.close('all')


def plot_threshold_graph(doc, x_data, y_data, lower_threshold, upper_threshold, variable, year, month_str,
                         x_limits=None):
    if x_limits is None:
        x_limits = [x_data[0], x_data[-1]]
    try:
        title_name = variable.long_name + ' (' + variable.name + ')'
    except AttributeError:
        logger.info(variable.name + ' has no long name.')
        try:
            title_name = variable.standard_name + ' (' + variable.name + ')'
        except AttributeError:
            logger.info(variable.name + ' has no standard name.')
            title_name = variable.name
    with doc.create(Figure(position='htbp')) as plot:
        plt.figure()
        plt.plot(x_data, y_data, 'k', linewidth=1, )
        if lower_threshold is not None:
            lower_line = np.ones((1, len(x_data)))[0] * lower_threshold
            plt.fill_between(x_data, lower_line, y_data, facecolor='gray', alpha=0.5)
            # plt.plot(x_data, lower_line, '-r', linewidth=1, )
        if upper_threshold is not None:
            upper_line = np.ones((1, len(x_data)))[0] * upper_threshold
            # plt.plot(x_data, upper_line, '-r', linewidth=1, )
            plt.fill_between(x_data, y_data, upper_line, facecolor='gray', alpha=0.5)
        data_axis = plt.gca()
        data_axis.xaxis.set_major_formatter(c.settings.xfmt)
        if variable.units == '1':
            cur_units = ''
        else:
            cur_units = variable.units
        data_axis.set_ylabel(cur_units, rotation=0, horizontalalignment='right')
        data_axis.grid(b=False, which='major', color='k', linestyle='--', linewidth=0.25)
        # data_axis.set_ylim([data_axis.get_ylim()[0] - 5, data_axis.get_ylim()[1] + 5])
        plt.margins(y=0.1)
        data_axis.set_xlim(x_limits)
        data_axis.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.title('\n'.join(wrap(title_name, 70)))
        plt.gcf().autofmt_xdate()
        plt.gcf().set_size_inches(11, 3)
        data_axis.set_xticks(np.arange(x_limits[0], x_limits[1]+1, 5.0))
        data_axis.locator_params(axis='y', nbins=6)
        plt.tight_layout()
        plt.subplots_adjust(top=0.8)
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption('Evolution of ' + title_name + ' and indicating valid thresholds in ' + month_str + ' ' +
                         str(year) + '.')
        plt.clf()
        plt.close('all')


def get_hf_simple_statistics(data, data_qc=None):
    data_mean = np.nanmean(data)
    data_std = np.nanstd(data)
    data_min = np.nanmin(data)
    data_max = np.nanmax(data)
    if data_qc is None:
        good_percent = None
    else:
        good_idx = data_qc == 1
        good_amount = np.where(good_idx == 1)[0]
        good_percent = float(len(good_amount))/len(data_qc)*100
    return data_mean, data_std, data_min, data_max, good_percent


def get_hf_radial_site_file_size(dir_path, elements, threshold):
    file_counter = 0
    good_counter = 0
    for e in elements:
        if e[0:4] == 'RDLi':
            continue
        file_counter += 1
        cur_path = dir_path + '' + e
        cur_size = os.stat(cur_path)
        cur_size = cur_size.st_size/1024.0
        if cur_size > threshold:
            good_counter += 1
    good_percent = float(good_counter)/file_counter*100
    return good_percent, file_counter, good_counter


def get_hf_radial_sites_file_sizes(galf_path, galf_elements, form_path, form_elements, totals_path, totals_elements,
                                   galf_threshold, form_threshold, total_threshold):
    galf_good_percent, _, _ = get_hf_radial_site_file_size(galf_path, galf_elements, galf_threshold)
    form_good_percent, _, _ = get_hf_radial_site_file_size(form_path, form_elements, form_threshold)
    totals_good_percent, _, _ = get_hf_radial_site_file_size(totals_path, totals_elements, total_threshold)
    logger.info('Galfi good percent: ' + str(galf_good_percent))
    logger.info('Formentera good percent: ' + str(form_good_percent))
    logger.info('Total good percent: ' + str(totals_good_percent))
    return galf_good_percent, form_good_percent, totals_good_percent


def plot_temporal_and_spatial_availability(doc, temporal_avail_percent, non_zero_sorted_spatial_avail, temp_avail):
    red_x = np.arange(0, 100, 0.01)
    red_y = np.zeros((1, 10000))[0]
    red_y[0:8000] = 80
    red_x[8001] = red_x[7999]
    cur_x = temporal_avail_percent
    cur_x = np.append(cur_x, [100.])
    cur_y = non_zero_sorted_spatial_avail
    cur_y = np.append(cur_y, [0.])
    f = interpolate.interp1d(cur_x, cur_y)
    interpolated_temporal_data = np.arange(0, 100, 0.01)
    interpolated_spatial_data = f(interpolated_temporal_data)
    idx = interpolated_temporal_data <= 80
    with doc.create(Figure(position='htbp')) as plot:
        plt.figure()
        plt.plot(interpolated_temporal_data, interpolated_spatial_data)
        plt.plot(red_x, red_y, color='red')
        plt.fill_between(interpolated_temporal_data, interpolated_spatial_data, red_y,
                         where=red_y >= interpolated_spatial_data, facecolor='salmon', interpolate=True)
        data_axis = plt.gca()
        data_axis.set_ylabel('%', rotation=0, horizontalalignment='right')
        data_axis.set_xlabel('%')
        data_axis.grid(b=False, which='major', color='k', linestyle='--', linewidth=0.25)
        data_axis.get_yaxis().get_major_formatter().set_useOffset(False)
        cur_title = 'Spatial vs Temporal Availability'
        if np.all(interpolated_spatial_data[idx]>=80):
            cur_title += '. The MARACOOS metric for long range data coverage and availability has surpassed the 80% coverage 80% of the time during this period.'
        plt.title('\n'.join(wrap(cur_title, 90)))
        plt.gcf().set_size_inches(11, 3)
        data_axis.locator_params(axis='y', nbins=6)
        plt.tight_layout()
        plt.subplots_adjust(top=0.8)
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption('Spatial (x-axis) vs Temporal (y-axis) Availability')
        plt.clf()
        plt.close('all')


def plot_spatial_availability(doc, lon2, lat2, masked_array, basemap1):
    mycmap = mp.cm.Greens
    mycmap.set_bad('w', 0.0)
    lon2, lat2 = basemap1(lon2, lat2)

    with doc.create(Figure(position='htbp')) as plot:
        plt.figure()
        basemap1.drawcoastlines(linewidth=.25, zorder=3)
        basemap1.drawmapboundary(fill_color='white')
        basemap1.drawparallels(np.arange(30., 40., 1.), labels=[True, False, True, False], zorder=2)
        basemap1.drawmeridians(np.arange(-10., 2., 1.), labels=[False, True, False, True], zorder=2)
        q = basemap1.pcolor(lon2, lat2, masked_array, cmap=mycmap, zorder=4)
        plt.clim(0, 100)
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="5%", pad=0.05)

        bounds = np.linspace(0,100,11)
        norm = mp.colors.BoundaryNorm(bounds, mycmap.N)

        # cb = plt.colorbar(q, cax=cax)
        cb = mp.colorbar.ColorbarBase(cax, cmap=mycmap, norm=norm, spacing='proportional', ticks=bounds,
                                      boundaries=bounds, format='%1i')
        cb.set_label(r'%', rotation=0, labelpad=15)
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption('Gridded Availabilty of Data')
        plt.clf()
        plt.close('all')


def plot_hf_temporal_availability(doc, form_path, galf_path, hf_time, year, month):
    form_elements = sorted(os.listdir(form_path))
    galf_elements = sorted(os.listdir(galf_path))
    logger.info('Formentera ' + str(len(form_elements)) + ' entries. Galfi ' + str(len(galf_elements)) + '.')
    pattern_search = str(year) + '_' + str(month).zfill(2) + '_'
    last_day = 1
    last_hour = 0
    form_missing_times = get_missing_times(form_elements, pattern_search, last_day, last_hour, year, month)
    galf_missing_times = get_missing_times(galf_elements, pattern_search, last_day, last_hour, year, month)
    non_nan_hf_time_form = np.asarray(md.date2num([datetime.fromtimestamp(ts, tz=pytz.utc) for ts in np.copy(hf_time)]))
    non_nan_hf_time_galf = np.asarray(md.date2num([datetime.fromtimestamp(ts, tz=pytz.utc) for ts in np.copy(hf_time)]))

    for m_t in form_missing_times:
        if m_t in hf_time:
            temp_idx = np.where(hf_time == m_t)[0]
            non_nan_hf_time_form[temp_idx] = np.nan

    ibiz_nan_idx = np.isnan(non_nan_hf_time_form)

    for m_t in galf_missing_times:
        if m_t in hf_time:
            temp_idx = np.where(hf_time == m_t)[0]
            non_nan_hf_time_galf[temp_idx] = np.nan

    data_form = np.ones((len(non_nan_hf_time_form), 1)) - 0.5
    data_galf = np.ones((len(non_nan_hf_time_galf), 1))

    plot_horizontal_fat_lines(doc, non_nan_hf_time_form, data_form, non_nan_hf_time_galf, data_galf)
    return np.logical_or(ibiz_nan_idx, np.isnan(non_nan_hf_time_galf))


def plot_horizontal_fat_lines(doc, time_one, data_one, time_two, data_two, x_limits=None):
    if x_limits is None:
        if np.nanmin(time_one) < np.nanmin(time_two):
            start_time = np.nanmin(time_one)
        else:
            start_time = np.nanmin(time_two)
        if np.nanmax(time_one) > np.nanmax(time_two):
            end_time = np.nanmax(time_one)
        else:
            end_time = np.nanmax(time_two)
        x_limits = [start_time, end_time]
    with doc.create(Figure(position='htbp')) as plot:
        f = plt.figure()
        plt.plot(time_one, data_one, color='k', linewidth=4.0)
        plt.plot(time_two, data_two, color='g', linewidth=4.0)
        data_axis = plt.gca()
        data_axis.xaxis.set_major_formatter(c.settings.xfmt)
        data_axis.set_ylabel('', rotation=0, horizontalalignment='right')
        data_axis.set_xlim(x_limits)
        plt.title('\n'.join(wrap('Temporal availability', 70)))
        f.set_size_inches(11, 3)
        f.autofmt_xdate()
        plt.yticks([0.5, 1], ['Form', 'Galf'])
        plt.ylim([0, 1.5])
        plt.grid()
        data_axis.set_xticks(np.arange(x_limits[0], x_limits[1] + 1, 5.0))
        plt.tight_layout()
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption('Form and Galf temporal availability')
        plt.clf()
        plt.close('all')


def totimestamp(dt, epoch=datetime(1970, 1, 1)):
    td = dt - epoch
    return (td.microseconds + (td.seconds + td.days * 86400) * 10**6) / 10**6


def get_missing_times(elements, pattern, last_day, last_hour, year, month):
    missing_times = []
    for cur_file in elements:
        if last_hour == 23:
            last_hour = 0
            last_day += 1
        else:
            last_hour += 1
        temp_idx = cur_file.find(pattern)
        if temp_idx == -1:
            # fix invalid file name
            logger.info('Invalid file found {0}'.format(cur_file))
            continue
        cur_day_hour = cur_file[temp_idx + len(pattern):-4]
        cur_day = int(cur_day_hour[0:2])
        cur_hour = int(cur_day_hour[3:-2])
        if (cur_hour != last_hour) and (cur_day != last_day):
            day_difference = cur_day - last_day
            hour_difference = day_difference * 24 + cur_hour - last_hour
            missing_times.extend(create_missing_times(year, month, last_day, last_hour, hour_difference))
        last_day = cur_day
        last_hour = cur_hour
    return np.unique(missing_times)


def create_missing_times(year, month, start_day, start_hour, hour_difference):
    temp_times = []
    for i in range(1, hour_difference + 1):
        if start_hour == 24:
            start_day += 1
            start_hour = 0
        temp_times.append(totimestamp(datetime(year, month, start_day, start_hour)))
        start_hour += 1
    return temp_times


def compare_angles(data1, data2):
    return (data2 - data1 + 180) % 360 - 180


def plot_overlapping_1d_graphs(doc, cur_time, data_one, data_two, year, month_str, x_limits=None, title_str=None,
                               y_label=None, is_degree_0_360_y_limit=False):
    if x_limits is None:
        x_limits = [cur_time[0], cur_time[-1]]
    if y_label is None:
        y_label == ''
    elif y_label == '1':
        y_label = ''
    if title_str is None:
        title_str = ''
    with doc.create(Figure(position='htbp')) as plot:
        plt.figure()
        plt.plot(cur_time, data_one, color='k', linewidth=0.5)
        plt.plot(cur_time, data_two, color='b', linewidth=0.5)
        data_axis = plt.gca()
        data_axis.xaxis.set_major_formatter(c.settings.xfmt)
        y_label = transform_y_label(y_label)
        data_axis.set_ylabel(y_label, rotation=0, horizontalalignment='right')
        data_axis.grid(b=False, which='major', color='k', linestyle='--', linewidth=0.25)
        data_axis.set_xlim(x_limits)
        data_axis.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.title('\n'.join(wrap(title_str, 80)))
        plt.gcf().autofmt_xdate()
        plt.gcf().set_size_inches(11, 3)
        data_axis.set_xticks(np.arange(x_limits[0], x_limits[1]+1, 5.0))
        has_set_ticks = False
        if is_degree_0_360_y_limit:
            data_axis.set_ylim([0, 360])
            data_axis.set_yticks([0, 90, 180, 270, 360])
            has_set_ticks = True
        if not has_set_ticks:
            data_axis.yaxis.set_major_locator(mp.ticker.MaxNLocator(nbins=6, trim=False))
            # data_axis.locator_params(axis='y', nbins=6)
        plt.tight_layout()
        plt.subplots_adjust(top=0.8)
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption(title_str + ' in ' + month_str + ' ' + str(year) + '.')
        plt.clf()
        plt.close('all')


def get_idx_closest_grid_point(lat_1d_vector, lon_1d_vector, reference_lat, reference_lon):
    return np.searchsorted(lat_1d_vector, reference_lat), np.searchsorted(lon_1d_vector, reference_lon) - 1


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mp.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def hf_monthly_mean_direction_plot(doc, u_mean, v_mean, np_longrid, np_latgrid, basemap, buoy_lat, buoy_lon):
    # Could be refined with more parameters
    uv_mean_norm = (u_mean * u_mean + v_mean * v_mean) ** 0.5
    with doc.create(Figure(position='htbp')) as plot:
        plt.figure()
        longrid, latgrid = basemap(np_longrid, np_latgrid)
        basemap.drawcoastlines(linewidth=.25, zorder=4)
        basemap.drawmapboundary(fill_color='white')
        basemap.drawparallels(np.arange(30., 40., 1.), labels=[True, False, True, False], zorder=2)
        basemap.drawmeridians(np.arange(-10., 2., 1.), labels=[False, True, False, True], zorder=2)
        cmap = plt.cm.winter_r
        # cmap = truncate_colormap(cmap, 0.2, 1.0)
        q = basemap.quiver(longrid, latgrid, u_mean, v_mean, uv_mean_norm, cmap=cmap)
        buoy_x, buoy_y = basemap(buoy_lon, buoy_lat)
        buoy_point = basemap.plot(buoy_x, buoy_y, 'bs', markersize=6)
        plt.legend(buoy_point, ['Ibiza Channel Buoy'], loc='upper left', numpoints=1)
        plt.clim(0, 0.4)
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = plt.colorbar(q, cax=cax, format='%0.2f', extend='max')
        # cb.set_label(r'in $\frac{m}{s}$', rotation=0)
        cb.set_label(r'm/s', rotation=0, labelpad=15)
        tick_locator = ticker.MaxNLocator(nbins=8)
        cb.locator = tick_locator
        cb.ax.yaxis.set_major_locator(ticker.AutoLocator())
        cb.update_ticks()
        plt.quiverkey(q, 0.1, 0.05, 0.2, r'$0.2 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption('Monthly averaged velocity vectors.')
        plt.clf()
        plt.close('all')


def hf_plot_grid_and_buoy_position(np_langrid, np_longrid, buoy_lat, buoy_lon, lat_grid, lon_grid):
    m = get_basemap('h')
    x, y = m(np_longrid, np_langrid)
    buoy_x, buoy_y = m(buoy_lon, buoy_lat)
    grid_x, grid_y = m(lon_grid, lat_grid)
    m.drawcoastlines(linewidth=.25, zorder=4)
    m.drawmapboundary(fill_color='white')
    m.drawparallels(np.arange(30., 40., 1.), labels=[True, False, True, False], zorder=2)
    m.drawmeridians(np.arange(-10., 2., 1.), labels=[False, True, False, True], zorder=2)
    m.scatter(x, y, marker='o', color='k', facecolors='none', s=3)
    m.scatter(buoy_x, buoy_y, marker='o', color='r', facecolors='none', s=10)
    m.scatter(grid_x, grid_y, marker='s', color='r', facecolors='none', s=20)


def get_basemap(_resolution):
    logger.info('Creating basemap...')
    return Basemap(projection='merc', llcrnrlat=38.30, urcrnrlat=39.50, llcrnrlon=-0.35, urcrnrlon=1.80, lat_ts=35.,
                   resolution=_resolution)


def get_temporal_mean_from_grid(variable):
    try:
        temp_nanmean = np.nanmean(variable[0::, :, :], axis=0)
    except RuntimeWarning:
        logger.debug('Warning using nanmean.', exc_info=True)
    finally:
        return temp_nanmean


def get_differences(data1, data2, time1=None, time2=None):
    """
    Will compute the differences between two data sets. Will return np.NaN values if one entry is np.NaN in one of the
    data sets. Proposed interpolation is intended to match eventual time mismatches. Requires numpy arrays for this.
    Philosophy: data2 - data1 >> return
    :param data1:
    :param data2:
    :param time1:
    :param time2:
    :param interpolation:
    :return:
    """
    # TODO: implement different times support
    if time1 is None or time2 is None:
        logger.debug('get_differences time not set.')
        time1 = None
        time2 = None
    else:
        same_idx = get_same_idx(time1, time2)
        data1, data2 = get_data_from_same_idx(same_idx, data1, data2)
    return data2 - data1


def get_data_from_same_idx(same_idx, data1, data2):
    if len(same_idx) == len(data1):
        return data1[same_idx], data2
    elif len(same_idx) == len(data2):
        return data1, data2[same_idx]
    else:
        logger.error('No dimension matches same_idx.')


def get_thredds_opendap_link(folder, sub_folder, processing_level, prefix, station_name, year, month):
    """

    :param station_name: string
    :param folder: string
    :param sub_folder: string
    :param processing_level: int
    :param prefix: string
    :param year: int
    :param month: int
    :return: string

    example:
        get_thredds_opendap_link('hf_radar', 'hf_radar_ibiza-scb_codarssproc001', 1, 'dep0001',
        'hf-radar-ibiza-scb_codarssproc001', 2016, 4)
    """
    str_processing_level = 'L' + str(processing_level)
    link = 'http://thredds.socib.es/thredds/dodsC/' + folder + '/' + sub_folder + '/' + str_processing_level + '/' + \
           str(year) + '/' + prefix + '_' + station_name + '_' + str_processing_level + '_' + str(year) + '-' + \
           str(month).zfill(2) + '.nc'
    return link


def create_netcdf_link(jwebchart_link):
    start_idx = jwebchart_link.find('?file=')
    start_idx += 6
    end_idx = jwebchart_link.find('/L1/')
    end_idx += 4
    appendix = jwebchart_link[end_idx::]
    link = jwebchart_link[start_idx:end_idx] + str(c.settings.year) + '/'
    rep_idx = appendix.find('_latest.nc')
    appendix = appendix[0:rep_idx] + '_' + str(c.settings.year) + '-' + str(c.settings.month).zfill(2) + '.nc'
    out = link + appendix
    return out


def check_other_instrument_idx(link, instrument_name):
    # dirty workaround to not discard instruments from further deployments... -.-
    out_links = []
    names = [instrument_name.lower(), instrument_name.lower().replace('-', '_')]
    for cur_link in link:
        cur_all_idx = []
        for cur_name in names:
            cur_all_idx.extend([n for n in xrange(len(cur_link)) if cur_link.find(cur_name, n) == n])
        cur_new_link = cur_link
        number_arr = range(0, 9)
        for number in number_arr:
            for idx in cur_all_idx:
                temp_len = len(cur_name)
                cur_new_link = cur_new_link[0:idx+temp_len-1] + str(number) + cur_new_link[idx+temp_len::]
            out_links.append(cur_new_link)
    return out_links


def check_other_deps(link):
    out_links = []
    dep_idx = link.find('dep000')
    number_arr = range(1, 6)
    for number in number_arr:
        cur_link = link
        cur_link = cur_link[0:dep_idx+6] + str(number) + cur_link[dep_idx+7:]
        out_links.append(cur_link)
    return out_links


def get_instrument_list(input_id):
    instrument_list_out = []
    instrument_types = []
    instrument_names = []
    link = 'http://apps.socib.es/DataDiscovery/mooring-services?id_platform=' + str(input_id)
    r = requests.Response()
    try:
        r = requests.get(link)
    except ValueError:
        print 'Cannot open ', link
    json = r.json()
    platform_name = json["name"]
    platform_type = json["platformType"]
    instrument_list = json["jsonInstrumentList"]
    for instrument in instrument_list:
        instrument_types.append(instrument["instrumentType"])
        instrument_list_out.append(create_netcdf_link(instrument["jwebChartLink"]))
        instrument_names.append(instrument["name"])
    return instrument_list_out, instrument_types, instrument_names, platform_name, platform_type


def get_data_array(data_array):
    """
    returns pure data in NetCDF variable (without mask)
    :param data_array: NetCDF Variable
    :return: data array (just [xxx])
    """
    if type(data_array.__array__()) is np.ma.masked_array:
        return data_array.__array__().data
    else:
        return data_array.__array__()


def convert_time(cur_time):
    """

    :param cur_time: posixtime from 1970...
    :return: string
    example: get_converted_time(1441065600.0)
        --> '2015-09-01 00:00:00'
    """
    return str(datetime.utcfromtimestamp(cur_time))


def get_qc_variable_name(variable):
    try:
        qc_variable_name = variable.ancillary_variables
    except AttributeError:
        logger.info("No QC variable found for " + variable.name)
        qc_variable_name = None
    return qc_variable_name


def replace_qc(data_array, logical_qc):
    data_array[logical_qc] = np.nan
    return data_array


def check_direction_variable(variable):
    try:
        cur_std_name = variable.standard_name
    except AttributeError:
        return False
    if cur_std_name in c.settings.quiver_direction_variables:
        return True
    else:
        return False


def get_standard_name(variable):
    try:
        std_name = variable.standard_name
        return std_name
    except AttributeError:
        try:
            std_name = variable.long_name
            return std_name
        except AttributeError:
            std_name = variable.name
            return std_name


def get_dimension(variable):
    return len(variable.dimensions)


def get_sampling_interval(time_array):
    if len(time_array) <= 1:
        return np.nan
    else:
        return time_array[1] - time_array[0]


def transform_y_label(y_label):
    if y_label == 'C':
        y_label = '$^{\circ}$C'
    elif y_label == 'degree':
        # y_label = '$^{\circ}$'
        y_label = 'degrees'
    elif y_label == 'm s-1':
        y_label = r'$m \, s^{-1}$'
    elif y_label == 'cm s-1':
        y_label = r'$cm \, s^{-1}$'
    elif y_label == '1':
        y_label = ''
    elif y_label == 'm/s':
        y_label = r'$m \, s^{-1}$'
    elif y_label == 'cm/s':
        y_label = r'$cm \, s^{-1}$'
    elif '-1' in y_label:
        y_label = y_label.replace('-1', '$^{-1}$')
    return y_label


def get_x_limits(timestamp):
    cur_year = datetime.utcfromtimestamp(timestamp).year
    cur_month = datetime.utcfromtimestamp(timestamp).month
    start_d = datetime(cur_year, cur_month, 1)
    if cur_month == 12:
        cur_year += 1
        cur_month = 1
    else:
        cur_month += 1
    end_d = datetime(cur_year, cur_month, 1)
    return [md.date2num(start_d), md.date2num(end_d)]


def plot_1d(doc, x_axis, data, y_label, plot_title, qc_data=None, s_name=None, x_limits=None, input_month_title=None,
            same_y_limits=False, is_degree_0_360_y_limit=False, is_degree_180_180_y_limit=False):
    if x_limits is None:
        x_limits = [x_axis[0], x_axis[-1]]
    plt.rcParams.update({'font.size': 13})
    y_label = transform_y_label(y_label)
    try:
        month_title = c.settings.month_str + ' ' + str(c.settings.year)
    except AttributeError:
        if input_month_title is None:
            logger.warning('Trying to access non-set setting variable.', exc_info=True)
            month_title = ''
        else:
            month_title = input_month_title
    if s_name is None:
        my_title = 'Evolution of ' + plot_title + ' in ' + month_title + '.'
    else:
        my_title = 'Evolution of ' + plot_title + ' at ' + s_name + ' in ' + month_title + '.'
    if np.all(np.isnan(data)):
        print 'Only NaN occured at 1D ' + plot_title
        doc.append(('NaN values at ' + plot_title + '.'))
        return
    with doc.create(Figure(position='htbp')) as plot:
        plt.figure()
        plt.plot(x_axis, data, color='k', linewidth=0.5)
        data_axis = plt.gca()
        data_axis.xaxis.set_major_formatter(c.settings.xfmt)
        data_axis.set_ylabel(y_label, rotation=0, horizontalalignment='right')
        data_axis.grid(b=False, which='major', color='k', linestyle='--', linewidth=0.25)
        if qc_data is None:
            pass
        else:
            if np.all(np.isnan(qc_data)):
                print 'QC data is set, but everything is NaN. Here we have a processing problem. ' + plot_title
                doc.append(('NaN values at QC of ' + plot_title + '. Should not happen.'))
                plt.gcf().autofmt_xdate()
                plt.gcf().set_size_inches(11, 3)
                plt.clf()
                plt.close('all')
                return
            qc_axis = plt.gca().twinx()
            qc_axis.plot(x_axis, qc_data, color='g', linewidth=0.75)
            qc_axis.set_ylabel('QC', rotation=0, horizontalalignment='left')
            qc_axis.set_ylim([0, 10])
        data_axis.set_xlim(x_limits)
        data_axis.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.title('\n'.join(wrap(plot_title, 70)))
        plt.gcf().autofmt_xdate()
        plt.gcf().set_size_inches(11, 3)
        data_axis.set_xticks(np.arange(x_limits[0], x_limits[1]+1, 5.0))
        has_set_ticks = False
        if is_degree_0_360_y_limit:
            # dirty, very dirty
            data_axis.set_ylim([0, 360])
            data_axis.set_yticks([0, 90, 180, 270, 360])
            has_set_ticks = True
        if is_degree_180_180_y_limit:
            data_axis.set_ylim([-180, 180])
            data_axis.set_yticks([-180, -90, 0, 90, 180])
            has_set_ticks = True
        if same_y_limits:
            cur_limits = data_axis.get_ylim()
            new_limit = np.max(np.abs(cur_limits))
            plt.ylim([-new_limit, new_limit])
        if not has_set_ticks:
            # data_axis.locator_params(tight=True, axis='y', nbins=6)
            data_axis.yaxis.set_major_locator(mp.ticker.MaxNLocator(nbins=7, symmetric=True, trim=False))
        plt.tight_layout()
        plt.subplots_adjust(top=0.8)
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption(my_title)
        plt.clf()
        plt.close('all')


def plot_2d(doc, x_axis, data, second_dimension, y_label, plot_title, qc_data=None, s_name=None, x_limits=None,
            month_title=None):
    if x_limits is None:
        x_limits = [x_axis[0], x_axis[-1]]
    plt.rcParams.update({'font.size': 13})
    y_label = transform_y_label(y_label)
    month_title = c.settings.month_str + ' ' + str(c.settings.year)
    if s_name is None:
        my_title = 'Evolution of ' + plot_title + ' in ' + month_title + '.'
    else:
        my_title = 'Evolution of ' + plot_title + ' at ' + s_name + ' in ' + month_title + '.'
    with doc.create(Figure(position='htbp')) as plot:
        plt.figure()
        for x in range(len(second_dimension)):
            if np.all(np.isnan(data)):
                print 'Only NaN occured at 2D ' + plot_title
                doc.append(('Only NaN values at ' + plot_title + '.'))
                continue
            plt.plot(x_axis, data[:, x], color='k', linewidth=0.5)
        plt.gca().set_ylabel(y_label, rotation=0, horizontalalignment='right')
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
        plt.gca().locator_params(axis='y', nbins=6)
        plt.gca().set_xticks(np.arange(x_limits[0], x_limits[1]+1, 5.0))
        plt.gca().grid(b=False, which='major', color='k', linestyle='--', linewidth=0.25)
        if qc_data is None:
            pass
        else:
            qc_axis = plt.gca().twinx()
            for x in range(len(second_dimension)):
                plt.plot(x_axis, qc_data[:, x], color='g', linewidth=0.75)
            qc_axis.set_ylabel('QC', rotation=0, horizontalalignment='left')
            qc_axis.set_ylim([0, 10])
        plt.gca().xaxis.set_major_formatter(c.settings.xfmt)
        plt.gca().set_xlim(x_limits)
        plt.title("\n".join(wrap(plot_title, 70)))
        plt.gcf().autofmt_xdate()
        plt.gcf().set_size_inches(11, 3)
        plt.tight_layout()
        plt.subplots_adjust(top=0.8)
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption(my_title)
        plt.clf()
        plt.close('all')


def plot_quiver_direction(doc, cur_time, direction, plot_title, amplifier=None, s_name=None, x_limits=None,
                          input_month_title=None):
    if x_limits is None:
        x_limits = [np.nanmin(cur_time) - 1, np.nanmax(cur_time) + 1]
    else:
        x_limits = [x_limits[0] - 1, x_limits[1] + 1]
    plt.rcParams.update({'font.size': 13})
    try:
        month_title = c.settings.month_str + ' ' + str(c.settings.year)
    except AttributeError:
        if input_month_title is None:
            logger.warning('Trying to access non-set setting variable.', exc_info=True)
            month_title = ''
        else:
            month_title = input_month_title
    if s_name is None:
        my_title = 'Evolution of ' + plot_title + ' in ' + month_title + '.'
    else:
        my_title = 'Evolution of ' + plot_title + ' at ' + s_name + ' in ' + month_title + '.'
    if np.all(np.isnan(direction)):
        print 'Only NaNs found at ' + plot_title
        doc.append('Only NaNs encountered at ' + plot_title)
        return
    if amplifier is not None:
        plot_title += ' normalized and scaled with speed'
        amplifier = normalize_data(amplifier, 'mean')
    u, v = compute_u_v_components(direction, amplifier)
    # u_raw = np.cos(np.deg2rad(direction))
    # v_raw = np.sin(np.deg2rad(direction))
    if amplifier is None:
        vector_length = np.sqrt((u * u) + (v * v))
        u = u / vector_length
        v = v / vector_length
    with doc.create(Figure(position='htbp')) as plot:
        plt.figure()
        try:
            plt.quiver(cur_time, 0, u, v, color='r', alpha=0.5, width=0.004, units='width', scale=1, scale_units='y',
                       headlength=0, headwidth=1)
        except RuntimeWarning:
            logger.debug('Problem at plotting quiver.', exc_info=True)
        plt.gca().set_xlim(x_limits)
        plt.gca().set_ylim(-1, 1)
        plt.gca().xaxis.set_major_formatter(c.settings.xfmt)
        plt.gca().set_xticks(np.arange(x_limits[0]+1, x_limits[1], 5.0))
        plt.title("\n".join(wrap(plot_title, 70)))
        plt.gcf().autofmt_xdate()
        plt.gcf().set_size_inches(11, 3)
        plt.tight_layout()
        plt.subplots_adjust(top=0.8)
        plot.add_plot(width=NoEscape(r'1\textwidth'))
        plot.add_caption(my_title)
        plt.clf()
        plt.close('all')


def normalize_data(amplifier, method=None):
    cur_min = np.nanmin(amplifier)
    cur_max = np.nanmax(amplifier)
    cur_mean = np.nanmean(amplifier)
    cur_std = np.std(amplifier)
    if method is None or method == 'mean':
        out = [((v-cur_min)/(cur_max - cur_min)) for v in np.nditer(amplifier)]
    elif method == 'std':
        out = [(v-cur_mean)/cur_std for v in np.nditer(amplifier)]
    else:
        print 'Normalization method not defined.'
        return 1
    out = np.asarray(out)
    return out


def get_sums_and_percents(cur_qc):
    good_idx = cur_qc == 1
    nan_idx = cur_qc == 9
    rest_idx = np.logical_or(good_idx, nan_idx)
    rest_idx = ~rest_idx
    sum_all = float(len(cur_qc))
    if np.any(rest_idx):
        prob_good_idx = cur_qc == 2
        sum_prob_good = np.count_nonzero(prob_good_idx)
        prob_bad_idx = cur_qc == 3
        sum_prob_bad = np.count_nonzero(prob_bad_idx)
        bad_idx = cur_qc == 4
        sum_bad = np.count_nonzero(bad_idx)
        spike_idx = cur_qc == 6
        sum_spike = np.count_nonzero(spike_idx)
    else:
        sum_prob_good = 0
        sum_prob_bad = 0
        sum_bad = 0
        sum_spike = 0
    sum_good = np.count_nonzero(good_idx)
    sum_nan = np.count_nonzero(nan_idx)
    percent_rest = "%.1f" % ((np.count_nonzero(rest_idx) / sum_all) * 100)
    percent_good = "%.1f" % ((sum_good / sum_all) * 100)
    percent_nan = "%.1f" % ((sum_nan / sum_all) * 100)
    return sum_good, sum_prob_good, sum_prob_bad, sum_bad, sum_spike, sum_nan, percent_rest, percent_good, percent_nan


def get_percents_only(level, cur_full_qc, counter):
    cur_qc = cur_full_qc
    if level != '':
        level = str(level) + ' m'
        cur_qc = cur_full_qc[:, counter]
        counter += 1
    (sum_good, sum_prob_good, sum_prob_bad, sum_bad, sum_spike, sum_nan, percent_rest, percent_good,
     percent_nan) = get_sums_and_percents(cur_qc)
    return (
        level, counter, sum_good, sum_prob_good, sum_prob_bad, sum_bad, sum_spike, sum_nan, percent_rest,
        percent_good, percent_nan)


class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print '[%s]' % self.name,
        print 'Elapsed: %s' % (time.time() - self.tstart)
