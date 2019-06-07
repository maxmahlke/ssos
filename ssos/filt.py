'''
    Python implementations of filter functions to detect SSOs in SExtractor and SCAMP catalogues

    Part of the ssos pipeline.
'''

import collections
import os

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np
import pandas as pd
from statsmodels import robust
import statsmodels.api as sm



def detections(sources, settings):
    '''
        Remove sources based on number of detections

        input
        ------
        sources - pandas dataframe, table of source data
        detections - list, containing integers responding to number of detections

        return
        ------
        sources - pandas dataframe, sources data without the sources filtered in this step
    '''
    measurements = collections.Counter(sources['SOURCE_NUMBER'])
    sources_to_remove = {key for key, value in measurements.items()
                         if value in settings['DETECTIONS']}

    sources = sources[~sources['SOURCE_NUMBER'].isin(sources_to_remove)]
    return sources


def proper_motion(sources, settings):
    '''
    Remove sources based on the proper motion range and signal-to-noise

    input
    ------
    sources - pandas dataframe, SCAMP full catalog
    settings - dictionary, settings parameter:value pairs

    return
    ------
    sources - pandas dataframe, sources data without the sources filtered in this step
    '''
    pm_low = settings['PM_LOW']
    pm_up = settings['PM_UP']
    SN = settings['PM_SNR']

    # Filter sources for proper motion range
    sources = sources[(sources.PM >= pm_low) & (sources.PM <= pm_up)]

    # Filter by proper motion SN
    sources = sources[sources.PM / sources.PMERR >= SN]
    return sources


def pixel(sources, settings):
    '''
    Remove sources based on pixel coordinates

    input
    ------
    sources - pandas dataframe, SCAMP full catalog
    settings - dictionary, contains the specified filter parameters

    return
    ------
    sources - pandas dataframe, sources data without the sources filtered in this step
    '''
    delta_pixel = settings['DELTA_PIXEL']

    # For all sources, check if the largest XWIN_IMAGE !and! YWIN_IMAGE differences
    # are smaller than delta_pixel.
    for source_number, group in sources.copy().groupby('SOURCE_NUMBER'):

        x_img = group['XWIN_IMAGE']
        y_img = group['YWIN_IMAGE']

        if x_img.max() - x_img.min() < delta_pixel:
            if y_img.max() - y_img.min() < delta_pixel:
                sources = sources[sources.SOURCE_NUMBER != source_number]

    return sources


def motion_is_linear(epoch, coords, sigma, r_squared):
    '''
    Motion fitting function

    input
    ------
    epoch - np.array, epoch of source
    ra - pd.Series, right ascension of source
    dec - pd.Series, declination of source
    r_squared - float, minimum r_squared for good fit

    return
    ------
    True if motion is linear else False
    '''
    epoch = sm.add_constant(epoch)
    model = sm.WLS(coords, epoch, weights=1/sigma**2)
    results = model.fit()

    if results.rsquared >= r_squared:
        return True
    else:
        return False

def motion_is_linear_in_ecliptic_coordinates(epoch, ra, dec, r_squ):
    # Transform to ecliptic coordinate system
    ecliptic_coords = SkyCoord(ra, dec, frame='icrs', unit='deg').barycentrictrueecliptic
    if not motion_is_linear(epoch, ecliptic_coords.lon.deg, 1, r_squ):
        return False
    else:
        if not motion_is_linear(epoch, ecliptic_coords.lat.deg, 1, r_squ):
            return False
        else:
            return True


def linear_motion(sources, settings):
    '''
        Remove sources based on linearity of motion and the Median Absolute Deviation

        input
        ------
        sources - pandas dataframe, table of source data
        settings - dict, 'parameter': value

        return
        ------
        sources - pandas dataframe, sources data without the sources filtered in this step
    '''

    r_squ, out_thres, id_out = settings['R_SQU_M'], settings['OUTLIER_THRESHOLD'],\
                               settings['IDENTIFY_OUTLIER']

     # Lower limit of subgroup size
    min_detections = 3

    sources_to_remove = []

    for source_number, group in sources.copy().groupby('SOURCE_NUMBER'):
        # GROUP contains the detections of a single source
        epoch = group['EPOCH'].values
        ra = group['ALPHA_J2000'].values
        dec = group['DELTA_J2000'].values

        if id_out:
            # Detect outliers using the Median Absolute Deviation of the EPOCHS
            outlier_index = np.where(abs(np.diff(epoch)) > robust.mad(epoch) *
                                                           out_thres)[0]
            if len(outlier_index) == 0:  # no outlier in EPOCH. Check edges
                outlier_in_middle = np.where(abs(np.diff(epoch[:-1])) > robust.mad(epoch[:-1]) *
                                                                         out_thres)[0]
                if len(outlier_in_middle) == 0:
                    # Linear in right ascencsion?
                    if not motion_is_linear(epoch, ra, group['ERRA_WORLD'].values, r_squ):
                        sources_to_remove.append(source_number)
                    else: # Linear motion in right ascension
                        # Linear motion in declination?
                        if not motion_is_linear(epoch, dec, group['ERRB_WORLD'].values, r_squ):
                            if not motion_is_linear_in_ecliptic_coordinates(epoch, ra, dec, r_squ):
                                sources_to_remove.append(source_number)
                    continue # source has linear motion and is kept

                else: # outlier found in middle (jump in epochs)
                    outlier_index = outlier_in_middle

            # Found outlier in EPOCH
            subgroups = np.split(group[['EPOCH', 'ALPHA_J2000', 'DELTA_J2000',
                                        'ERRA_WORLD', 'ERRB_WORLD']], outlier_index + 1)

            for subgroup in subgroups:

                if len(subgroup) < min_detections: # if there aren't sufficient detections
                    sources.loc[subgroup.index, 'FLAGS_SSOS'] +=1  # Add outlier flag

                else:
                    epoch = subgroup['EPOCH'].values
                    ra = subgroup['ALPHA_J2000'].values
                    dec = subgroup['DELTA_J2000'].values
                    # Linear in right ascencsion?
                    if not motion_is_linear(epoch, ra, subgroup['ERRA_WORLD'].values, r_squ):
                        sources.loc[subgroup.index, 'FLAGS_SSOS'] +=1  # Add outlier flag
                    else: # Linear motion in right ascension
                        # Linear motion in declination?
                        if not motion_is_linear(epoch, dec, subgroup['ERRB_WORLD'].values, r_squ):
                            if not motion_is_linear_in_ecliptic_coordinates(epoch, ra, dec, r_squ):
                                sources.loc[subgroup.index, 'FLAGS_SSOS'] +=1  # Add outlier flag

            # FLAG_SSOS is only uneven if detection is outlier
            # If all source detections were flagged, remove the source
            if all(sources[sources.SOURCE_NUMBER == source_number]['FLAGS_SSOS'] % 2 == 1):
                sources_to_remove.append(source_number)

        else:
            # Linear in right ascencsion?
            if not motion_is_linear(epoch, ra, group['ERRA_WORLD'].values, r_squ):
                sources_to_remove.append(source_number)

            else: # Linear motion in right ascension
                # Linear motion in declination?
                if not motion_is_linear(epoch, dec, group['ERRB_WORLD'].values, r_squ):
                    if not motion_is_linear_in_ecliptic_coordinates(epoch, ra, dec, r_squ):
                        sources_to_remove.append(source_number)

    sources = sources[~sources.SOURCE_NUMBER.isin(sources_to_remove)]
    return sources


def compare_std_and_sigma_of_average(y, sigma):
    ''' Compares weighted average uncertainty to standard deviation '''
    weights = 1 / sigma**2
    weighted_average_uncertainty = 1 / np.sqrt(np.sum(weights))
    standard_deviation_of_data = np.std(y)
    ratio = weighted_average_uncertainty / standard_deviation_of_data

    return ratio


def constant_trail(sources, settings):
    '''
    Remove sources based on constant trail size

    input
    ------
    sources - pd.DataFrame, source data
    settings - dict, 'parameter': value

    return
    ------
    sources - pandas dataframe, sources data without the sources filtered in this step
    '''

    ratio = settings['RATIO']

    for source_number, group in sources.copy().groupby('SOURCE_NUMBER'):

        awin = group['AWIN_IMAGE']
        bwin = group['BWIN_IMAGE']

        for i, dimension in enumerate([awin, bwin]):

            sigma = group['ERRAWIN_IMAGE'] if i == 0 else group['ERRBWIN_IMAGE']

            fitted_ratio = compare_std_and_sigma_of_average(dimension, sigma)

            if fitted_ratio >= ratio and i == 0:
                pass
            elif fitted_ratio >= ratio and i == 1:
                pass
            else:
                sources = sources[sources.SOURCE_NUMBER != source_number]

    return sources


def bright_sources_catalog(sources, settings):
    '''
        Remove sources based on proximity to stars. Requires the HYG stellar
        catalog.

        input
        ------
        sources - pandas dataframe, table of source data
        settings - dict, 'parameter': value

        return
        ------
        sources - pandas dataframe, sources data without the sources filtered in this step
    '''

    # If the bright sources cat is the SCAMP reference cat, we have to find the filename
    if type(settings['BRIGHT_SOURCES_CAT']) == list:
        refout_catpath, astref_catalog = settings['BRIGHT_SOURCES_CAT']
        try:
            settings['BRIGHT_SOURCES_CAT'] = [file for file in os.listdir(refout_catpath)
                                                       if file.startswith(astref_catalog)][0]
        except IndexError:
            from ssos.core import PipelineSettingsException
            raise PipelineSettingsException('Could not find SCAMP reference catalogue %s\n' \
                                             'Run SCAMP with SAVE_REFCATALOG set to Y to retrieve it to file.' % astref_catalog)
    try:
        with fits.open(settings['BRIGHT_SOURCES_CAT']) as ref_cat:
            bright_sources = Table(ref_cat[2].data).to_pandas()
    except OSError:
        try:
            bright_sources = pd.read_csv(settings['BRIGHT_SOURCES_CAT'])
        except ParserError:
            from ssos.core import PipelineSettingsException
            raise PipelineSettingsException('Could not read bright-sources catalogue %s ' % settings['BRIGHT_SOURCES_CAT'])

    bright_sources.rename(index=str, columns={'X_WORLD': 'RA', 'Y_WORLD': 'DEC'}, inplace=True)

    distance = settings['DISTANCE'] / 3600 # convert distance to degree
    lower_mag_limit, upper_mag_limt = [int(limit) for limit in settings['MAG_LIMITS'].split(',')]

    # Remove reference sources which are too faint or too bright
    bright_sources = bright_sources[(lower_mag_limit <= bright_sources.MAG) &
                                    (bright_sources.MAG <= upper_mag_limt)]

    # Remove bright_sources which are outside region of interest, given by the minimum and maximum
    # RA and Dec values of the source candiates plus the distance of the stellar regions
    upper_ra, lower_ra = (sources['ALPHA_J2000'].max() + distance,
                          sources['ALPHA_J2000'].min() - distance)
    upper_dec, lower_dec = (sources['DELTA_J2000'].max() + distance,
                            sources['DELTA_J2000'].min() - distance)

    bright_sources = bright_sources[(bright_sources['RA'] > lower_ra) &\
                                    (bright_sources['RA'] < upper_ra)]
    bright_sources = bright_sources[(bright_sources['DEC'] > lower_dec) &\
                                    (bright_sources['DEC'] < upper_dec)]

    for source_number, group in sources.copy().groupby('SOURCE_NUMBER'):

            mean_ra = np.mean(group['ALPHA_J2000'])
            mean_dec = np.mean(group['DELTA_J2000'])

            try:
                min_distance_to_star = min(np.sqrt((bright_sources['RA'] - mean_ra)**2 +
                                                   (bright_sources['DEC'] - mean_dec)**2).abs())

            except ValueError:  # reports "min() arg is an empty sequence"
                continue

            if min_distance_to_star < distance:
                sources = sources[sources.SOURCE_NUMBER != source_number]

    return sources