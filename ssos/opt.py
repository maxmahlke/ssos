#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
    Author: Max Mahlke
    Date: 16 December 2019

    Extraction of cutouts of SSOs using SWarp,
    cross-matching with the SkyBoT database, and other analyses.
'''

import os

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

import ssos.utils as utils


def query_skybot(images, settings, sci_ext, date_obs_fmt):
    '''Queries the VO SkyBoT service for known objects
    in the images

    :images: list - paths to images
    :settings: dict - pipeline settings
    :sci_ext: int or False - number of science extension
    :date_obs_fmt: str - format of observation date string
    :returns: pd.DataFrame - SkyBoT response for all images
    '''
    skybot = pd.DataFrame()

    for img in tqdm(images, unit='imgs',
                    desc='Querying SkyBoT for known SSOs in FoV'):
        # ------
        # Construct the query string
        with fits.open(img) as exp:
            ra, dec, = utils.compute_image_center(exp[sci_ext].header)
            date_obs, texp = \
                utils.unpack_header_kw(exp, [settings['DATE-OBS'],
                                             settings['EXPTIME']], sci_ext)
            mid_epoch = (Time(date_obs, format=date_obs_fmt) +
                         (float(texp) / 2) * u.second).isot

            if settings['FOV_DIMENSIONS'] == '0x0':  # we determine the FoV

                naxis1, naxis2, cdelt1, cdelt2 = \
                    utils.unpack_header_kw(exp, ['NAXIS1', 'NAXIS2',
                                                 'CDELT1', 'CDELT2'],
                                           sci_ext)

                if cdelt1 is False or cdelt2 is False:
                    cd11, cd12, cd21, cd22 = \
                        utils.unpack_header_kw(exp, ['CD1_1', 'CD1_2',
                                                     'CD2_1', 'CD2_2'],
                                               sci_ext)
                    dra = naxis1 * abs(cd11) + naxis2 * abs(cd12)
                    ddec = naxis1 * abs(cd21) + naxis2 * abs(cd22)

                else:
                    dra = naxis1 * cdelt1
                    ddec = naxis2 * cdelt2

                # Ensure full coverage by adding CROSSMATCH_RADIUS
                dra = round(dra  + settings['CROSSMATCH_RADIUS'] / 60, 1)
                ddec = round(ddec + settings['CROSSMATCH_RADIUS'] / 60, 1)
                fov = f'{dra:.1f}x{ddec:.1f}'

            else:
                fov = settings['FOV_DIMENSIONS']

        result = _query_skybot_per_image(mid_epoch, ra, dec, fov,
                                         settings['OBSERVATORY_CODE'])
        if not result.empty:
            skybot = skybot.append(result)

    return skybot


def _query_skybot_per_image(epoch, ra, dec, fov, obs_code):
    ''' Retrieves the SkyBoT query for the specified parameters

    input
    ------
    epoch - float, mid-exposure epoch in MJD
    ra - float, right ascension in degree
    dec - float, declination in degree
    fov - str, value of FOV_DIMENSIONS parameter
    obs_code - str, value of OBSERVATORY_CODE parameter

    return
    -------
    skybot - pd.DataFrame, skybot sources in the FoV
    '''

    obs = {
        '-ep': str(epoch), '-ra': str(ra),
        '-dec': str(dec),  '-bd': fov,
        '-mime': 'text',  '-loc': obs_code,
        '-from': 'ssos', '-refsystem': 'EQJ2000',
        '-output': 'obs'
    }

    url = 'http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php?'
    url = '&'.join([url, *['='.join([k, v]) for k, v in obs.items()]])

    # log.debug(url)

    skybot = pd.read_csv(url, skiprows=2, delimiter='|')

    if skybot.empty:
        return skybot

    # Fix column names
    skybot = skybot.rename({c: c.replace(' ', '') for c in skybot.columns},
                           axis=1)

    # Convert RA and DEC columns to degree. SCAMP works with ICRS,
    # SkyBoT returns FK5 coordinates
    skybot['RA'] = skybot.apply(lambda x: SkyCoord(x['RA(h)'], x['DE(deg)'],
                                                   frame='fk5',
                                                   unit=(u.hourangle, u.deg)
                                                   ).icrs.ra.deg,
                                axis=1)
    skybot['DEC'] = skybot.apply(lambda x: SkyCoord(x['RA(h)'], x['DE(deg)'],
                                                    frame='fk5',
                                                    unit=(u.hourangle, u.deg)
                                                    ).icrs.dec.deg,
                                 axis=1)

    skybot.drop(['RA(h)', 'DE(deg)'], inplace=True, axis=1)
    skybot['EPOCH'] = epoch

    return skybot


def _call_swarp(row, settings, log, paths, args):

    image_dir    = args.fields[0]
    cutout_dir   = paths['cutouts']
    cutout_size  = settings['CUTOUT_SIZE']
    swarp_config = settings['SWARP_CONFIG']

    img_ext = settings['SCI_EXTENSION'][row['EXTENSION'] - 1] if settings['SCI_EXTENSION'] else row['EXTENSION'] - 1
    image_file = os.path.join(image_dir, row['IMAGE_FILENAME']) + '[%i]' % img_ext

    cutout_filename = os.path.join(cutout_dir, '{:.0f}_{:02d}.fits'.format(row['SOURCE_NUMBER'],
                                                                           row['CATALOG_NUMBER']))

    if not args.swarp and os.path.isfile(cutout_filename):
        log.debug('Cutout %s exists, skipping..' % cutout_filename)
        return False

    swarp_args = {
                'file': image_file,
                'config': swarp_config,
                'overwrite_params': {
                    'WEIGHTOUT_NAME': os.path.join(paths['tmp'], 'tmp_weight'),
                    'CENTER'        : '{:.5f},{:.5f}'.format(row['ALPHA_J2000'], row['DELTA_J2000']),
                    'IMAGEOUT_NAME' : '{:s}'.format(cutout_filename),
                    'IMAGE_SIZE'    : str(cutout_size),
                    'COPY_KEYWORDS' : 'OBJECT,RA,DEC,FILTER,DATE-OBS'},
                }
    # ------
    # Build command and execute SWARP
    cmd = ' '.join(['swarp', swarp_args['file'], '-c', swarp_args['config']])
    for param, value in swarp_args['overwrite_params'].items():
        cmd += ' '.join([' -' + param, value])


    log.debug('\nExecuting SWARP command:\n%s' % cmd)
    if log.level > 10:  # Shown SWARP warnings only when debugging
        cmd += ' >/dev/null 2>&1'
    os.system(cmd)

    # Clean bad image keywords
    utils.create_clean_image(cutout_filename, cutout_filename, update=True)


def extract_cutouts(sources, settings, log, paths, args):
    '''
        Extract cutouts of sources found in SCAMP catalogues using SWARP

        input
        ------
        sources - pd.DataFrame, table of source data
        settings - dict, dictionary of pipeline settings with PARAM:VALUE pairs
        log - logging object
        paths - dict, absolute paths of target directories
        args - command line argument parser

        return
        ------
        sources - pd.DataFrame, table of source data (unchanged)
    '''

    log.info('\nExtracting cutouts with SWARP..')

    tqdm.pandas(desc='Creating cutouts', unit='cutouts')
    sources.progress_apply(lambda row: _call_swarp(row, settings, log, paths, args), axis=1)
    return sources


def compute_aperture_magnitudes(sources, settings, log, paths, args):
    '''
    Executes SExtractor in dual image mode and computes fixed aperture magnitudes

    input
    ------
    sources - pd.DataFrame, table of source data
    settings - dict, dictionary of pipeline settings with PARAM:VALUE pairs
    log - logging object
    paths - dict, absolute paths of target directories
    args - command line argument parser

    return
    ------
    sources - pd.DataFrame, table of source data with new magnitudes
    '''

    log.info('Computing fixed aperture magnitudes with SExtractor..')

    cutout_dir = paths['cutouts']
    cutouts = [os.path.join(cutout_dir, cutout) for cutout in os.listdir(cutout_dir)
               if cutout.endswith('.fits')]
    reference_filter = settings['REFERENCE_FILTER'].split(',')

    for number, source in sources.groupby('SOURCE_NUMBER'):

        # Select detection image based on reference filter
        filter_found = False
        for filter_ in reference_filter:
            try:
                reference_detection = source.loc[source['FILTER'] == filter_].iloc[0]
                filter_found = True
            except IndexError:
                continue
            if filter_found:
                sources.loc[(sources.SOURCE_NUMBER == reference_detection.SOURCE_NUMBER) &
                            (sources.CATALOG_NUMBER == reference_detection.CATALOG_NUMBER), 'FLAGS_SSOS'] += 2
                break
        else:
            from ssos.core import PipelineSettingsException
            raise PipelineSettingsException('No detection of candidate %i found in reference filter %s.'\
                                            % (number, reference_filter))

        detection_image = os.path.join(cutout_dir, '{:.0f}_{:02d}.fits'.format(number,
                                                                               reference_detection['CATALOG_NUMBER']))

        # ------
        # Call SExtractor in dual image mode
        for image in cutouts:

            if str(number) not in image:
                continue

            cat = os.path.join(paths['cutouts'], os.path.basename(image.replace('.fits', '.cat')))

            if image == detection_image:
                detection_cat = cat

            sex_args = {
                'file': ','.join([detection_image, image]),
                'config': settings['SEX_CONFIG'],
                'overwrite_params': {
                    'PARAMETERS_NAME': settings['SEX_PARAMS'],
                    'CATALOG_NAME'   : cat,
                    'WEIGHT_TYPE'    : 'BACKGROUND',
                    'FILTER_NAME'    : settings['SEX_FILTER'],
                    'STARNNW_NAME'   : settings['SEX_NNW'],
                },
            }

            # ------
            # Exectue SExtractor
            cmd = ' '.join(['sex', sex_args['file'], '-c', sex_args['config']])
            for param, value in sex_args['overwrite_params'].items():
                cmd += ' '.join([' -' + param, value])
            os.system(cmd)

        # Identify the SSO in the detection cat to track it among the other ones
        cutout_size = settings['CUTOUT_SIZE']

        with fits.open(detection_cat) as detection_cat:

            # ------
            # Find the source closes to the center coordinate
            x_images = detection_cat[2].data.field('XWIN_IMAGE')
            y_images = detection_cat[2].data.field('YWIN_IMAGE')

            # Get 2d distances of sources to image center
            distances = [np.linalg.norm(np.array([cutout_size, cutout_size]) / 2 - np.array([x_img, y_img]))
                         for x_img, y_img in zip(x_images, y_images)]

            # Find index of source closest to center -> SSO candidate
            center_source_ind = np.where(distances == np.min(distances))[0]

        for ind, detection in source.iterrows():

            catalog = os.path.join(cutout_dir, '{:.0f}_{:02d}.cat'.format(detection.SOURCE_NUMBER,
                                                                           detection.CATALOG_NUMBER))

            with fits.open(catalog) as cat:
                for prop in ['MAG', 'FLUX']:
                    sources.loc[ind, prop + '_APER'] = cat[2].data.field(prop + '_AUTO')[center_source_ind]
                    sources.loc[ind, prop + '_APER_ERR'] = cat[2].data.field(prop + 'ERR_AUTO')[center_source_ind]

    log.info('\tDone.\n')

    return sources


def _cross_match(source, skybot, radius):
    '''
        Cross-matching is applied to single rows in database.
        Checks for matches and adds SkyBoT source parameters
        to database if a match is found.

        input
        ------
        source - pd.Series, single source detection parameters
        skybot - pd.DataFrame, parameters of all SkyBoT sources in exposure
        radius - float, crossmatch radius in arcsec

        returns
        source - pd.Series, single source detection parameters
                 with SkyBot parameters added
    '''

    # Find SkyBoT source that is closest in angular separation to source
    dists = abs(np.sqrt((np.cos(source['DELTA_J2000']) *\
                       (skybot['RA'] - source['ALPHA_J2000']))**2 +\
                       (skybot['DEC'] - source['DELTA_J2000'])**2))

    if dists.min() <= radius:
        match = skybot.iloc[dists.values.argmin()]

        source['MATCHED']         = True
        source['SKYBOT_NUMBER']   = int(match['#Num']) \
                                    if match['#Num'].strip().isdigit()\
                                    else np.nan
        source['SKYBOT_NAME']     = match['Name']
        source['SKYBOT_MAG']      = match['Mv']
        source['SKYBOT_RA']       = match['RA']
        source['SKYBOT_DEC']      = match['DEC']
        source['SKYBOT_PMRA']     = match['dRA(arcsec/h)']
        source['SKYBOT_PMDEC']    = match['dDEC(arcsec/h)']
        source['SKYBOT_CLASS']    = match['Class']
        source['SKYBOT_DELTARA']  = 3600 * (source.ALPHA_J2000 - match['RA']) *\
                                           np.cos(source.DELTA_J2000)
        source['SKYBOT_DELTADEC'] = 3600 * (source.DELTA_J2000 - match['DEC'])
    else:
        source['MATCHED'] = False

    return source


def _compute_pm_difference_angle(row):
    # returns the angle between the skybot and source proper motion vectors
    skybot_pm = np.array([row.SKYBOT_PMALPHA, row.SKYBOT_PMDELTA])

    try:
        skybot_pm /= np.linalg.norm(skybot_pm)

    except ValueError: # no skybot match
        return 100.    # should never be the largest angle

    source_pm = np.array([row.PMALPHA_J2000, row.PMDELTA_J2000])
    source_pm /= np.linalg.norm(source_pm)

    return np.arccos(np.clip(np.dot(source_pm, skybot_pm), -1.0, 1.0))


def crossmatch_skybot(sources, settings, log, paths, args):
    '''
    Query SkyBoT for each epoch and cross-match the expected SSOs
    with the sources in the database.  Requires working internet connection.

    input
    ------
    sources - pd.DataFrame, table of source data
    settings - dict, dictionary of pipeline settings with PARAM:VALUE pairs
    log - logging object
    paths - dict, absolute paths of target directories
    args - command line argument parser

    return
    ------
    sources - pd.DataFrame, table of source data with SkyBoT matches added
    '''
    log.info('\nCross-matching SSO candidates with SkyBoT..')

    crossmatch_radius = settings['CROSSMATCH_RADIUS'] / 3600

    skybot_columns = ['MATCHED', 'SKYBOT_NAME', 'SKYBOT_NUMBER',
                      'SKYBOT_CLASS', 'SKYBOT_MAG', 'SKYBOT_RA',
                      'SKYBOT_DEC', 'SKYBOT_PMRA', 'SKYBOT_PMDEC',
                      'SKYBOT_DELTARA', 'SKYBOT_DELTADEC']

    # The default value for the SkyBoT columns are empty strings
    sources = pd.concat([sources, pd.DataFrame('', sources.index,
                         skybot_columns)], axis=1)

    try:
        skybot = pd.read_csv(os.path.join(paths['cats'], 'skybot_all.csv'))
    except pd.errors.EmptyDataError:
        # If no SkyBoT objects where found in FoV
        return sources

    # ------
    # Cross-match for each mid-exposure epoch
    for mid_exp_mjd, obs in sources.groupby('MID_EXPOSURE_MJD'):

        skybot_obs = skybot[skybot.EPOCH == Time(mid_exp_mjd,
                                                 format='mjd').isot]
        sources.loc[obs.index,
                    skybot_columns] = obs.apply(_cross_match,
                                                skybot=skybot_obs,
                                                radius=crossmatch_radius,
                                                axis=1)

    # For each SSO, check if multiple matches were found
    # If so, choose the match which is most co-linear
    #  with the measured proper motion (minimize the angles)
    for number, source in sources.groupby('SOURCE_NUMBER'):

        if len(set(source.SKYBOT_NAME)) > 1:
            # if one of the matches is NaN,
            #  remove it and see if there's only one left
            names = set([name for name in source.SKYBOT_NAME if name != ''])

            if len(names) == 1: # some detections were not matched,
                                #  add the designation to their columns

                numbers = [numb for numb in source.SKYBOT_NUMBER if numb != '']
                sources.loc[source.index, 'SKYBOT_NAME'] = list(names)[0]
                sources.loc[source.index, 'SKYBOT_NUMBER'] = list(numbers)[0]

            else: # compute the best match by
                  # minimizing the angular difference in proper motion

                pm_difference = source.apply(_compute_pm_difference_angle,
                                             axis=1)
                best_match = pm_difference.idxmin()

                best_match_name = sources.loc[best_match, 'SKYBOT_NAME']
                best_match_number = sources.loc[best_match, 'SKYBOT_NUMBER']

                # Overwrite metadata of detections not matched to best match
                sources.loc[source[source.SKYBOT_NAME !=
                                   best_match_name].index, skybot_columns] =\
                                       [best_match_name, best_match_number,
                                        '', '', '', '', '', '']

    n_matches = len(set(sources[sources.MATCHED]['SOURCE_NUMBER']))
    log.info(f' {n_matches} matched.\n')
    return sources


def create_checkplots(sources, settings, log, paths, args):
    ''' Creates checkplots and saves them to checkplot dir

    input
    ------
    sources - pd.DataFrame, table of source data
    settings - dict, dictionary of pipeline settings with PARAM:VALUE pairs
    log - logging object
    paths - dict, absolute paths of target directories
    args - command line argument parser
    '''

    for cp in settings['CHECKPLOTS']:

        skybot = sources[sources.MATCHED]

        if cp == 'SKYBOT_RESIDUALS':

            if skybot.empty: # when SkyBot reply is empty
                log.info(f'No SkyBoT matches. Skipping checkplot {cp}.\n')
                return sources

            fig, ax = plt.subplots()

            ax.scatter(skybot.SKYBOT_DELTARA, skybot.SKYBOT_DELTADEC,
                       color='black', s=2)

            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()

            ax.plot([xmin, xmax], [0, 0], ls='--', lw=0.5, color='gray')
            ax.plot([0, 0], [ymin, ymax], ls='--', lw=0.5, color='gray')

            ax.set(xlabel='(RA - SKYBOT_RA) * cos(DEC) / arcsec',
                   ylabel='DEC - SKYBOT_DEC / arcsec', xlim=(xmin, xmax),
                   ylim=(ymin, ymax))

            plt.tight_layout()
            fig.savefig(os.path.join(paths['checkplots'], f'{cp.lower()}.png'),
                        dpi=400)

        elif cp == 'SKYBOT_PM':

            if skybot.empty: # when SkyBot reply is empty
                log.info(f'No SkyBoT matches. Skipping checkplot {cp}.\n')
                return sources

            fig, ax = plt.subplots()

            skybot_PM = np.sqrt(
                            np.power(skybot.SKYBOT_PMRA.astype('float'), 2)\
                          + np.power(skybot.SKYBOT_PMDEC.astype('float'), 2)
                               )

            ax.scatter(skybot.PM, skybot_PM,
                       color='black', s=2)

            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()

            ax.plot([xmin, xmax], [xmin, xmax], ls='--', lw=0.5, color='gray')

            ax.set(xlabel='PM / (arcsec / h)',
                   ylabel=f'sqrt((SKYBOT_PMRA)**2 + (SKYBOT_PMDEC)**2) / '
                          f'(arcsec / h)', xlim=(xmin, xmax),
                   ylim=(ymin, ymax))

            plt.tight_layout()
            fig.savefig(os.path.join(paths['checkplots'], f'{cp.lower()}.png'),
                        dpi=400)

        else:
            log.info(f'Did not understand checkplot {cp}. Skipping.\n')

    return sources
