'''
    Extraction of cutouts of SSOs using SWarp,
    cross-matching with the SkyBoT database, and other analyses.

    Part of the SSO recovery pipeline.
'''

import os

from astropy.io import fits
from astropy.time import Time
import numpy as np
import pandas as pd


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

    image_dir    = args.fields[0]
    cutout_dir   = paths['cutouts']
    cutout_size  = settings['CUTOUT_SIZE']
    swarp_config = settings['SWARP_CONFIG']

    for index, row in sources.iterrows():
        img_ext = settings['SCI_EXTENSION'][row['EXTENSION'] - 1] if settings['SCI_EXTENSION'] else row['EXTENSION'] - 1
        image_file = os.path.join(image_dir, row['IMAGE_FILENAME']) + '[%i]' % img_ext

        cutout_filename = os.path.join(cutout_dir, '{:.0f}_{:02d}.fits'.format(row['SOURCE_NUMBER'],
                                                                               row['CATALOG_NUMBER']))

        if not args.swarp and os.path.isfile(cutout_filename):
            log.debug('Cutout %s exists, skipping..' % cutout_filename)
            continue

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

    log.info('\tDone.\n')
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

        detection_image = os.path.join(cutout_dir, '_'.join([str(number),
                                                             str(reference_detection['CATALOG_NUMBER'])])
                                                             + '.fits')

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

            catalog = os.path.join(cutout_dir, '_'.join([str(detection.SOURCE_NUMBER),
                                                         str(detection.CATALOG_NUMBER)]) + '.cat')

            with fits.open(catalog) as cat:
                for prop in ['MAG', 'FLUX']:
                    sources.loc[ind, prop + '_CI'] = cat[2].data.field(prop + '_AUTO')[center_source_ind]
                    sources.loc[ind, prop + 'ERR_CI'] = cat[2].data.field(prop + 'ERR_AUTO')[center_source_ind]

    log.info('\tDone.\n')

    return sources


def _cross_match(source, skybot, radius):
    '''
        Cross-matching is applied to single rows in database.  Checks for matches
        and adds SkyBoT source parameters to database if a match is found.

        input
        ------
        source - pd.Series, single source detection parameters
        skybot - pd.DataFrame, parameters of all SkyBoT sources in exposure
        radius - float, crossmatch radius in arcsec

        returns
        source - pd.Series, single source detection parameters with SkyBot parameters added
    '''

    # Find SkyBoT source that is closest in angular separation to source
    distances_to_source = abs(np.sqrt((skybot['RA'] - source['ALPHA_J2000'])**2 +\
                                      (skybot['DEC'] - source['DELTA_J2000'])**2))
    if distances_to_source.min() <= radius:
        skybot_match = skybot.iloc[distances_to_source.values.argmin()]

        source['SKYBOT']         = True
        source['SKYBOT_NUMBER']  = skybot_match['#Num']
        source['SKYBOT_NAME']    = skybot_match['Name']
        source['SKYBOT_MAG']     = skybot_match['Mv']
        source['SKYBOT_ALPHA']   = skybot_match['RA']
        source['SKYBOT_DELTA']   = skybot_match['DEC']
        source['SKYBOT_PMALPHA'] = skybot_match['dRA(arcsec/h)']
        source['SKYBOT_PMDELTA'] = skybot_match['dDEC(arcsec/h)']
        source['SKYBOT_CLASS']   = skybot_match['Class']
    else:
        source['SKYBOT'] = False

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

    skybot = pd.read_csv(os.path.join(paths['skybot'], 'skybot_all.csv'))
    crossmatch_radius = settings['CROSSMATCH_RADIUS'] / 3600 # convert from arcsecond to degrees

    skybot_columns = ['SKYBOT', 'SKYBOT_NAME', 'SKYBOT_NUMBER', 'SKYBOT_CLASS', 'SKYBOT_MAG',
                      'SKYBOT_ALPHA', 'SKYBOT_DELTA', 'SKYBOT_PMALPHA', 'SKYBOT_PMDELTA']

    # The default value for the SkyBoT columns are empty strings
    sources = pd.concat([sources, pd.DataFrame('', sources.index,
                        skybot_columns)], axis=1)


    # Cross-match for each mid-exposure epoch
    for mid_exp_mjd, detections in sources.groupby('MID_EXPOSURE_MJD'):

        skybot_detections = skybot[skybot.EPOCH == Time(mid_exp_mjd, format='mjd').isot]
        sources.loc[detections.index, skybot_columns] = detections.apply(_cross_match, skybot=skybot_detections,
                                                                                       radius=crossmatch_radius, axis=1)


    # For each SSO, check if multiple matches were found
    # If so, choose the match which is most co-linear with the measured proper motion (minimize the angles)
    for number, source in sources.groupby('SOURCE_NUMBER'):

        if len(set(source.SKYBOT_NAME)) > 1:
            # if one of the matches is NaN, remove it and see if there's only one left
            names = set([name for name in source.SKYBOT_NAME if name != ''])

            if len(names) == 1: # some detections were not matched, add the designation to their columns

                numbers = [numb for numb in source.SKYBOT_NUMBER if numb != '']
                sources.loc[source.index, 'SKYBOT_NAME'] = list(names)[0]
                sources.loc[source.index, 'SKYBOT_NUMBER'] = list(numbers)[0]

            else: # compute the best match by minimizing the angular difference in proper motion

                pm_difference = source.apply(_compute_pm_difference_angle, axis=1)
                best_match = pm_difference.idxmin()

                best_match_name = sources.loc[best_match, 'SKYBOT_NAME']
                best_match_number = sources.loc[best_match, 'SKYBOT_NUMBER']

                sources.loc[source[source.SKYBOT_NAME != best_match_name].index,
                                  skybot_columns] = [best_match_name, best_match_number,
                                                                         '', '', '', '', '', '']

    log.info(' %i matched.\n' % len(set(sources[sources.SKYBOT]['SOURCE_NUMBER'])))
    return sources