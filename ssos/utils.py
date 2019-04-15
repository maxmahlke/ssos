import argparse
import logging
import os
import sys
import time
import warnings

from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd


def init_argparse():
    '''
    Executes the setup of the arguments parser for the command-line api

    return
    ------
    args - command line argument parser
    '''

    parser = argparse.ArgumentParser(description='Pipeline to search for Solar System\
                                                  objects in wide-field imaging surveys')

    parser.add_argument('fields', nargs='+', help='Path to directory of field to be searched')

    parser.add_argument('-c', '--config', action='store', dest='set_up', nargs=1,
                        help='Path to configuration file')

    parser.add_argument('-d', '--default', action='store_true',
                        help='Copy default setup files to CWD')

    parser.add_argument('-i', '--inspect', action='store_true',
                        help='Launch visual inspection mode on SSO candidates sample')

    parser.add_argument('-l', '--log', action='store', dest='log', nargs=1,
                        help='Set the logging level. Valid arguments are DEBUG, INFO, WARNING, ERROR,\
                        CRITICAl.', default=['INFO'])

    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress logging to console')

    parser.add_argument('-t', '--target', action='store', dest='target', nargs=1,
                        help='Target directory to save fits files. If no target given, writing to CWD')

    parser.add_argument('--sex', action='store_true',
                        help='Force SExtractor runs')

    parser.add_argument('--scamp', action='store_true',
                        help='Force SCAMP runs')

    parser.add_argument('--swarp', action='store_true',
                        help='Force SWARP runs')

    parser.add_argument('--skybot', action='store_true',
                        help='Force SkyBoT query')

    group = parser.add_argument_group('Filter Settings')

    for prop in ['REMOVE_REF_SOURCES', 'FILTER_DETEC', 'FILTER_PM', 'FILTER_PIXEL', 'FILTER_MOTION', 'IDENTIFY_OUTLIER',
                 'FILTER_TRAIL', 'FILTER_BRIGHT_SOURCES', 'CROSSMATCH_SKYBOT',
                 'EXTRACT_CUTOUTS', 'FIXED_APER_MAGS']:

        group.add_argument('-' + prop, metavar='bool', action='store',
                           help='Override ' + prop + ' setting. Must be True or False.',
                           choices=['True', 'False'])


    for val in ['SCI_EXTENSION', 'WEIGHT_IMAGES', 'RA', 'DEC', 'OBJECT', 'DATE_OBS', 'FILTER', 'EXPTIME',
                'SEX_CONFIG', 'SEX_PARAMS', 'SEX_FILTER', 'SEX_NNW', 'SCAMP_CONFIG', 'ASTREF_CATALOG',
                'CROSSID_RADIUS', 'SWARP_CONFIG', 'DETECTIONS', 'PM_LOW', 'PM_UP', 'PM_SNR', 'DELTA_PIXEL',
                'OUTLIER_THRESHOLD', 'R_SQU_M', 'RATIO', 'BRIGHT_SOURCES_CAT', 'DISTANCE', 'MAG_LIMITS', 'CROSSMATCH_RADIUS',
                'CUTOUT_SIZE', 'REFERENCE_FILTER', 'OBSERVATORY_CODE', 'FOV_DIMENSIONS']:

        group.add_argument('-' + val, metavar='value', action='store',
                           help='Override ' + val + ' setting.')

    args = parser.parse_args()

    return args


def init_logger(args, log_dir):
    '''
    Initiate logging. By default, the log is written to file.
    The filename consists of the string 'sso_' followed by the current date and time.

    input
    ------
    args - command line argument parser
    log_dir - str, absolute path to log directory

    return
    ------
    log - logger instance
    log_file - str, absolute path to log file
    '''

    # See if logging level was set via argument. Else, default to INFO
    numeric_level = getattr(logging, args.log[0].upper(), None)

    if not isinstance(numeric_level, int):
        numeric_level = 20

    log = logging.getLogger(__name__)
    log.setLevel(numeric_level)

    start_time = time.struct_time(time.localtime())
    log_file = 'ssos_{:s}.log'.format(time.strftime('%Y%m%d%H%M%S', start_time))

    file_log = logging.FileHandler(os.path.join(log_dir, log_file))
    file_log.setLevel(numeric_level)
    file_log.terminator = '' # suppress newline characters
    log.addHandler(file_log)

    # If the quiet flag is set, logs to console are supressed
    if not args.quiet:
        stdout_log = logging.StreamHandler()
        stdout_log.setLevel(numeric_level)
        stdout_log.terminator = '' # suppress newline characters
        log.addHandler(stdout_log)

    def log_excepthook(execption_type, exception_value, traceback, logger=log):
        ''' Function to pipe tracebacks into the logfile '''
        log.error('Exception occurred during analysis: \n',
                  exc_info=(execption_type, exception_value, traceback))
        log.error('\n')

    sys.excepthook = log_excepthook
    return log, log_file, start_time


def create_target_dir(args):
    '''
    Creates the target directory file structure

    input
    ------
    args - command line argument parser

    return
    ------
    target_dir - str, absolute path to target directory
    paths - dict, dictionary containing the absolute paths of the subdirectories
    '''

    if not args.target:
        target_dir = os.getcwd()
    else:
        target_dir = os.path.abspath(args.target[0])

    paths = {
        'cats':    os.path.join(target_dir, 'cats'),
        'cutouts': os.path.join(target_dir, 'cutouts'),
        'logs':     os.path.join(target_dir, 'logs'),
        'images':  os.path.abspath(args.fields[0]),
        'skybot':  os.path.join(target_dir, 'skybot'),
        'tmp':  os.path.join(target_dir, 'tmp'),
        'weights':  os.path.join(target_dir, 'weights')
    }

    for _, path in paths.items():
        os.makedirs(path, exist_ok=True)

    return target_dir, paths


def query_skybot(epoch, ra, dec, fov, obs_code):
    '''
    Retrieves the SkyBoT query for the specified parameters

    input
    ------
    epoch - float, mid-exposure epoch in MJD
    ra - float, right ascension in degree
    dec - float, declination in degree
    fov - str, value of FOV_DIMENSIONS parameter
    obs_code - str, value of OBSERVATORY_CODE parameter

    return
    -------
    skybot - pd.DataFrame, skybot sources in the FoV. Empty dataframe if no sources
    '''

    obs = {
        '-ep': str(epoch), '-ra': str(ra),
        '-dec': str(dec),  '-bd': fov,
        '-mime': 'text',  '-loc': obs_code,
        '-output': 'obs'
    }

    url = 'http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php?'
    url = '&'.join([url, *['='.join([k, v]) for k, v in obs.items()]])

    skybot = pd.read_csv(url, skiprows=2, delimiter='|')

    if skybot.empty:
        print(skybot.columns)
        return skybot
    # Fix column names
    skybot = skybot.rename({col: col.replace(' ', '') for col in skybot.columns},
                           axis=1)

    # Convert RA and DEC columns to degree
    skybot['RA'] = skybot.apply(lambda x: SkyCoord(x['RA(h)'], x['DE(deg)'], frame='icrs', unit=(u.hourangle, u.deg)).ra.deg, axis=1)
    skybot['DEC'] = skybot.apply(lambda x: SkyCoord(x['RA(h)'], x['DE(deg)'], frame='icrs', unit=(u.hourangle, u.deg)).dec.deg, axis=1)
    skybot.drop(['RA(h)', 'DE(deg)'], inplace=True, axis=1)

    skybot['EPOCH'] = epoch

    return skybot


def unpack_header_kw(hdus, keywords, try_first=0):
    '''
    Looks up keyword values in FITS headers. Iterates over
    extensions until keyword is found, else returns False

    hdus - HDUList, exposure headers
    keywords - string or list, keywords to look up
    try_first - header index to look up first

    returns
    ----
    vals - single parameter value or list of parameter values
    '''
    vals = []

    kws = [keywords] if not isinstance(keywords, list) else keywords

    for kw in kws:
        for hdu_ind in set([try_first, *range(len(hdus))]):
            try:
                hdu = hdus[hdu_ind]
                vals.append(hdu.header[kw])
                break
            except (IndexError, KeyError):
                pass
        else:
            vals.append(False)

    if not isinstance(keywords, list):
        return vals[0]

    return vals