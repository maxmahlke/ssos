import argparse
import logging
import os
import sys
import time


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

    parser.add_argument('-t', '--target', action='store', dest='target', nargs=1,
                        help='Target directory to save fits files. If no target given, writing to CWD')

    parser.add_argument('-l', '--log', action='store', dest='log', nargs=1,
                        help='Set the logging level. Valid arguments are DEBUG, INFO, WARNING, ERROR,\
                        CRITICAl.', default=['INFO'])

    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress logging to console')

    parser.add_argument('--sex', action='store_true',
                        help='Force SExtractor runs')

    parser.add_argument('--scamp', action='store_true',
                        help='Force SCAMP runs')

    parser.add_argument('--swarp', action='store_true',
                        help='Force SWARP runs')
    parser.add_argument('--skybot', action='store_true',
                        help='Force SkyBoT query')


    group = parser.add_argument_group('Filter Settings')

    for prop in ['FILTER_DETEC', 'FILTER_PM', 'FILTER_PIXEL', 'FILTER_MOTION', 'IDENTIFY_OUTLIERS',
                 'FILTER_TRAIL', 'FILTER_T_DIST', 'FILTER_STAR_REGIONS', 'CROSSMATCH_SKYBOT',
                 'EXTRACT_CUTOUTS', 'FIXED_APER_MAGS']:

        group.add_argument('-'+prop, metavar='bool', action='store',
                           help='Override ' + prop + ' setting. Must be True or False.',
                           choices=['True', 'False'])


    for val in ['SEX_CONFIG', 'SEX_PARAMS', 'SEX_FILTER', 'SEX_NNW', 'SCAMP_CONFIG', 'SWARP_CONFIG',
                'SCI_EXTENSION', 'WEIGHT_IMAGES', 'DETECTIONS', 'PM_LOW', 'PM_UP', 'PM_SNR', 'DELTA_PIXEL',
                'OUTLIER_THRESHOLD', 'R_SQU_M', 'R_SQU_T', 'SIGMA', 'DISTANCE', 'HYGCAT',
                'CROSSMATCH_RADIUS', 'CUTOUT_SIZE', 'REFERENCE_FILTER', 'OBSERVATORY_CODE',
                'FOV_DIMENSIONS']:

        group.add_argument('-'+val, metavar='value', action='store',
                           help='Override ' + val + ' setting.')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit()

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

    log_file = 'sso_{}.log'.format(time.strftime('%Y%m%d%H%M%S'))

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
    return log, log_file


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

