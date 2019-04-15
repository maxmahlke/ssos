from collections import Counter
import os
import time
import sys

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd

import ssos.filt as filt
from ssos.opt import extract_cutouts
from ssos.opt import compute_aperture_magnitudes
from ssos.opt import crossmatch_skybot
from ssos.utils import create_target_dir
from ssos.utils import init_argparse
from ssos.utils import init_logger
from ssos.utils import query_skybot
from ssos.utils import unpack_header_kw



FILTER_STEPS = {
        'FILTER_DETEC':        filt.detections,
        'FILTER_PM':           filt.proper_motion,
        'FILTER_MOTION':       filt.linear_motion,
        'FILTER_PIXEL':        filt.pixel,
        'FILTER_TRAIL':        filt.constant_trail,
        'FILTER_BRIGHT_SOURCES': filt.bright_sources_catalog
        }


ANALYSIS_STEPS = {
    'CROSSMATCH_SKYBOT': crossmatch_skybot,
    'EXTRACT_CUTOUTS':   extract_cutouts,
    'FIXED_APER_MAGS':   compute_aperture_magnitudes
    }


class Pipeline:
    ''' Class to handle the pipeline execution '''

    def __init__(self):

        # Preparations: Handling of command line arguments
        self.args = init_argparse()


        # Assert that images are found and contain the required header keywords
        self.images = [os.path.join(os.path.abspath(self.args.fields[0]), image) for image in
                       os.listdir(os.path.abspath(self.args.fields[0])) if image.endswith('.fits')]
        assert len(self.images) > 0, 'No images found in %s! Ensure that they have a .fits'\
                                     ' extension' % os.path.abspath(self.args.fields[0])

        self.target_dir, self.paths = create_target_dir(self.args)
        self.log, self.log_file, self.start_time = init_logger(self.args, self.paths['logs'])
        self.log.info('\n\t--- The ssos Pipeline ---\t\t--- {:s} ---\n\n'
                      .format(time.strftime('%Y/%m/%d %H:%M:%S', self.start_time)))


        # Reading and checking the settings
        self.settings = self._set_settings()
        self.settings = self._check_settings(self.settings)

        self._print_field_info()

        self.steps = [param for param in FILTER_STEPS.keys() if self.settings[param]]
        self.analysis_steps = [param for param in ANALYSIS_STEPS.keys() if self.settings[param]]
        self.added_proper_motion = False
        self.added_SExtractor_data = False

        # Little hack for fixed apertuer magnitudes "without" cutout extraction
        if not self.settings['EXTRACT_CUTOUTS'] and self.settings['FIXED_APER_MAGS']:
            self.log.debug('FIXED_APER_MAGS is True but EXTRACT_CUTOUTS is False. Will'\
                           ' delete cutouts upon finishing.')
            self.paths['cutouts'] = self.paths['tmp']
            self.analysis_steps[0:0] = ['EXTRACT_CUTOUTS']


    def _set_settings(self):
        '''
        Set pipeline settings

        input
        ------
        args - command line arguments parser
        log - logging instance

        return
        ------
        settings - dict, dictionary containing PARAMETER: VALUE pairs
        '''

        settings = {}

        try:
            if self.args.set_up:
                set_up = open(self.args.set_up[0], 'r')
            else:
                set_up = open('default.ssos', 'r')

        except IOError:
            self.log.debug('No config file provided in CWD or via -c flag.\
                            Using default settings.\n')
            try:
                set_up = open(os.path.join(os.path.dirname(__file__),
                                           'default.ssos'))
            except IOError:
                raise PipelineSettingsException('No configuration file provided '\
                       'and none found in package directory.')

        with set_up:
            for line in set_up:
                if line == '\n' or line[0] == '#':
                    continue

                parameter, value = line.split()[:2]
                settings[parameter] = value

        # Override values from command line
        for arg in vars(self.args):
            if arg in settings.keys() and vars(self.args)[arg] != None:
                settings[arg] = vars(self.args)[arg]

        # Expand $HOME to home directory
        settings = {k: v.replace('$HOME', os.path.expanduser('~')) for k, v in settings.items()}

        # Filter setting finished. Write them to the logfile
        with open(os.path.join(self.paths['logs'], self.log_file), 'a') as logfile:

            for setting, value in settings.items():
                line = setting.ljust(20) + str(value) + '\n'
                self.log.debug(line)
                logfile.write(line)

            logfile.write('\n')
            self.log.debug('\n')

        return settings


    def _check_settings(self, settings):
        '''
        Converts parameter values to expected formats

        input
        -----
        settings - dict, dictionary containing PARAMETER: VALUE pairs

        return
        -----
        settings - dict, dictionary containing PARAMETER: VALUE pairs
        '''

        # Unpack provided DETECTIONS and SCI extension(s) into list
        for param in ['SCI_EXTENSION', 'DETECTIONS']:
            try:
                settings[param] = [int(character) for character in settings[param].split(',')]
            except ValueError:
                if param == 'SCI_EXTENSION' and settings[param] == 'All':
                    settings[param] = False # treat it as 'None provided'
                else:
                    raise PipelineSettingsException('%s value invalid' % param)
        if self.settings['SCI_EXTENSION']:
            self.sci_ext = self.settings['SCI_EXTENSION'][0]
        else:
            self.sci_ext = 0

        # Check one image header for keyword presence
        with fits.open(self.images[0]) as exposure:
            kws = ['RA', 'DEC', 'OBJECT', 'DATE-OBS', 'FILTER', 'EXPTIME']
            for kw in kws:
                if not unpack_header_kw(exposure, self.settings[kw], self.sci_ext):
                    raise PipelineSettingsException('Could not find keyword %s in FITS header.'\
                                                    ' Is the SCI_EXTENSION correct?' % self.settings[kw])
                if kw == 'DATE-OBS':
                    # Find format of DATE-OBS keyword
                    try:
                        Time(unpack_header_kw(exposure, kw, try_first=self.sci_ext), format='isot')
                        self.date_obs_fmt = 'isot'
                    except ValueError:
                        self.date_obs_fmt = 'mjd'


        # Check that config files exist
        for file in ['SEX_CONFIG', 'SEX_PARAMS', 'SEX_NNW', 'SEX_FILTER',
                     'SCAMP_CONFIG', 'SWARP_CONFIG']:
            if not os.path.isfile(settings[file]):
                raise PipelineSettingsException('Could not find %s in %s' %
                                                (file, settings[file]) )

        # Check if weight images are provided and exist
        if not settings['WEIGHT_IMAGES'].upper() == 'FALSE':
            if not os.path.isdir(settings['WEIGHT_IMAGES']):
                raise PipelineSettingsException('Could not find weight images directory %s.' %
                                                 settings['WEIGHT_IMAGES'])
        else:
            settings['WEIGHT_IMAGES'] = False

        # Convert filter strings to booleans
        for param in ['REMOVE_REF_SOURCES', 'FILTER_DETEC', 'FILTER_PM', 'FILTER_PIXEL', 'FILTER_MOTION',
                      'IDENTIFY_OUTLIER', 'FILTER_TRAIL',
                      'FILTER_BRIGHT_SOURCES', 'CROSSMATCH_SKYBOT', 'EXTRACT_CUTOUTS',
                      'FIXED_APER_MAGS']:

            if settings[param].upper() == 'TRUE':
                settings[param] = True
            elif settings[param].upper() == 'FALSE':
                settings[param] = False
            else:
                raise PipelineSettingsException('Could not evaluate %s value, has to be True or\
                                                 False' % param)

        if settings['FILTER_MOTION']:
            if not settings['FILTER_DETEC'] or\
              (not 1 in settings['DETECTIONS'] and not 2 in settings['DETECTIONS']):
                raise PipelineSettingsException('When FILTER_MOTION is True, DETECTIONS needs\
                                                 to contain "1,2".')

        # Convert numeric values to float
        for param in ['PM_LOW', 'PM_UP', 'PM_SNR', 'DELTA_PIXEL', 'OUTLIER_THRESHOLD',
                      'R_SQU_M', 'RATIO', 'DISTANCE', 'CROSSMATCH_RADIUS',
                      'CUTOUT_SIZE']:

            try:
                settings[param] = float(settings[param])
            except ValueError:
                raise PipelineSettingsException('Could not convert %s value to float'
                                                % param)
            if param == ['R_SQU_M']:
                assert 0 <= settings[param] <= 1, 'The %s parameter has to be in range [0 - 1]' % param

        if settings['FILTER_BRIGHT_SOURCES']:
            # Evaluate bright-sources catalogue path
            if settings['BRIGHT_SOURCES_CAT'] == 'REFCAT':
                # get the filename and column names from the scamp config file
                with open(settings['SCAMP_CONFIG'], 'r') as file:
                    for line in file:
                        if 'REFOUT_CATPATH' in line:
                            refout_catpath = line.split()[1]
                        if 'ASTREF_CATALOG' in line:
                            astref_catalog = line.split()[1]
                if self.args.ASTREF_CATALOG: # if overwritten via command line
                    astref_catalog = self.args.ASTREF_CATALOG

                settings['BRIGHT_SOURCES_CAT'] = [refout_catpath, astref_catalog]

        return settings


    def _print_field_info(self):
        ''' Prints RA, DEC, and OBJECT keywords to log '''
        with fits.open(self.images[0]) as exposure:
            ra, dec, object = unpack_header_kw(exposure, ['RA', 'DEC', 'OBJECT'], try_first=self.sci_ext)

        ecli_lat = SkyCoord(ra, dec, frame='icrs', unit='deg').barycentrictrueecliptic.lat.deg

        self.log.info('\t|\t'.join(['%i Exposures' % len(self.images),
                                    '%s' % object,
                                    '%.2fdeg Ecliptic Latitude\n\n' % ecli_lat]))

        # Get number of known SSOs in FoV from SkyBoT
        if self.settings['CROSSMATCH_SKYBOT']:
            self.log.info('Querying SkyBoT for known SSOs in FoV..\n')

            if not self.args.skybot and \
               os.path.isfile(os.path.join(self.paths['skybot'], 'skybot_all.csv')):

               self.skybot = pd.read_csv(os.path.join(self.paths['skybot'], 'skybot_all.csv'))

            else:
                self.skybot = pd.DataFrame()

                for img in self.images:
                    with fits.open(img) as exp:
                        ra, dec, date_obs, texp = unpack_header_kw(exp, [self.settings['RA'], self.settings['DEC'],
                                                                         self.settings['DATE-OBS'], self.settings['EXPTIME']], self.sci_ext)

                        mid_epoch = (Time(date_obs, format=self.date_obs_fmt) + float(texp) / 2 * u.second).isot

                        if self.settings['FOV_DIMENSIONS'] == '0x0':

                            # Determine exposure FoV
                            naxis1, naxis2, cdelt1, cdelt2 = unpack_header_kw(exp, ['NAXIS1', 'NAXIS2', 'CDELT1', 'CDELT2'], self.sci_ext)

                            if not cdelt1 or not cdelt2:
                                cd11, cd12, cd21, cd22 = unpack_header_kw(exp, ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'], self.sci_ext)
                                dra  = naxis1 * abs(cd11) + naxis2 * abs(cd12)
                                ddec = naxis1 * abs(cd21) + naxis2 * abs(cd22)

                            else:
                                dra = naxis1 * cdelt1
                                ddec = naxis2 * cdelt2

                            # Ensure coverage by adding CROSSMATCH_RADIUS
                            fov = '{:.1f}x{:.1f}'.format(round(dra  + self.settings['CROSSMATCH_RADIUS'] / 60, 1),
                                                         round(ddec + self.settings['CROSSMATCH_RADIUS'] / 60, 1))

                        else:
                            fov = self.settings['FOV_DIMENSIONS']

                        result = query_skybot(mid_epoch, ra, dec, fov,
                                              self.settings['OBSERVATORY_CODE'])
                        if not result.empty:
                            self.skybot = self.skybot.append(result)

                if not self.skybot.empty:
                    self.skybot.to_csv(os.path.join(self.paths['skybot'], 'skybot_all.csv'), index=False)

            if not self.skybot.empty:
                # Print gimmicky SkyBoT info
                bins = np.arange(np.floor(min(self.skybot.Mv)),
                                 np.ceil(max(self.skybot.Mv)), 0.5)

                counts = pd.cut(self.skybot['Mv'], bins).value_counts()
                counts.sort_index(inplace=True)

                magnitudes = '  '.join([str(i) for i in range(int(min(bins)),
                                                             int(np.ceil(max(bins))) + 1
                                                             )])
                # construct bar chart top-down
                bars = []
                bar_height = 7

                for i in range(1, bar_height):
                    bar = ' '
                    for count in counts:
                        if count >= max(counts.values) * (bar_height-i)/bar_height:
                            bar += '##'
                        else:
                            bar += '  '
                    bars.append(bar)
                # bottom bar
                bar = ' '
                for count in counts:
                    if count > 0:
                        bar += '##'
                    else:
                        bar += '  '
                bars.append(bar)

                self.log.info('\nSkyBoT: %i SSOs with %i detections\n'
                              '\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n %s\n%sMv\n'
                               % (len(set(self.skybot.Name)),
                                  len(self.skybot),
                                  *bars, magnitudes,
                                  ' ' * int(len(magnitudes) / 2)))
            else:
                self.log.info('No known SSOs returned by SkyBoT.\n')

    def _run_SExtractor_on_single_image(self, image, extension):
        '''
        Run SExtractor on individual images.

        input
        ------
        image - str, absolute path to image file
        extension - int, index of FITS file extension to be SExtracted

        return
        ------
        cat, str - absolute path to SExtractor output catalog
        date_obs, str - observation epoch from exposure
        '''

        # Name of output catalog
        cat = os.path.join(self.paths['cats'],
                           os.path.splitext(os.path.basename(image))[0])

        # Select extension or run whole image
        if extension is not False:
            image_ext = image + '[%i]' % extension
            cat += '_%i.cat' % extension
            ext = extension

        else:
            image_ext = image
            cat += '.cat'
            ext = 1

        if not self.args.sex and os.path.isfile(cat):
            self.log.debug('SExtractor catalog %s already exists! Skipping this sextraction..\n' % cat)

            with fits.open(image) as exposure:
                date_obs = unpack_header_kw(exposure, self.settings['DATE-OBS'], self.sci_ext)

            return cat, date_obs

        sex_args = {

            'file': image_ext,
            'config': self.settings['SEX_CONFIG'],
            'overwrite_params': { # Arguments to initialize Astromatic class
                    'CATALOG_NAME': cat,
                    'PARAMETERS_NAME': self.settings['SEX_PARAMS'],
                    'FILTER_NAME': self.settings['SEX_FILTER'],
                    'STARNNW_NAME': self.settings['SEX_NNW'],
            },
        }

        if self.settings['WEIGHT_IMAGES']:
            sex_args['overwrite_params']['WEIGHT_IMAGE'] = 'MAP_WEIGHT'

            if extension is not False:
                weight_suffix = '_%i.weight' % extension
            else:
                weight_suffix = '.weight'

            sex_args['overwrite_params']['WEIGHT_IMAGE'] = os.path.join(self.settings['WEIGHT_IMAGES'],
                                                           os.path.basename(image).replace(
                                                           '.fits', weight_suffix))

        if self.log.level <= 10:  # if we're at DEBUG log level, print SExtractor output
            sex_args['overwrite_params']['VERBOSE_TYPE'] = 'NORMAL'

        # ------
        # Exectue SExtractor
        cmd = ' '.join(['sex', sex_args['file'], '-c', sex_args['config']])
        for param, value in sex_args['overwrite_params'].items():
            cmd += ' '.join([' -' + param, value])

        self.log.debug('\nExecuting SExtractor command:\n%s\n\n' % cmd)
        os.system(cmd)

        # ------
        with fits.open(image) as exposure:
            # Make .ahead file with the EPOCH in MJD for SCAMP
            # Following http://www.astromatic.net/forum/showthread.php?tid=501
            date_obs = unpack_header_kw(exposure, self.settings['DATE-OBS'], self.sci_ext)
            mjd = Time(date_obs, format=self.date_obs_fmt).mjd

            with open(os.path.splitext(cat)[0] + '.ahead', 'w+') as file:
                for hdu in exposure:
                    file.write('MJD-OBS = %.6f\nEND\n' % mjd)
        return cat, date_obs


    def run_SExtractor(self):
        ''' Wrapper for SExtractor run.  Adds progress bar.'''

        self.SExtractor_catalogues = []

        observation_epochs = []
        # Add progress bar
        if not self.args.quiet:
            sys.stdout.write('[%s]' % (' ' * len(self.images)))
            sys.stdout.flush()
            sys.stdout.write('\b' * (len(self.images) + 1))

        for image in self.images:

            if self.settings['SCI_EXTENSION']:

                for extension in self.settings['SCI_EXTENSION']:
                    cat, date_obs = self._run_SExtractor_on_single_image(image, extension)
                    observation_epochs.append(date_obs)
                    self.SExtractor_catalogues.append(cat)

            else:
                cat, date_obs = self._run_SExtractor_on_single_image(image, False)
                observation_epochs.append(date_obs)
                self.SExtractor_catalogues.append(cat)

            if not self.args.quiet:
                sys.stdout.write('-')
                sys.stdout.flush()

        # Sort the SExtractor catalogues by their observation epochs, important for SCAMP
        self.SExtractor_catalogues = [cat for _, cat in sorted(zip(observation_epochs,
                                                               self.SExtractor_catalogues))]


    def run_SCAMP(self, crossid_radius=None, full_name='full.cat',
                  merged_name='merged.cat', keep_refcat=False,
                  solve_astronomy=True, pattern_matching=True,
                  adjust_SExtractor_and_aheader=False):
        ''' Run SCAMP on SExtractor catalogs '''

        # ------
        # Set-up SCAMP run
        self.merged_cat = os.path.join(self.paths['cats'], merged_name)
        self.full_cat = os.path.join(self.paths['cats'], full_name)

        scamp_args = {
            'files':  self.SExtractor_catalogues,
            'config': self.settings['SCAMP_CONFIG'],

            'overwrite_params': {
                'MERGEDOUTCAT_NAME': self.merged_cat,
                'FULLOUTCAT_NAME':   self.full_cat,
                },
            }

        # SCAMP appends an index for the Field Group to the output catalogues
        self.merged_cat = self.merged_cat.replace('.cat', '_1.cat')
        self.full_cat = self.full_cat.replace('.cat', '_1.cat')

        if self.args.ASTREF_CATALOG is not None:
            scamp_args['overwrite_params']['ASTREF_CATALOG'] = self.args.ASTREF_CATALOG

        if self.settings['FILTER_BRIGHT_SOURCES'] and \
           self.settings['BRIGHT_SOURCES_CAT'] == 'REFCAT':
            scamp_args['overwrite_params']['SAVE_REFCATALOG'] = 'Y'

        if not solve_astronomy:
            scamp_args['overwrite_params']['SOLVE_ASTROM'] = 'N'
            scamp_args['overwrite_params']['INCLUDE_ASTREFCATALOG'] = 'Y'

        if not pattern_matching:
            scamp_args['overwrite_params']['MATCH'] = 'N'

        if crossid_radius is not None:
            scamp_args['overwrite_params']['CROSSID_RADIUS'] = str(crossid_radius)

        if not self.args.scamp and os.path.isfile(self.merged_cat) and os.path.isfile(self.full_cat):
            if not pattern_matching: # suppress this message in case SCAMP is run twice
                self.log.info('\nReading SCAMP catalogues from file..\n')
            scamp_from_file = True

        else:
            self.log.info('\nRunning SCAMP..\t')

            # ------
            # Exectue SCAMP
            catalogs = ' '.join(scamp_args['files'])
            cmd = ' '.join(['scamp', catalogs, '-c', scamp_args['config']])
            for param, value in scamp_args['overwrite_params'].items():
                cmd += ' '.join([' -' + param, value])

            self.log.debug('\nExecuting SCAMP command:\n%s\n' % cmd)
            os.system(cmd)
            scamp_from_file = False

        # Create the full catalog as pipeline property
        with fits.open(self.full_cat) as full:
            data = Table(full[2].data)

        self.sources = data.to_pandas()

        if not keep_refcat:
            # Catalogue number 0 are the reference sources
            self.sources = self.sources[self.sources['CATALOG_NUMBER'] != 0]

        if (self.args.scamp or not scamp_from_file) \
            and adjust_SExtractor_and_aheader:
            self.adjust_SExtractor_catalogues()
            self.adjust_ahead_files()

        # Ensure that EPOCHs are in order
        self.sources = self.sources.sort_values(['SOURCE_NUMBER', 'EPOCH'],
                                                ascending=[True, True])
        self.sources = self.sources.reset_index()

        # Add flag to dataframe
        self.sources['FLAGS_SSOS'] = 0


    def adjust_SExtractor_catalogues(self):
        # Flags source detections which were associated to reference sources
        # These will subsequently be ignored by SCAMP
        self.log.info('Adjusting SExtractor catalogues for pattern matching..\n')
        REMOVE_DETECTIONS = 2 # 1 image + reference, lower limit

        detections = Counter(self.sources.SOURCE_NUMBER) # this includes references
        flag_count = [source for source, count in detections.items() if count >= REMOVE_DETECTIONS]

        # only remove sources which were associated to a reference
        flag = []
        for source_number in flag_count:
            if 0 in self.sources[self.sources.SOURCE_NUMBER == source_number]['CATALOG_NUMBER'].values:
                flag.append(source_number)

        self.sources = self.sources[self.sources.SOURCE_NUMBER.isin(flag)]
        self.sources = self.sources[self.sources.CATALOG_NUMBER != 0]

        for cat_number, catalogue in self.sources.groupby('CATALOG_NUMBER'):
            with fits.open(self.SExtractor_catalogues[cat_number-1]) as cat:

                for extension, detections in catalogue.groupby('EXTENSION'):
                    extension *= 2 # extensions are HEADER and DATA

                    xwin_img = cat[extension].data.field('XWIN_IMAGE')
                    detections_to_flag = np.in1d(xwin_img, detections.X_IMAGE)
                    cat[extension].data['FLAGS'][detections_to_flag] = 128

                cat.writeto(self.SExtractor_catalogues[cat_number-1], overwrite=True)

    def adjust_ahead_files(self):
        self.log.info('Copying astrometric solution to aheader files..\n')
        aheaders = [cat.replace('.cat', '.ahead') for cat in self.SExtractor_catalogues]

        for ahead in aheaders:
            with open(ahead, 'r') as ahead_file:
                mjd_obs = ahead_file.readline()

            with open(ahead, 'w') as ahead_file:
                with open(ahead.replace('.ahead', '.head'), 'r') as head_file:
                    for line in head_file:
                        if 'END' in line:
                            line = mjd_obs + line
                        ahead_file.write(line)


    def number_of_sources(self):
        return len(set(self.sources.SOURCE_NUMBER))


    def check_remaining(func):
        ''' Wrapper to check if source candidates remain after filter step'''
        def wrapped_func(*args, **kwargs):

            result = func(*args, **kwargs)

            if not args[0].number_of_sources() > 0:
                raise NoSourcesRemainingException('\nNo source candidates left after %s.'
                                                   ' Stopping pipeline.\n' % args[1])
                sys.exit()

            return result
        return wrapped_func


    @check_remaining
    def execute_filter(self, step_name):
        '''
        Calls the proper filter function

        input
        ------
        step_name - str, name of filter step

        return
        ------
        sources - pd.DataFrame, remaining source candidates
        '''

        # The proper motion filter has to access the merged catalogue first
        if step_name == 'FILTER_PM':
            self.add_proper_motion()

        # The trail filters have to access the SExtractor catalogues first
        if step_name in ['FILTER_PIXEL', 'FILTER_TRAIL', 'FILTER_T_DIST']:
            self.add_SExtractor_data()

        self.sources = FILTER_STEPS[step_name](self.sources, self.settings)


    def add_proper_motion(self):
        '''
        Adds proper motion values from merged to full catalog.  Has to be called
        either before filtering the proper motion or at the end of the pipeline.
        '''
        with fits.open(self.merged_cat) as merged:
            data = Table(merged[2].data)

        proper_motion_columns = ['PMALPHA_J2000', 'PMALPHAERR_J2000',
                                 'PMDELTA_J2000', 'PMDELTAERR_J2000']

        merged = data[['SOURCE_NUMBER'] + proper_motion_columns].to_pandas()

        # Add proper motions column from merged catalog to source database
        self.sources = self.sources.merge(merged, on='SOURCE_NUMBER')

        for column in proper_motion_columns:
            self.sources[column] = self.sources[column].apply(lambda x: x / 8.75e6)

        # Add proper motion column.  Factor 8.75e6 converts to "/h
        self.sources['PM'] = self.sources.apply(lambda x:
                               np.sqrt(x['PMALPHA_J2000']**2 + x['PMDELTA_J2000']**2),
                               axis=1)
        self.sources['PMERR'] = self.sources.apply(lambda x:
                               np.sqrt(x['PMALPHAERR_J2000']**2 + x['PMDELTAERR_J2000']**2),
                               axis=1)

        self.added_proper_motion = True


    def add_SExtractor_data(self):
        '''
        Adds the trail morphology parameters from the SExtractor catalogues to
        the source database.  Has to be called either before filtering on trail
        morphology or at the end of the pipeline.
        '''

        # Add SExtractor data
        morphometry_columns = ['AWIN_IMAGE', 'BWIN_IMAGE', 'THETAWIN_IMAGE',
                               'XWIN_IMAGE', 'YWIN_IMAGE', 'ERRAWIN_IMAGE',
                               'ERRBWIN_IMAGE', 'ERRTHETAWIN_IMAGE', 'FLUX_AUTO',
                               'FLUXERR_AUTO']
        self.sources['XWIN_IMAGE'] = self.sources['X_IMAGE']

        sextractor_data = pd.DataFrame(index=[], columns=morphometry_columns)


        for catalogue in self.SExtractor_catalogues:
            with fits.open(catalogue) as cat:
                for j, header in enumerate(cat):
                    if j % 2 == 0 and j != 0:
                        data = Table(header.data)[:][morphometry_columns]
                        sextractor_data = sextractor_data.append(data.to_pandas(),
                                                     ignore_index=True, sort=True)

        self.sources = pd.merge(self.sources, sextractor_data[morphometry_columns],
                                on='XWIN_IMAGE', how='left')
        self.added_SExtractor_data = True


    def add_image_metadata(self):

        # Derive image filename from SExtractor catalogue
        for cat_number, group in self.sources.groupby('CATALOG_NUMBER'):

            if self.settings['SCI_EXTENSION']:
                image_filename = '_'.join(os.path.splitext(
                                            os.path.basename(self.SExtractor_catalogues[cat_number-1])
                                                          )[0].split('_')[:-1]) + '.fits'
            else:
                image_filename = os.path.basename(self.SExtractor_catalogues[cat_number-1]).replace('.cat', '.fits')

            self.sources.loc[group.index, 'IMAGE_FILENAME'] = image_filename

        # Add image metadata. Have to use the correct header extension
        for image_filename, group in self.sources.groupby(['IMAGE_FILENAME']):

            extension = group.EXTENSION.values[0]
            # Add exposure keywords
            with fits.open(os.path.join(self.paths['images'], image_filename)) as exposure:
                for prop in ['OBJECT', 'DATE-OBS', 'FILTER', 'EXPTIME', 'RA_IMAGE', 'DEC_IMAGE']:
                    self.sources.loc[group.index, prop] = date_obs = unpack_header_kw(exposure, self.settings[prop.split('_')[0]], extension)

        self.sources['MID_EXPOSURE_MJD'] = self.sources.apply(lambda x: ( Time(x['DATE-OBS'], format=self.date_obs_fmt) +
                                                                       float(x['EXPTIME']) / 2 * u.second).mjd, axis=1)


    def execute_analysis(self, step_name):
        '''
        Calls the proper analysis function

        input
        ------
        step_name - str, name of analysis step
        '''
        self.sources = ANALYSIS_STEPS[step_name](self.sources, self.settings, self.log,
                                      self.paths, self.args)


    def save_and_cleanup(self):
        ''' Clean up the final database.  Rename columns, remove unnecessary columns. '''

        self.sources.drop(columns=['ASTR_INSTRUM', 'PHOT_INSTRUM'], inplace=True)

        self.sources.rename(index=str, columns={'ALPHA_J2000': 'RA', 'DELTA_J2000': 'DEC', 'PMALPHA_J2000': 'PMRA',
                                                'PMALPHAERR_J2000': 'PMRA_ERR', 'PMDELTA_J2000': 'PMDEC',
                                                'PMDELTAERR_J2000': 'PMDEC_ERR', 'SKYBOT_PMALPHA': 'SKYBOT_PMRA',
                                                'SKYBOT_PMDELTA': 'SKYBOT_PMDEC', 'SKYBOT_ALPHA': 'SKYBOT_RA',
                                                'SKYBOT_DELTA': 'SKYBOT_DEC', 'FLUX_AUTO': 'FLUX',
                                                'FLUXERR_AUTO': 'FLUXERR'}, inplace=True)

        # Rearrange columns
        ordered_columns = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'RA', 'DEC', 'EPOCH', 'MAG', 'MAGERR', 'FLUX', 'FLUXERR',
                           'PM', 'PMERR', 'PMRA',  'PMRA_ERR', 'PMDEC', 'PMDEC_ERR', 'MID_EXPOSURE_MJD', 'DATE-OBS',
                           'EXPTIME', 'OBJECT', 'FILTER', 'RA_IMAGE', 'DEC_IMAGE', 'IMAGE_FILENAME',
                           'EXTENSION', 'XWIN_IMAGE', 'YWIN_IMAGE', 'AWIN_IMAGE', 'ERRAWIN_IMAGE', 'BWIN_IMAGE',
                           'ERRBWIN_IMAGE', 'THETAWIN_IMAGE', 'ERRTHETAWIN_IMAGE', 'ERRA_WORLD', 'ERRB_WORLD',
                           'ERRTHETA_WORLD', 'FLAGS_EXTRACTION', 'FLAGS_SCAMP', 'FLAGS_IMA', 'FLAGS_SSOS']

        if self.settings['FIXED_APER_MAGS']:
            ordered_columns[7:7] = ['MAG_CI', 'MAGERR_CI', 'FLUX_CI', 'FLUXERR_CI']

        if self.settings['CROSSMATCH_SKYBOT']:
            ordered_columns[22:22] = ['SKYBOT_NUMBER', 'SKYBOT_NAME', 'SKYBOT_CLASS', 'SKYBOT_MAG',
                                      'SKYBOT_RA', 'SKYBOT_DEC', 'SKYBOT_PMRA', 'SKYBOT_PMDEC']

        self.sources = self.sources[ordered_columns]

        # Save final database
        output = Table.from_pandas(self.sources)
        output_filename = os.path.join(self.paths['cats'], 'ssos_{:s}.csv'.format(time.strftime('%Y%m%d%H%M%S', self.start_time)))
        output.write(output_filename, format='csv', overwrite=True)

        os.system('rm -r %s' % self.paths['tmp'])
        self.run_time = time.time() - time.mktime(self.start_time)


class PipelineSettingsException(Exception):
    pass


class NoSourcesRemainingException(Exception):
    pass
