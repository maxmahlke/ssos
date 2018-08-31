###############
Getting started
###############

.. role:: bash(code)
   :language: bash


.. role:: python(code)
   :language: python

Installing the Pipeline
=======================

Clone the `GitHub Repository <https://github.com/maxmahlke/SSO_Pipeline>`_ or download the `zip file <https://github.com/maxmahlke/SSO_Pipeline/archive/master.zip>`_.

The pipeline requires the following additional `python` packages: `astropy`, `numpy`, `pandas`, `scipy`, and `statsmodels`. You can quickly install them using the following command within the package directory

.. code-block:: bash

    $ pip install -r requirements.txt

Make sure that the astrOmatic binaries :bash:`sex`, :bash:`scamp`, and :bash:`swarp` are in your :bash:`PATH` shell variable. Try it with e.g.

.. code-block:: bash

    $ sex --version
    SExtractor version 2.19.5 (2015-06-19)
    $ scamp --version
    SCAMP version 2.0.4 (2015-06-19)
    $ swarp --version
    SWarp version 2.38.1 (2018-06-28)




Pipeline Setting Files
======================

The default ``pipeline_settings.ssos`` file is can be `found here <https://github.com/maxmahlke/SSO_Pipeline/blob/master/sso_pipeline/pipeline_settings.ssos>`_. It is a plain ASCII, designed very similar to the configuration files of SExtractor and SCAMP in order to make the user feel right at home. Short descriptions and expected values of the parameters are below, for more detailed descriptions refer to the `Pipeline <pipelins.rst>`_ and `Implementation <implementation.rst>`_ pages.

.. _Guide to SExtractor: http://astroa.physics.metu.edu.tr/MANUALS/sextractor/Guide2source_extractor.pdf


+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| Parameter             | Values  | Examples                |Description                                                                |
+=======================+=========+=========================+===========================================================================+
| `SEX_CONFIG`          | string  | semp/sso.sex            | SExtractor configuration file for source detection in the survey images.  |
|                       |         |                         | For details, see :ref:`sextractor_section`.                               |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `SEX_PARAMS`          | string  | semp/sso.param          | SExtractor output parameter for source detection in the survey images.    |
|                       |         |                         | For details, see :ref:`sextractor_section`.                               |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `SEX_FILTER`          | string  |semp/gauss_2.5_5x5 .conv | SExtractor convolution filter file for source detection in the survey     |
|                       |         |                         | images. For details, see :ref:`sextractor_section` and the                |
|                       |         |                         | `Guide to SExtractor`_.                                                   |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `SEX_NNW`             | string  | semp/sso.nnw            | SExtractor neural network for galaxy-star differentiation. For details,   |
|                       |         |                         | see :ref:`sextractor_section` and the `Guide to SExtractor`_.             |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `SCI_EXTENSION`       | integer | 1 |  2 | 1,2            | Index of science extension of FITS images. For details, see               |
|                       |         |                         | :ref:`sextractor_section`.                                                |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `SCAMP_CONFIG`        | string  | semp/sso.scamp          | SCAMP configuration file to link source detections at different epochs,   |
|                       |         |                         | see :ref:`scamp_section`.                                                 |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `SWARP_CONFIG`        | string  | semp/sso.swarp          | SWARP configuration file for creation of cutout images of SSO candidates, |
|                       |         |                         | see :ref:`swarp_section`.                                                 |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `FILTER_DETEC`        | bool    | True | False            | See  :ref:`filter_section`                                                |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `DETECTIONS`          | integer |  1,2 |  1,2,3,4 | 1,5   |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `FILTER_PM`           | bool    |   True | False          |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `PM_LOW`              | float   |     0.                  |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `PM_UP`               | float   |     200.                |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `PM_SNR`              | float   |      20.                |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `FILTER_PIXEL`        | bool    |   True | False          |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `DELTA_PIXEL`         | float   |      2.                 |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `FILTER_MOTION`       | bool    |    True | False         |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `IDENTIFY_OUTLIER`    | bool    |    True | False         |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `OUTLIER_THRESHOLD`   | float   |     2.                  |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `R_SQU_M`             | float   |     0.95                |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `FILTER_TRAIL`        | bool    |      True | False       |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `RATIO`               | float   |      0.25               |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `FILTER_T_DIST`       | bool    |     True | False        |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `SIGMA`               | float   |         2.              |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `FILTER_STAR_REGIONS` | bool    |      True | False       |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `DISTANCE`            | float   |        300.             |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `HYGCAT`              | string  | semp/hygdata_v3.csv     |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `CROSSMATCH_SKYBOT`   | bool    |     True | False        |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `CROSSMATCH_RADIUS`   | float   |        10.              |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `OBSERVATORY_CODE`    | string  |        500              |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `FOV_DIMENSIONS`      | string  |       1x1.5             |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `EXTRACT_CUTOUTS`     | bool    |     True | False        |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `CUTOUT_SIZE`         | integer |        256              |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `FIXED_APER_MAGS`     | bool    |    True | False         |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
| `REFERENCE_FILTER`    | string  |         gSDSS           |                                                                           |
+-----------------------+---------+-------------------------+---------------------------------------------------------------------------+

The configuration file can be formatted with tabs and spaces. Comments are marked with `#`. Lines beginning with # or newline characters are ignored.

.. note:: The pipeline script first checks if the `-c` flag is pointing to a configuration file. If not, a file called `pipeline_settings.ssos` is looked for in the current working directory. If no file is found, the hard-coded default values are used. Any parameter can be overwritten temporarily by using the appropriate flag, see :ref:`Command-Line API <Command-Line API>`.


Survey-specific changes
=======================

It is highly unlikely that the pipeline will give you the optimum result (clean and complete sample of SSOs) right out-of-the-box. You likely have to adjust the following files and parameters before running it the first time, mostly by setting them to the appropriate FITS header keywords of your images:



``sso.sex``

* SATUR_KEY

* GAIN_KEY

* SEEING_FWHM

* MAG_ZEROPOINT


``semp/sso.scamp``

* ASTRINSTRU_KEY

* ASTRACCURACY_KEY

* PHOTINSTRU_KEY

* MAGZERO_KEY

* EXPOTIME_KEY

* AIRMASS_KEY

* EXTINCT_KEY

* PHOTOMFLAG_KEY


``semp/sso.swarp``

* GAIN_KEYWORD



``pipeline_settings.ssos``

* SEX_CONFIG

* SEX_PARAMS

* SEX_FILTER

* SEX_NNW

* SCAMP_CONFIG

* SWARP_CONFIG

* HYGCAT

* OBSERVATORY CODE

* FOV SIZE


Furthermore, the implementation of weight images for the SExtractor runs can be beneficial.

After these initial changes, you should experiment with the different SExtractor, SCAMP, and pipeline settings, adjusting e.g. the filter chain parameters. A good way to fine-tune is to pick a test field with several SSOs and run the pipeline with different configurations. The cutout images will tell you what types of artifacts are remaining and whether you accidentally filtered out SSOs by restricting the candidate filters too much.


Input Files
===========

Input files are the survey images, passed to the pipeline in one directory:

.. code-block :: bash

    $ ssos path/to/images

The images **must** have a ``.fits`` file ending to be recongized by the script.

Output Files
============

The script creates several directories in the target directory (CWD by default). The cats directory contains the SExtractor and SCAMP catalogues and the final output ssos.csv. For every SExtractor catalog, there is also one .ahead file with the same filename. This file contains the observation date as MJD-OBS keyword, which is required for the subsequent run of SCAMP. SCAMP looks for extensions of catalog headers in .ahead files.

The cutouts directory contains the cutouts made by SWarp of the SSO detections. In the logs directory, you can find the log file of the pipeline, following the naming scheme sso_$CURRENTDATETIME.log. The skybot directory stores the SkyBoT query results.

To judge the results of the pipeline, first go into the logfile. It looks like this:

.. code-block :: bash

    $ cat sso_20180720103807.log

    --- SSO Recovery Pipeline ---   --- 10:32:04 2018/07/20 ---

    SEX_CONFIG          jvar_semp/jvar.sex
    SEX_PARAMS          jvar_semp/jvar.param
    SEX_NNW             jvar_semp/default.nnw
    SEX_FILTER          jvar_semp/gauss_2.5_5x5.conv
    SCAMP_CONFIG        jvar_semp/jvar.scamp
    SWARP_CONFIG        jvar_semp/jvar.swarp
    FILTER_DETEC        True
    DETECTIONS          123
    FILTER_PM           True
    PM_LOW              0
    PM_UP               225
    PM_SN               20
    FILTER_MOTION       True
    R_SQU_M             0.95
    FILTER_TRAIL        False
    RATIO               0.25
    FILTER_T_DIST       False
    SIGMA               2
    FILTER_STAR_REGIONS True
    MAGNITUDE           30
    DISTANCE            300
    HYGCAT              jvar_semp/hygdata_v3.csv
    CROSSMATCH_SKYBOT   True
    CROSSMATCH_RADIUS   10
    EXTRACT_CUTOUTS     True

    21 Exposures  |  JVAR00596epoch1  |  35.92deg Ecliptic Latitude

    Running SExtractor..  [---------------------]
    Running SCAMP..

     --- Starting Filter Pipeline ---

    All Sources         15326
    # of Detections     7267
    Proper Motion       26
    Linear Motion       2
    Star Catalog        2

    Cross-matching SSO candidates with SkyBoT.. 0 matches found
    Extracting cutouts with SWARP.. Done.

    All done!    |    2 SSOs found in 15326 Sources    |    The analysis ran in 364 seconds

    Output File: /tmp/cats/ssos.csv

In case an SSO was detected, you should look at the cats/ssos.csv file and the cutouts to verify the detection. In cats/ssos.csv, you can also find basic SkyBoT parameters, if the object was successfully matched. For more detailed information on the possible match, look into the SkyBoT queries in skybot.

The final database contains the following columns

.. code-block:: bash

    SOURCE_NUMBER - SCAMP groups single detections into one source by giving them the same SOURCE_NUMBER
    CATALOG_NUMBER - Number of SExtractor catalog containing this source detections
    RA - Right Asecension of source in degree
    DEC - Declincation of source in degree
    EPOCH - Beginning of observation in decimalyear
    MAG - Magnitude
    MAGERR - Error of magnitude as derived by SExtractor
    PM - Proper motion of source in "/h
    PMERR - Error on proper motion
    PMRA - Proper motion in RA in "/h
    PMRA_ERR - Error
    PMDEC - Proper motion in Dec in "/h
    PMDEC_ERR - Error
    MID_EXP_MJD - Mid-exposure time in MJD
    DATE-OBS_EXP - Beginning of observation in ISOT format
    EXPTIME_EXP - Exposure time
    OBJECT_EXP - Object ID of J_VAR field
    FILTER_EXP - Name of filter that the field was imaged in
    RA_EXP - Center right ascension coordinate of expsoure
    DEC_EXP - Center declination coordinate of expsoure
    FILENAME_EXP - Filename of input image that the SSO was detected in
    SKYBOT_NAME - If matched, name of matching SSO
    SKYBOT_CLASS - Class of SkyBoT match
    SKYBOT_MAG - Predicted magnitude
    SKYBOT_RA - Predicted RA
    SKYBOT_DEC - Predicted Dec
    SKYBOT_PMRA - Predicted PM in RA
    SKYBOT_PMDEC - Predicted PM in Dec
    SKYBOT_NUMBER - SkyBoT match number
    X_IMAGE - position of source in exposure in X in pixel
    Y_IMAGE - position of source in exposure in Y in pixel
    AWIN_IMAGE - Semi-major axis of fitted source ellipse in pixel
    ERRA_IMAGE - Error
    BWIN_IMAGE - Semi-minor axis of fitted source ellipse in pixel
    ERRB_IMAGE - Error
    THETAWIN_IMAGE - Angle of source, see SExtractor guide
    ERRTHETA_IMAGE - Error
    ERRA_WORLD - Error of AWIN in degree
    ERRB_WORLD - Error of BWIN in degree
    ERRTHETA_WORLD - Error of THETAWIN in degree
    ERRX2_WORLD - Variance of RA in degree
    ERRY2_WORLD - Variance of Dec in degree
    ERRXY_WORLD - Covariance of RA/Dec in degree
    FLAGS_EXTRACTION - SCAMP parameter
    FLAGS_SCAMP - SCAMP parameter
    FLAGS_IMA - SCAMP parameter

More information on these parameters can be found in the `SExtractor <https://readthedocs.org/projects/sextractor/>`_ and `SCAMP <https://scamp.readthedocs.io/en/latest/>`_ manuals.


.. _clapi:

Command-Line API
================

Again, the command-line API is heavily inspired by the SExtractor and SCAMP softwares. The following help is printed when the pipeline is called without arguments at all or with the ``-h`` or ``--help`` flag set:

.. code-block:: bash

    $ ssos --help
    usage: ssos [-h] [-c SET_UP] [-t TARGET] [-l LOG] [-v] [--sex] [--scamp]
                [--swarp] [-FILTER_DETEC bool] [-FILTER_PM bool]
                [-FILTER_PIXEL bool] [-FILTER_MOTION bool]
                [-IDENTIFY_OUTLIERS bool] [-FILTER_TRAIL bool]
                [-FILTER_T_DIST bool] [-FILTER_STAR_REGIONS bool]
                [-CROSSMATCH_SKYBOT bool] [-EXTRACT_CUTOUTS bool]
                [-FIXED_APER_MAGS bool] [-SEX_CONFIG value] [-SEX_PARAMS value]
                [-SEX_FILTER value] [-SEX_NNW value] [-SCAMP_CONFIG value]
                [-SWARP_CONFIG value] [-SCI_EXTENSION value] [-DETECTIONS value]
                [-PM_LOW value] [-PM_UP value] [-PM_SNR value]
                [-DELTA_PIXEL value] [-OUTLIER_THRESHOLD value] [-R_SQU_M value]
                [-R_SQU_T value] [-SIGMA value] [-DISTANCE value] [-HYGCAT value]
                [-CROSSMATCH_RADIUS value] [-CUTOUT_SIZE value]
                [-FIXED_UPER_MAG value] [-OBSERVATORY_CODE value]
                [-FOV_DIMENSIONS value]
                fields [fields ...]

    Pipeline to search for Solar System objects in wide-field imaging surveys

    positional arguments:
      fields                Path to directory of field to be searched

    optional arguments:
      -h, --help            show this help message and exit
      -c SET_UP, --config SET_UP
                            Path to configuration file
      -t TARGET, --target TARGET
                            Target directory to save fits files. If no target
                            given, writing to CWD
      -l LOG, --log LOG     Set the logging level. Validarguments are DEBUG, INFO,
                            WARNING, ERROR, CRITICAl.
      -v, --verbose         Print logging to console
      --sex                 Force SExtractor runs
      --scamp               Force SCAMP runs
      --swarp               Force SWARP runs

    Filter Settings:
      -FILTER_DETEC bool    Override FILTER_DETEC setting. Must be True or False.
      -FILTER_PM bool       Override FILTER_PM setting. Must be True or False.
      -FILTER_PIXEL bool    Override FILTER_PIXEL setting. Must be True or False.
      -FILTER_MOTION bool   Override FILTER_MOTION setting. Must be True or False.
      -IDENTIFY_OUTLIERS bool
                            Override IDENTIFY_OUTLIERS setting. Must be True or
                            False.
      -FILTER_TRAIL bool    Override FILTER_TRAIL setting. Must be True or False.
      -FILTER_T_DIST bool   Override FILTER_T_DIST setting. Must be True or False.
      -FILTER_STAR_REGIONS bool
                            Override FILTER_STAR_REGIONS setting. Must be True or
                            False.
      -CROSSMATCH_SKYBOT bool
                            Override CROSSMATCH_SKYBOT setting. Must be True or
                            False.
      -EXTRACT_CUTOUTS bool
                            Override EXTRACT_CUTOUTS setting. Must be True or
                            False.
      -FIXED_APER_MAGS bool
                            Override FIXED_APER_MAGS setting. Must be True or
                            False.
      -SEX_CONFIG value     Override SEX_CONFIG setting.
      -SEX_PARAMS value     Override SEX_PARAMS setting.
      -SEX_FILTER value     Override SEX_FILTER setting.
      -SEX_NNW value        Override SEX_NNW setting.
      -SCAMP_CONFIG value   Override SCAMP_CONFIG setting.
      -SWARP_CONFIG value   Override SWARP_CONFIG setting.
      -SCI_EXTENSION value  Override SCI_EXTENSION setting.
      -DETECTIONS value     Override DETECTIONS setting.
      -PM_LOW value         Override PM_LOW setting.
      -PM_UP value          Override PM_UP setting.
      -PM_SNR value         Override PM_SNR setting.
      -DELTA_PIXEL value    Override DELTA_PIXEL setting.
      -OUTLIER_THRESHOLD value
                            Override OUTLIER_THRESHOLD setting.
      -R_SQU_M value        Override R_SQU_M setting.
      -R_SQU_T value        Override R_SQU_T setting.
      -SIGMA value          Override SIGMA setting.
      -DISTANCE value       Override DISTANCE setting.
      -HYGCAT value         Override HYGCAT setting.
      -CROSSMATCH_RADIUS value
                            Override CROSSMATCH_RADIUS setting.
      -CUTOUT_SIZE value    Override CUTOUT_SIZE setting.
      -FIXED_UPER_MAG value
                            Override FIXED_UPER_MAG setting.
      -OBSERVATORY_CODE value
                            Override OBSERVATORY_CODE setting.
      -FOV_DIMENSIONS value
                            Override FOV_DIMENSIONS setting.


The verbose flag ``-v`` is recommended at first to understand what the script does. It prints the logging output to console (as well as to the logfile).