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

or by using the setup script:

.. code-block:: bash

    $ [sudo] python setup.py install

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

The default ``pipeline_settings.ssos`` file can be `found here <https://github.com/maxmahlke/SSO_Pipeline/blob/master/sso_pipeline/pipeline_settings.ssos>`_. It is a plain ASCII, designed very similar to the configuration files of SExtractor and SCAMP in order to make the user feel right at home. Short descriptions and expected values of the parameters are below, for more detailed descriptions refer to the `Pipeline <pipeline.rst>`_ and `Implementation <implementation.rst>`_ pages.

.. _Guide to SExtractor: http://astroa.physics.metu.edu.tr/MANUALS/sextractor/Guide2source_extractor.pdf

.. _IAU Observatory Code: http://vo.imcce.fr/webservices/data/displayIAUObsCodes.php

.. _SkyBoT: http://vo.imcce.fr/webservices/skybot/?conesearch


.. table::
    :align: center

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
    | `WEIGHT_IMAGES`       | bool    | False | /tmp/weights    | Absolute path to weight images for SExtractor run. If False, SExtractor   |
    |                       |         |                         | runs with settings according to ``ssos.sex`` file.                        |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `SCAMP_CONFIG`        | string  | semp/sso.scamp          | SCAMP configuration file to link source detections at different epochs,   |
    |                       |         |                         | see :ref:`scamp_section`.                                                 |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `SWARP_CONFIG`        | string  | semp/sso.swarp          | SWARP configuration file for creation of cutout images of SSO candidates, |
    |                       |         |                         | see :ref:`optional`.                                                      |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `FILTER_DETEC`        | bool    | True | False            | Turn filter based on number of detections on or off.                      |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `DETECTIONS`          | integer |  1,2 |  1,2,3,4 | 1,5   | Sources with this number of detections are rejected.                      |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `FILTER_PM`           | bool    |   True | False          | Turn filter based on proper motion values on or off.                      |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `PM_LOW`              | float   |     0.                  | Lower limit on proper motion of sources. See :ref:`filter_section`.       |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `PM_UP`               | float   |     200.                | Upper limit on proper motion of sources. See :ref:`filter_section`.       |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `PM_SNR`              | float   |      20.                | Lower limit on signal-to-noise ratio of proper motion of sources.         |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `FILTER_PIXEL`        | bool    |   True | False          | Turn filter based on pixel positions on or off. See :ref:`filter_section`.|
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `DELTA_PIXEL`         | float   |      2.                 | Minimum number of pixel the centre position of the source has to shift by |
    |                       |         |                         | over all exposures in X and Y. See :ref:`filter_section`.                 |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `FILTER_MOTION`       | bool    |    True | False         | Turn filter based on linearity of motion on or off.                       |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `IDENTIFY_OUTLIER`    | bool    |    True | False         | Identify outliers in epoch-space and treat their motion separately.       |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `OUTLIER_THRESHOLD`   | float   |     2.                  | Threshold in Median Absolute Deviations for identification of outlier.    |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `R_SQU_M`             | float   |     0.95                | Lower limit of R-Squared goodness-of-fit parameter for linear motion fit. |
    |                       |         |                         | Must be between 0 and 1. See :ref:`filter_section`.                       |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `FILTER_TRAIL`        | bool    |      True | False       | Turn filter based on constant trail parameters on or off.                 |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `RATIO`               | float   |      0.25               | Lower limit on the ratio of the error on the weighted mean to the standard|
    |                       |         |                         | deviation of the source ellipse parameters. See :ref:`filter_section`     |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `FILTER_T_DIST`       | bool    |     True | False        | Turn filter based on distribution of trail sizes in image on or off.      |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `SIGMA`               | float   |         2.              | Upper limit in standard deviation to find outlier in source ellipse       |
    |                       |         |                         | parameters. See :ref:`filter_section`.                                    |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `FILTER_STAR_REGIONS` | bool    |      True | False       | Turn filter based on source distance to bright stars on or off.           |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `DISTANCE`            | float   |        300.             | Minimum distance of source to bright star in star catalogue in arcsecond. |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `HYGCAT`              | string  | semp/hygdata_v3.csv     | Absolute path to `HYG <http://www.astronexus.com/hyg>`_ star catalogue.   |
    |                       |         |                         | See :ref:`filter_section`.                                                |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `CROSSMATCH_SKYBOT`   | bool    |     True | False        | Turn cross-matching with SkyBoT database on or off. See :ref:`optional`.  |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `CROSSMATCH_RADIUS`   | float   |        10.              | Upper limit of distance between source candidate and SkyBoT source to     |
    |                       |         |                         | be considered a match, in arcsecond. See :ref:`optional`.                 |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `OBSERVATORY_CODE`    | string  |        500              | `IAU Observatory Code`_                                                   |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `FOV_DIMENSIONS`      | string  |       1x1.5             | Dimensions of exposure field-of-view in degrees, see `SkyBoT`_.           |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `EXTRACT_CUTOUTS`     | bool    |     True | False        | Turn cutout extraction with SWARP on or off. See :ref:`optional`.         |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `CUTOUT_SIZE`         | integer |        256              | Size of cutouts in pixel, each dimension, see :ref:`optional`.            |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `FIXED_APER_MAGS`     | bool    |    True | False         | Compute fixed aperture magnitudes for colours. See :ref:`optional`.       |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+
    | `REFERENCE_FILTER`    | string  |         gSDSS           | Filter to use as reference in SExtractor dual-image mode runs. Value has  |
    |                       |         |                         | to correspond to `FILTER` keyword in FITS header. See :ref:`optional`.    |
    +-----------------------+---------+-------------------------+---------------------------------------------------------------------------+

The configuration file can be formatted with tabs and spaces. Comments are marked with `#`. Lines beginning with # or newline characters are ignored.

.. note:: The pipeline script first checks if the `-c` flag is pointing to a configuration file. If not, it looks for a file called `pipeline_settings.ssos` in the current working directory. If no file is found, the hard-coded default values are used. Any parameter can be overwritten temporarily by using the appropriate flag, see :ref:`Command-Line API <Command-Line API>`.


Survey-specific changes
=======================

It is highly unlikely that the pipeline will give you the optimum result (clean and complete sample of SSOs) right out-of-the-box. You likely have to adjust the following files and parameters before running it the first time, mostly by setting them to the appropriate FITS header keywords of your images:



``ssos.sex``

    - `SATUR_KEY`

    - `GAIN_KEY`

    - `SEEING_FWHM`

    - `MAG_ZEROPOINT`


``semp/ssos.scamp``

    - `ASTRINSTRU_KEY`

    - `ASTRACCURACY_KEY`

    - `PHOTINSTRU_KEY`

    - `MAGZERO_KEY`

    - `EXPOTIME_KEY`

    - `AIRMASS_KEY`

    - `EXTINCT_KEY`

    - `PHOTOMFLAG_KEY`


``semp/ssos.swarp``

    - `GAIN_KEYWORD`



``pipeline_settings.ssos``

    - `SEX_CONFIG`

    - `SEX_PARAMS`

    - `SEX_FILTER`

    - `SEX_NNW`

    - `SCAMP_CONFIG`

    - `SWARP_CONFIG`

    - `HYGCAT`

    - `OBSERVATORY CODE`

    - `FOV SIZE`


After these initial changes, you should experiment with the different SExtractor, SCAMP, and pipeline settings, adjusting e.g. the filter chain parameters. A good way to fine-tune is to pick a test field with several SSOs and run the pipeline with different configurations. The cutout images will tell you what types of artifacts are remaining and whether you accidentally filtered out SSOs by restricting the candidate filters too much.


Input Files
===========

Input files are the survey images, passed to the pipeline in one directory:

.. code-block :: bash

    $ ssos path/to/images

The images **must** have a ``.fits`` file ending to be recognized by the script.

If weight images shall be passed to the SExtractor runs, specify the directory containing the weight images using the `WEIGHT_IMAGES` parameter in the ``pipeline_settings.ssos``. The weight images need to have the same filename as the exposures they shall be applied on, but with a ``_SCI_EXTENSION.weight`` file extension replacing the ``.fits``, where `SCI_EXTENSION` has to be replaced by the appropriate extension integer. [#]_

Output Files
============

The script creates several directories in the target directory (CWD by default). The cats directory contains the SExtractor and SCAMP catalogues and the final output ssos.csv. For every SExtractor catalogue, there is also one .ahead file with the same filename. This file contains the observation date as MJD-OBS keyword, which is required for the subsequent run of SCAMP. SCAMP looks for extensions of catalogue headers in .ahead files.

The cutouts directory contains the cutouts made by SWarp of the SSO detections. In the logs directory, you can find the log file of the pipeline, following the naming scheme sso_$CURRENTDATETIME.log. The skybot directory stores the SkyBoT query results.

To judge the results of the pipeline, first go into the logfile. It looks like this:

.. code-block :: bash

    $ cat /tmp/logs/sso_20180906142433.log

        --- SSO Recovery Pipeline ---       --- 14:24:33 2018/09/06 ---

    SEX_CONFIG          semp/ssos.sex
    SEX_PARAMS          semp/ssos.param
    SEX_NNW             semp/default.nnw
    SEX_FILTER          semp/gauss_2.5_5x5.conv
    SCI_EXTENSION       0
    SCAMP_CONFIG        semp/ssos.scamp
    SWARP_CONFIG        semp/ssos.swarp
    FILTER_DETEC        True
    DETECTIONS          1,2,3
    FILTER_PM           True
    PM_LOW              0
    PM_UP               225
    PM_SNR              20
    FILTER_PIXEL        True
    DELTA_PIXEL         2
    FILTER_MOTION       True
    IDENTIFY_OUTLIER    True
    OUTLIER_THRESHOLD   2
    R_SQU_M             0.95
    FILTER_TRAIL        False
    RATIO               0.25
    FILTER_T_DIST       False
    SIGMA               2
    FILTER_STAR_REGIONS True
    DISTANCE            300
    HYGCAT              semp/hygdata_v3.csv
    CROSSMATCH_SKYBOT   True
    CROSSMATCH_RADIUS   10
    OBSERVATORY_CODE    500
    FOV_DIMENSIONS      1.5x1.5
    EXTRACT_CUTOUTS     True
    FIXED_APER_MAGS     False
    REFERENCE_FILTER    gSDSS
    WEIGHT_IMAGES       False
    CUTOUT_SIZE         256

    21 Exposures    |   epoch1 |   35.92deg Ecliptic Latitude


    Running SExtractor..
    Running SCAMP.. Done.

     --- Starting Filter pipeline ---

    All Sources          15323
    FILTER_DETEC         7744
    FILTER_PM            30
    FILTER_PIXEL         30
    FILTER_MOTION        8
    FILTER_STAR_REGIONS  7

    Cross-matching SSO candidates with SkyBoT.. 2 matches found.

    Extracting cutouts with SWARP.. Done.

    All done!   |   7 SSO candidates found  |   The analysis ran in 4 seconds

    Output File: /tmp/cats/ssos.csv
    Log File:   /tmp/logs/sso_20180906142433.log

In case an SSO was detected, you should look at the ``cats/ssos.csv`` file and the cutouts to verify the detection. In ``cats/ssos.csv``, you can also find basic SkyBoT parameters, if the object was successfully matched. For more detailed information on the possible match, look into the SkyBoT queries in ``skybot/``.

The final database contains the following columns

.. code-block:: bash

    SOURCE_NUMBER - SCAMP groups detections into sources by giving them the same SOURCE_NUMBER
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
    RA_EXP - Center right ascension coordinate of exposure
    DEC_EXP - Center declination coordinate of exposure
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
    FLAS_SSOS - Flags added by ssos pipeline

More information on these parameters can be found in the `SExtractor <https://readthedocs.org/projects/sextractor/>`_ and `SCAMP <https://scamp.readthedocs.io/en/latest/>`_ manuals.


.. _clapi:

Command-Line API
================

Again, the command-line API is heavily inspired by the SExtractor and SCAMP softwares. The following help is printed when the pipeline is called without arguments or with the ``-h`` or ``--help`` flag set:

.. code-block:: bash

    $ ssos --help
    usage: ssos [-h] [-c SET_UP] [-t TARGET] [-l LOG] [-q] [--sex] [--scamp]
            [--swarp] [--skybot] [-FILTER_DETEC bool] [-FILTER_PM bool]
            [-FILTER_PIXEL bool] [-FILTER_MOTION bool]
            [-IDENTIFY_OUTLIERS bool] [-FILTER_TRAIL bool]
            [-FILTER_T_DIST bool] [-FILTER_STAR_REGIONS bool]
            [-CROSSMATCH_SKYBOT bool] [-EXTRACT_CUTOUTS bool]
            [-FIXED_APER_MAGS bool] [-SEX_CONFIG value] [-SEX_PARAMS value]
            [-SEX_FILTER value] [-SEX_NNW value] [-SCAMP_CONFIG value]
            [-SWARP_CONFIG value] [-SCI_EXTENSION value]
            [-WEIGHT_IMAGESDETECTIONS value] [-PM_LOW value] [-PM_UP value]
            [-PM_SNR value] [-DELTA_PIXEL value] [-OUTLIER_THRESHOLD value]
            [-R_SQU_M value] [-R_SQU_T value] [-SIGMA value] [-DISTANCE value]
            [-HYGCAT value] [-CROSSMATCH_RADIUS value] [-CUTOUT_SIZE value]
            [-REFERENCE_FILTER value] [-OBSERVATORY_CODE value]
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
      -l LOG, --log LOG     Set the logging level. Valid arguments are DEBUG,
                            INFO, WARNING, ERROR, CRITICAl.
      -q, --quiet           Suppress logging to console
      --sex                 Force SExtractor runs
      --scamp               Force SCAMP runs
      --swarp               Force SWARP runs
      --skybot              Force SkyBoT query

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
      -WEIGHT_IMAGESDETECTIONS value
                            Override WEIGHT_IMAGESDETECTIONS setting.
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
      -REFERENCE_FILTER value
                            Override REFERENCE_FILTER setting.
      -OBSERVATORY_CODE value
                            Override OBSERVATORY_CODE setting.
      -FOV_DIMENSIONS value
                            Override FOV_DIMENSIONS setting.

.. [#] Do not forget to change the `WEIGHT_TYPE` parameter in ``ssos.sex`` to activate the weight images, only supplying the path to the directory is not enough.
