######
How To
######

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


::

  target directory
  │
  └───cats
  │      exposure1.cat
  │      exposure1.head
  │      exposure1.ahead  *optional*
  │      exposure2.cat
  │      exposure2.head
  │      exposure2.ahead  *optional*
  │      ..
  │      full_1.cat
  │      merged_1.cat
  │      scamp.xml
  │      ssos.csv
  │
  └───cutouts
  │      SOURCE1_CATALOG1.fits
  │      SOURCE1_CATALOG2.fits
  │      SOURCE2_CATALOG1.fits
  │      ..
  │
  └───logs
  │      sso_$DATETIME.log
  │
  └───skybot
  │      skybot_query_string1.xml
  │      skybot_query_string2.xml
  │      ..
  │
  └───weights
  │      exposure1.weight  *optional*
  │      expsoure2.weight  *optional*
  │      ..


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


