###########################
Usage & Code Implementation
###########################

Input Files
===========

The `pipeline_settings.sso` file looks like this:

::

    # ------
    # SExtractor File Settings

    SEX_CONFIG           path/to/default.sex              # Path to SExtractor config file
    SEX_PARAMS           path/to/default.param            # Path to SExtractor param list
    SEX_FILTER           path/to/gauss_2.5_5x5.conv       # Path to SExtractor filter file
    SEX_NNW              path/to/default.nnw              # Path to SExtractor neural network file
    SCI_EXTENSION        1                                # Science Extension of FITS file

    # ------
    # SCAMP File Settings

    SCAMP_CONFIG         path/to/default.scamp            # Path to SCAMP config file

    # ------
    # SWARP File Settings

    SWARP_CONFIG         path/to/default.swarp            # Path to SWARP config file

    # ------
    # SSO Filter Settings


    FILTER_DETEC         True                             # Filter by number of detections
    DETECTIONS           123                              # Number of detections to filter

    FILTER_PM            True                             # Filter by proper motion range
    PM_LOW               0                                # Lower proper motion limit / "/h
    PM_UP                2000                             # Upper proper motion limit / "/h
    PM_SNR               20                               # Lower limit of proper motion SNR

    FILTER_PIXEL         True                             # Filter by pixel margin
    DELTA_PIXEL          2                                # Lower limit of pixel margin

    FILTER_MOTION        True                             # Filter by linear motion
    IDENTIFY_OUTLIER     True                             # Flag outlier and fit subgroups
    OUTLIER_THRESHOLD    2                                # Threshold in DeltaEpoch / MAD
    R_SQU_M              0.95                             # Minimum goodness-of-fit

    FILTER_TRAIL         False                            # Filter by constant trail size
    RATIO                0.25                             # Minimum ratio of mean/outlier

    FILTER_T_DIST        False                            # Filter by catalog trail distribution
    SIGMA                2                                # Upper limit of trail size in stddev

    FILTER_STAR_REGIONS  True                             # Filter by distance to bright star
    DISTANCE             300                              # Minimum distance to star in "

    CROSSMATCH_SKYBOT    True                             # Cross-match sources with SkyBoT
    CROSSMATCH_RADIUS    10                               # Upper distance in " to count as match
    OBSERVATORY_CODE     500                              # IAU Observatory code for SkyBoT query
    FOV_DIMENSIONS       0.55x0.55                        # Size of query region / deg

    EXTRACT_CUTOUTS      True                             # Save cutouts of source detections

    HYGCAT               path/to/hygdata_v3.csv           # Path to star catalog



Command-Line API
================

-a            command-line option "a"
-b file       options can have arguments
              and long descriptions
--long        options can be long also
--input=file  long options can also have
              arguments
/V            DOS/VMS-style options too

SExtractor and SCAMP wrapper
============================


Filter Chain
============

Output Files
============

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
  │      exposure1_weight.fits  *optional*
  │      expsoure2_weight.fits  *optional*
  │      ..
