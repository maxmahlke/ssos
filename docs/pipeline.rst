##################
The Pipeline Steps
##################

The three main steps of the pipeline are

* identifying all sources in all images using `Source Extractor <https://www.astromatic.net/software/sextractor>`_
* cross-matching all sources between all images using `SCAMP <https://www.astromatic.net/software/scamp>`_
* separating the SSOs from other sources such as stars, galaxies, artifacts, using a filter chain


.. _sextractor_section:

SExtractor
==========

SExtractor identifies sources in CCD images using the pixel values and outputs catalogue data. Parameters like the pixel positions, sky coordinates, instrumental magnitudes and more are retrieved for each source. Refer to `the official documentation <https://readthedocs.org/projects/sextractor/>`_ and the `Guide to SExtractor <http://astroa.physics.metu.edu.tr/MANUALS/sextractor/Guide2source_extractor.pdf>`_ for much better explanations.

SExtractor is highly configurable using a configuration file. The config file ``ssos.sex`` that is provided with this survey differs slightly from the default version ``default.sex``, specifically the deblending and photometry parameters were adjusted to better deblend SSOs close to other sources. The path to this config file has to be set in the ``pipeline_settings.ssos`` file using the `SEX_CONFIG` parameter. Likewise, `SEX_PARAMS` has to point to the output parameter file, by default ``semp/ssos.param``.

The SExtractor configuration as set in ``ssos.sex`` requires two more files, the convolution filter set via the `SEX_FILTER` parameter and the neural network for star-galaxy differentiation, set via the `SEX_NNW` keyword, located in ``semp/gauss_2.5_5x5.conv`` and ``semp/default.nnw`` respectively.

Unless the path parameters above are set using absolute paths, the pipeline will look for the files in the directory it is executed in.

Finally, using the `SCI_EXTENSION`, the user has to provided the index of the science extension of the FITS image. It is common case that data is provided in multi-extension FITS format, where besides the science data also weight images and other supplementary data is stored in the FITS file. Therefore, the user has to specify the science extension. Multiple extensions are allowed as well, if for example different CCD images were stored in the same file. The valid values are integers, separated by commas if multiple extensions should be analysed, e.g. **0**, **1**, or **1,2**.

If you are unsure which extension contains your image, you can trial run SExtractor with the following syntax and check the output catalogues:

.. code-block:: bash

    sex -c semp/ssos.sex image_file.fits[SCI_EXTENSION]


.. note::
    Searching multiple extensions at the same time only makes sense if the field-of-views overlap. Otherwise, running the pipeline on the extensions separately will yield better results.


The SEXtractor output catalogues will be called ``image_file_SCI_EXTENSION.cat`` [#]_.
After the SExtractor run, the input images are checked for the `MJD-OBS` header keyword.
If it does not exist, the `DATE-OBS keyword` is read in, converted to MJD, and saved in an ``.ahead`` file with the same filename structure as the SExtractor catalog. This additional header file is important for the subsequent SCAMP run.

.. _scamp_section:

SCAMP
=====

SCAMP provides the astrometric solution for the pipeline: The SExtractor catalogues from the first step sharing the same field-of-view (FoV) are matched against each other using a reference catalogue, and the source coordinates are solved for translation, rotation, and distortion of the original images. This way, sources detected in several images over different epochs can be traced from one image to another, their detections are linked. Again, find a better explanation in the `official documentation <https://scamp.readthedocs.io/en/latest/>`_.

Much like SExtractor, SCAMP is highly configurable and the supplied configuration file ``ssos.scamp`` has to be linked to using the `SCAMP_CONFIG` in the ``pipeline_settings.ssos`` file.

Once SCAMP has matched the SExtractor catalogues, it creates among others two catalogues, the ``full_cat_1.cat`` and the ``merged_cat_1.cat``. The former contains all source detections of all images, given positions in pixel and sky coordinates, fluxes, etc., while the latter lists the properties of the merged (linked) detections, i.e. proper motion and other properties of all sources which were traced over several images. Both catalogues will be used in the subsequent analysis, specifically the full catalogue, as it holds the important information on the movement of the sources over time.

.. note::

    The SExtractor and SCAMP runs are the computationally most challenging parts of this pipeline and therefore the bottlenecks in execution time. To allow for quick pipeline runs in order to find the optimal settings, the script checks for the existence of the output catalogues before running the software. If the catalogues already exist, these steps are skipped. This behaviour can be overruled by setting the ``--sex``, ``--scamp``, and ``--swarp`` flags in the pipeline call.

.. _filter_section:

Filter Chain
============

All filter steps are optional and can be turned on/off and adjusted via the ``pipeline_settings.ssos`` configuration file.

Filter by Number of Detections
------------------------------
Setting: `FILTER_DETECT`  |  Parameters: `DETECTIONS`

All sources with a number of detections equal to the numbers specified in the `DETECTIONS` parameter are removed. By default, `DETECTIONS` is `123`, removing all sources with fewer than 4 detections. Removing sources with only 1 or 2 detections is always recommended, as their motion cannot be judged. It is, however, not enforced by the pipeline.

As a rule of thumb, artifacts such as random CR associations tend to have fewer detections than SSOs, which in turn have fewer detections than stars and galaxies. Increasing the required number of detections is an effective way to clean the sample, though at the cost of possibly losing faint SSOs and SSOs in the edge regions of the images.

Filter by Proper Motion Range
-----------------------------
Setting: `FILTER_PM`  |  Parameters: `PM_LOW`, `PM_UP`, `PM_SNR`

All sources with proper motions lower than `PM_LOW` and larger than `PM_UP` are rejected. Furthermore, the lower limit of the signal-to-noise ratio (SNR) of the proper motion measurement performed by SCAMP can be set using `PM_SNR`. SCAMP performs a linear fit of the source coordinates over time to determine the proper motion. Large uncertainties signal sources which do not move with constant proper motions, as expected from SSOs.

Effectively, the SNR lower limit introduces a lower limit on the proper motion as well. If the proper motion of an SSO over the exposure time is within order of the seeing conditions, it will exhibit large fluctuations in position and therefore be assigned a large error in the proper motion measurement by SCAMP.

Filter by Bad Pixel
-------------------
Setting: `FILTER_PIXEL`  |  Parameters: `DELTA_PIXEL`

If all detections of a single source fall within the same pixel `DELTA_PIXEL` range (both `XWIN_IMAGE` and `YWIN_IMAGE` parameters), the source is rejected.
Bad CCD pixel can be falsely interpreted as sources by SExtractor and SCAMP. Due to the dithering patterns, they appear to move perfectly linear and with a constant proper motion. SExtractor parameters like `DETECT_MINAREA` can be used to clean these sources, but increasing the minimum pixel area per source can also reject faint SSOs. The filter chain therefore also offers this rudimentary bad pixel rejection.

Filter by Motion
----------------
Setting: `FILTER_MOTION`  |  Parameters: `IDENTIFY_OUTLIER`, `OUTLIER_THRESHOLD`, `R_SQU_M`

The motion filter is the most effective and strictest filter. A linear fit is applied to both the RA and the DEC coordinates against observation epochs. If the `R^2` goodness-of-fit parameter of both fits is equal or larger than the user-defined `R_SQU_M` parameter (0 <= R^2 <= 1), the source is accepted. If either fit is not within the limit, the source is rejected. If `R_SQU_M` is between 0.95 and 1, this imposes very strict rules on the motion. Slow moving SSOs (proper motion in the order of seeing) might be missed if `R_SQU_M` is too big, while a lower setting will increase the number of artifacts surviving the pipeline.

The filter is effective in sorting out stars and galaxies from the sample, as they are stationary over the period of time, and the centroid position found by SExtractor will randomly fluctuate within the order of the seeing.

Problems arise when the observations span multiple hours or nights. If the survey images for example cover one area of the sky for the whole night with 50 exposures, it may occur that an SSO is observed in the first and the last 5 exposures. Such a long baseline with no observations in between will almost always yield a perfect linear fit. The same is true for sources randomly associated by stars, e.g. two stars close together or a star and several CRs. Again, the linear motion filter will be fooled by the large baseline of observations.
To tackle this problem, the `IDENTIFY_OUTLIER` option was introduced. If `True`, the motion filter starts by detecting outliers in epoch-space within the detections of one source. This is achieved using the **Median Absolute Deviation** (MAD) of the observation epochs *E*.

.. math::

   \mathrm{MAD} = \mathrm{median}(|E_{i} - \mathrm{median}(E)|)

This calculates the median duration between one observation and the median observation epoch. The median is not affected by outliers, therefore it can be used to identify jumps in the epochs. If the time difference between any two observations is larger than `MAD*OUTLIER_THRESHOLD`, the source detections are split into subgroups. If more than one of the jumps is found, the detections are split into several subgroups.
As long as the number of detections in each subgroup is larger or equal to the lower limit defined by the `DETECTIONS`, the detections within the subgroup are then checked for linear motion by the fitting procedure described above. If any subgroup fails the linear motion test, the source is discarded. If a subgroup has too few detections, it is only discarded if the other subgroup fails the linear motion test or if all other subgroups do not contain the sufficient amount of observations either.

All source detections which were identified as outliers in epoch space get +1 added to their `FLAGS_SSOS` parameter. If a source contains "only outliers" (e.g. two pairs of two detections with a large gap in between), the source is removed.

.. todo::

    Add figures of fits, outliers, subgroup fits


Filter by Trail Consistency
---------------------------
Setting: `FILTER_TRAIL`  |  Parameters: `RATIO`

Assuming roughly constant exposure time and seeing conditions, the SSO trail in the images should not vary in size. This is expressed by the `AWIN_IMAGE` and `BWIN_IMAGE` SExtractor parameters, which are the semi-major- and semi-minor axes of the ellipse fitted to the source. Varying size parameters indicate an association of random sources (e.g. cosmic ray + star). This filter compares the standard deviation of both `AWIN_IMAGE` and `BWIN_IMAGE` of all detections of one source against the weighted average uncertainty,

.. math::

    \mathrm{\texttt{RATIO}} = \frac{\overline{\sigma_{w}} }{ \sigma_{x}}, \qquad x~\epsilon~\{\verb|AWIN_IMAGE|, \verb|BWIN_IMAGE|\}

.. math::

    \overline{\sigma_{w}} = \Big( \sqrt{\sum_i w_{i,x}} \Big)^{-1}

.. math::

    w_{x} = \frac{1}{\sqrt{var_{x}}} \qquad var~\epsilon~\{\verb|ERRAWIN_IMAGE|, \verb|ERRBWIN_IMAGE|\}

and removes sources which show standard deviations larger than the `RATIO` parameter allows for.


.. note::

    By default, this filter is disabled.


Filter by Trail Size Distribution
---------------------------------
Setting: `FILTER_T_DIST`  |  Parameters: `SIGMA`

This filter acts on the SExtractor source ellipse parameters `AWIN_IMAGE` and `BWIN_IMAGE`. The standard deviation of each of the semi-major- and semi-minor axes is calculated. Sources with size parameters larger than the mean plus `SIGMA` times the standard deviation are rejected. This filter was implemented against ghosts introduced by bright stars, which can perfectly imitate linear motion depending on the dithering pattern of the observations.

.. note::

    By default, this filter is disabled.


Filter by Star Region
---------------------
Setting: `FILTER_STAR_REGIONS`  |  Parameters: `DISTANCE`, `HYGCAT`

Bright stars tend to introduce numerous artifacts like refraction spikes and reflection ghosts into images. As the position of these artifacts depends on the camera geometry and pointing, they tend to follow the dithering pattern and display linear movement over all observation epochs. Sources close to bright stars therefore tend to contain a large fraction of these artifacts, and can be rejected with this filter. The `DISTANCE` parameter sets the radius around bright stars in arcsecond where all sources are cleared from. The `HYG database <http://www.astronexus.com/hyg>`_ is used to define the RA / DEC coordinate pairs of bright stars and is located in ``semp/hygdata_v3.csv``.

.. _optional:

Optional Analyses
=================

SkyBoT Cross-match
------------------
Setting: `CROSSMATCH_SKYBOT`  |  Parameters: `CROSSMATCH_RADIUS`, `OBSERVATORY_CODE`, `FOV_DIMENSIONS`

Query the `SkyBoT <http://vo.imcce.fr/webservices/skybot/?conesearch>`_ database for SSOs in the field-of-view defined by `FOV_DIMENSIONS` and the centre coordinates of each exposure for each observation epoch. The query result is saved as ``skybot/query_string.XML`` file. The positions of all SSO candidates are then compared to the predicted positions of known SSOs, and if a match is found within the `CROSSMATCH_RADIUS` (in arcsecond), the predicted SkyBoT parameters are added to the source metadata in the database.

The `FOV_DIMENSIONS` parameter has to be defined as described on the SkyBoT webpage, a string of format "YxZ", where Y and Z are the image dimensions (integer or floating value) in degree.


Cutout Extraction with SWARP
----------------------------
Setting: `EXTRACT_CUTOUTS`  |  Parameters: `SWARP_CONFIG`, `CUTOUT_SIZE`

Use SWARP to create cutout images with dimension size `CUTOUT_SIZE` in pixel of each SSO detection. The cutouts are saved in the format ``cutouts/SOURCE_NUMBER__CATALOG_NUMBER.fits``. Using e.g. `imagemagick <https://www.imagemagick.org/script/index.php>`_, these cutouts can be quickly turned into little movies for visual confirmation of their nature. The `SWARP_CONFIG` file is used to configure the cutout extraction.


Compute Fixed Aperture Magnitudes
---------------------------------
Setting: `FIXED_APER_MAG`  |  Parameters: `REFERENCE_FILTER`, `CUTOUT_SIZE`

To measure SSO colours, the magnitudes in different bands using fixed apertures has to be computed. In the mandatory SExtractor part of the pipeline, the magnitudes are measured with variable Kron-apertures. This step uses the cutout images of SSOs to apply SExtractor in dual-image mode: One exposure is used to detect objects and compute the apertures, whereas the other is used for flux measurements. The detection image should be the deepest exposure available for best results. This band can be chosen using the `REFERENCE_FILTER` parameter, which has to be equal to the `FILTER` keyword of the detection image.
As not all source candidates are necessarily observed in this band, the value can be set to multiple bands, separated by commas. The script will then prioritize the source detections according to the order specified in this value, e.g. `FILTER1,FILTER2,FILTER3`. The source detection which was chosen as reference detection this way gets flagged by adding 2 to the `FLAGS_SSOS` output parameter.

After the fixed aperture magnitudes are calculated, the columns `MAG_CI` and `MAGERR_CI` are added to the database.

If the cutout extraction with SWARP was set to False, the cutouts will be created in this step and saved to a temporary folder, which is deleted after the pipeline finishes.


Flags
=====

The `FLAGS_SSOS` parameter is used to highlight sources which pass the filter but might be artifacts. An example are sources with jumps (outliers) in their observation epochs, which fools the linear motion filter. The flag values are represented by powers of 2 and added together, allowing for multiple flags to be set at the same time. The flag values are:

    =============  =======================================
    Integer Value     Meaning
    -------------  ---------------------------------------
          1        Source detection is an outlier in EPOCH
          2        Source detection used as reference for
                   fixed aperture magnitude measurement
    =============  =======================================

.. [#] Appending the [SCI_EXTENSION] bit after .cat confuses the popular TOPCAT tool, so consistency in naming was neglected here.
