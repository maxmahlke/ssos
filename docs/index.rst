.. SSO Pipeline documentation master file, created on Thu Aug 16 17:02:11 2018

#####################
The ``ssos`` Pipeline
#####################

Recover astrometry and photometry of Solar System Objects (SSOs) in wide-field imaging surveys.

.. image:: imgs/sso_banner.png
    :align: center
    :scale: 40
    :alt: alternate text


Introduction to the SSO Pipeline
================================
The SSO Pipeline is a flexible tool to identify and characterize Solar System Objects in consecutive exposures of the sky. Best use cases are wide-field imaging surveys such as the `Kilo-Degree Survey <http://kids.strw.leidenuniv.nl/>`_, for which the pipeline was `originally designed <https://www.aanda.org/articles/aa/ref/2018/02/aa30924-17/aa30924-17.html>`_.

The pipeline is written in ``python 3``. The three main steps are

* identifying all sources in all images using `SExtractor <https://www.astromatic.net/software/sextractor>`_
* cross-matching all sources between all images using `SCAMP <https://www.astromatic.net/software/scamp>`_
* separating the SSOs from other sources such as stars, galaxies, artifacts, using a filter chain

Additional analyses can be done, e.g. cross-matching SSO candidates with the `SkyBoT <http://vo.imcce.fr/webservices/skybot/>`_ database or creating cutouts with `SWarp <https://www.astromatic.net/software/swarp>`_.

The strength of the pipeline lies in the large degree of flexibility. All steps can be adapted to the images at hand using configuration files. The install process and survey-specific set-ups are described in the `Getting Started <setup.html>`_.

----------------

.. toctree::
   :maxdepth: 2

   setup
   pipeline
   implementation
   misc

Acknowledgements
================

.. Supervisors, sextractor, scamp, swarp, hygcat, astrowise, skybot, astropy, pandas,