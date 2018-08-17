###############
Getting started
###############

Installing the Pipeline
=======================

The pipeline requires the following additional `python` packages: `astropy`, `numpy`, `pandas`, `scipy`, and `statsmodels`. You can quickly install them using the following command

``$ pip install -r requirements.txt``

Make sure that the astrOmatic binaries `sex`, `scamp`, and `swarp` are in your `PATH` shell variable. Try it with e.g.

`
$ sex --version

SExtractor version 2.19.5 (2015-06-19)

$ scamp --version

SCAMP version 2.0.4 (2015-06-19)

$ swarp --version

SWarp version 2.38.1 (2018-06-28)
`

In the ASCII text file `pipeline_settings.sso`, replace the `path/to/jvar_semp` with the absolute path to the folder. This has to be done for 7 settings in total.


Survey-specific changes
=======================

Pipeline Setting Files
======================