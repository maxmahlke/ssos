Tips & Tricks
=============

.. role:: python(code)
   :language: python


* In the `pipeline_settings.sso` file, paths containing `$HOME` will be expanded to the user home directory following :python:`os.path.expanduser('~')`.

* A quick sanity check of the SkyBoT matching can be done using TOPCAT, by comparing the proper motion of the SSOs as retrieved from the pipeline and as predicted by SkyBoT. If the matched SSOs are the ones recovered by the pipeline (and not a coincidental match), you expect ratios close to one.

* `OBSERVATORY_CODE`, `FOV_DIMENSIONS`, `REFERENCE_FILTER` are the only settings which are not checked for correct input upon the pipeline initialization.