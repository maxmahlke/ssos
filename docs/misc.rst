Tips & Tricks
=============

.. role:: python(code)
   :language: python


* In the `pipeline_settings.sso` file, paths containing `$HOME` will be expanded to the user home directory following :python:`os.path.expanduser('~')`.