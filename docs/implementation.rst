###########################
Code Implementation
###########################


SExtractor and SCAMP wrapper
============================


Filter
======

Target Directory
================

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
