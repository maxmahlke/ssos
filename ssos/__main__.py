#!/usr/bin/env python3
'''
  Pipeline to search for Solar System objects in wide-field imaging surveys
  Information on the project can be found in arXiv:1711.02780

  For questions, contact max.mahlke (at) cab.inta-csic.es

  Max Mahlke, August 2018
'''

import os
import sys
import time

if len(sys.argv) == 1:
    from ssos import GREETING
    print(GREETING)
    sys.exit()

if sys.argv[1] in ['-d', '--default']:
    path_to_module = os.path.dirname(__file__)
    os.system('cp -i -r {%s,%s} .' % (os.path.join(path_to_module, 'semp'),
                                   os.path.join(path_to_module, 'default.ssos')))
    sys.exit()

elif sys.argv[1] in ['-i', '--inspect']:
    from ssos.inspect import inspectCutouts
    inspectCutouts(sys.argv[2:])
    sys.exit()


from ssos.core import Pipeline


def main():
    # ========
    # Start the pipeline.  This step initializes the log, the target directory,
    # evaluates and checks the pipeline settings and verifies the input images.
    # ========
    pipeline = Pipeline()

    # ------
    # Run SExtractor on images
    pipeline.log.info('\nRunning SExtractor..\t')
    pipeline.run_SExtractor()

    # ------
    # Run SCAMP on SExtractor catalogues
    if pipeline.settings['REMOVE_REF_SOURCES']:
        pipeline.run_SCAMP(crossid_radius=1, full_name='full_stars.cat',
                           merged_name='merged_stars.cat', keep_refcat=True,
                           adjust_SExtractor_and_aheader=True)
        pipeline.run_SCAMP(crossid_radius=pipeline.args.CROSSID_RADIUS, solve_astronomy=False, pattern_matching=False)

    else:
        pipeline.run_SCAMP(crossid_radius=pipeline.args.CROSSID_RADIUS)

    # ========
    # Catalogue creation done
    # ========
    pipeline.log.info('\n --- Starting Filter pipeline ---\n\n')
    pipeline.log.info('%s %i\n' % ('All Sources'.ljust(20), pipeline.number_of_sources()))

    # Call pipeline filter steps
    for step in pipeline.steps:
        pipeline.execute_filter(step)
        pipeline.log.info('%s %i \n' % (step.ljust(20), pipeline.number_of_sources()))

    if not pipeline.added_proper_motion:
        pipeline.add_proper_motion()

    if not pipeline.added_SExtractor_data:
        pipeline.add_SExtractor_data()

    # ========
    # FILTERING COMPLETE
    # ========

    # ------
    # Optional analyses

    pipeline.add_image_metadata()

    for step in pipeline.analysis_steps:
        pipeline.execute_analysis(step)

    # ------
    # Rename columns, remove unnecessary data, save to file
    pipeline.save_and_cleanup()

    pipeline.log.info('\t|\t'.join(
                     ['\nAll done!',
                      '%i SSO candidates found' % pipeline.number_of_sources(),
                      'The analysis ran in %i seconds\n\n' % pipeline.run_time]))

    pipeline.log.info('Output File: %s\nLog File: %s\n\n' %
                     (os.path.join(pipeline.paths['cats'],
                                   'ssos_{:s}.csv'.format(time.strftime('%Y%m%d%H%M%S',
                                                          pipeline.start_time))),
                      os.path.join(pipeline.paths['logs'], pipeline.log_file)))


if __name__ == '__main__':
    main()
