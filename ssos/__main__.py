#!/usr/bin/env python
"""
    Pipeline to search for Solar System objects in wide-field imaging surveys
    Information on the project can be found in arXiv:1711.02780
    and arXiv:1906.03673

    For questions, contact max.mahlke (at) oca.eu

    Call as:    ssos
"""

import os
import shutil
import sys
import time


if len(sys.argv) == 1:
    from ssos import GREETING

    print(GREETING)
    sys.exit()

# Copy default configuratin files to CWD
if sys.argv[1] in ["-d", "--default"]:

    path_to_module = os.path.dirname(__file__)

    try:
        shutil.copytree(
            os.path.join(path_to_module, "semp"), os.path.join(os.getcwd(), "semp")
        )
    except FileExistsError:
        print(
            f'Path {os.path.join(os.getcwd(), "semp")} exists in CWD,'
            f" not overwriting. Exiting."
        )
        sys.exit()

    if not os.path.isfile(os.path.join(os.getcwd(), "default.ssos")):
        shutil.copy(
            os.path.join(path_to_module, "default.ssos"),
            os.path.join(os.getcwd(), "default.ssos"),
        )
    else:
        print(f"File default.ssos exists in CWD, not overwriting. Exiting.")
    sys.exit()

elif sys.argv[1] in ["-i", "--inspect"]:
    from ssos.ins import inspectCutouts

    inspectCutouts(sys.argv[2:])
    sys.exit()

elif sys.argv[1] in ["-m", "--mpc"]:
    from ssos.utils import convert_to_mpc

    convert_to_mpc(sys.argv[2])
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
    pipeline.run_SExtractor()

    # ------
    # Run SCAMP on SExtractor catalogues
    if pipeline.settings["REMOVE_REF_SOURCES"]:
        pipeline.run_SCAMP(
            crossid_radius=1,
            full_name="full_stars.cat",
            merged_name="merged_stars.cat",
            keep_refcat=True,
            adjust_SExtractor_and_aheader=True,
        )
        pipeline.run_SCAMP(
            crossid_radius=pipeline.args.CROSSID_RADIUS,
            solve_astronomy=False,
            pattern_matching=False,
        )

    else:
        pipeline.run_SCAMP(crossid_radius=pipeline.args.CROSSID_RADIUS)

    # ========
    # Catalogue creation done
    # ========
    pipeline.log.info(f'{pipeline.term_size * "-"}\n')

    # Call pipeline filter steps
    init_source_numb = pipeline.number_of_sources()
    pipeline.log.info(
        f'{"Initial sample size".ljust(22)}'
        f" {str(init_source_numb).ljust(len(str(init_source_numb)))} | "
        f'{"#" * (pipeline.term_size - 26 - len(str(init_source_numb)))}'
    )
    for step in pipeline.steps:
        pipeline.execute_filter(step)
        pipeline.log.info(
            f"{step.ljust(22)} {str(pipeline.number_of_sources()).ljust(len(str(init_source_numb)))} | "
            f'{"#" * int(((pipeline.term_size - 26 - len(str(init_source_numb)))) * pipeline.number_of_sources()/init_source_numb)}\n'
        )

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
    print(f'{"-" * pipeline.term_size}')

    for step in pipeline.analysis_steps:
        pipeline.execute_analysis(step)

    # ------
    # Rename columns, remove unnecessary data, save to file
    pipeline.save_and_cleanup()
    print(f'{"-" * pipeline.term_size}')

    pipeline.log.info(
        "\t|\t".join(
            [
                "\nAll done!",
                "%i SSO candidates found" % pipeline.number_of_sources(),
                "The analysis ran in %i seconds\n\n" % pipeline.run_time,
            ]
        )
    )

    pipeline.log.info(
        "Output File: %s\nLog File: %s\n\n"
        % (
            os.path.join(
                pipeline.paths["cats"],
                "ssos_{:s}.csv".format(
                    time.strftime("%Y%m%d%H%M%S", pipeline.start_time)
                ),
            ),
            os.path.join(pipeline.paths["logs"], pipeline.log_file),
        )
    )


if __name__ == "__main__":
    main()
