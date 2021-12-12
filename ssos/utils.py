#!/usr/bin/env python
""" Utility functions for the ssos pipeline. """

import argparse
import logging
import os
import re
import sys
import time
import warnings

from astropy.coordinates import SkyCoord
from astropy.coordinates import IllegalSecondWarning
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.time import Time
import astropy.units as u
from astropy import wcs
import pandas as pd
from sbpy.data import Names
from tqdm import tqdm

warnings.simplefilter("ignore", category=VerifyWarning)


def init_argparse():
    """Executes the setup of the arguments parser for the command-line api

    :return: args - command line argument parser
    """
    parser = argparse.ArgumentParser(description=f"The ssos pipeline")

    parser.add_argument(
        "fields", nargs="+", help=f"Path to directory of field to be searched"
    )
    parser.add_argument(
        "-c",
        "--config",
        action="store",
        dest="set_up",
        help=f"Path to configuration file",
        nargs=1,
    )
    parser.add_argument(
        "-d", "--default", action="store_true", help=f"Copy default setup files to CWD"
    )
    parser.add_argument(
        "-i",
        "--inspect",
        action="store_true",
        help=f"Launch visual inspection of SSO candidates",
    )
    parser.add_argument(
        "-l",
        "--log",
        action="store",
        dest="log",
        nargs=1,
        help=f"Set the logging level. Valid arguments are "
        f"DEBUG, INFO, WARNING, ERROR, CRITICAl",
        default=["INFO"],
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help=f"Suppress logging to console"
    )
    parser.add_argument(
        "-t",
        "--target",
        action="store",
        dest="target",
        help=f"Target directory for analysis results",
        nargs=1,
    )
    parser.add_argument("--sex", action="store_true", help=f"Force SExtractor runs")
    parser.add_argument("--scamp", action="store_true", help=f"Force SCAMP runs")
    parser.add_argument("--swarp", action="store_true", help=f"Force SWARP runs")
    parser.add_argument("--skybot", action="store_true", help="Force SkyBoT query")

    # Override filter parameters using the CLI
    group = parser.add_argument_group("Pipeline Settings")

    with open(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "default.ssos")
    ) as default:
        for line in default:

            # Skip comments and newlines
            if not line[0].isalpha():
                continue

            param, val, *_ = line.split()

            # If boolean argument
            if val in ["True", "False"] and param not in ["WEIGHT_IMAGES"]:
                group.add_argument(
                    f"--{param}",
                    metavar="bool",
                    action="store",
                    help=f"Override {param} setting. Must be " f"True or False.",
                    choices=["True", "False"],
                )
            # If parameter argument
            else:
                group.add_argument(
                    f"--{param}",
                    metavar="value",
                    action="store",
                    help=f"Override {param} setting.",
                )

    group = parser.add_argument_group("SCAMP Settings")
    for param in ["ASTREF_CATALOG", "CROSSID_RADIUS"]:
        group.add_argument(
            f"--{param}",
            metavar="value",
            action="store",
            help=f"Override {param} setting.",
        )
    args = parser.parse_args()
    return args


def create_target_dir(args):
    """Creates the target directory file structure

    :args: command line argument parser

    return
    ------
    target_dir - str, absolute path to target directory
    paths - dict, dict containing the absolute paths of the subdirectories
    """

    if args.target:
        target_dir = os.path.abspath(args.target[0])
    else:
        target_dir = os.getcwd()

    paths = {
        "cats": os.path.join(target_dir, "cats"),
        "checkplots": os.path.join(target_dir, "checkplots"),
        "cutouts": os.path.join(target_dir, "cutouts"),
        "logs": os.path.join(target_dir, "logs"),
        "images": os.path.abspath(args.fields[0]),
        "weights": os.path.join(target_dir, "weights"),
        "tmp": os.path.join(target_dir, "tmp"),
    }

    for path in paths.values():
        os.makedirs(path, exist_ok=True)

    return target_dir, paths


def init_logger(args, log_dir):
    """Initiate logging. By default, the log is written to file.
    The filename consists of the string 'sso_' followed by
    the current date and time.

    input
    ------
    args - command line argument parser
    log_dir - str, absolute path to log directory

    return
    ------
    log - logger instance
    log_file - str, absolute path to log file
    """

    # See if logging level was set via argument. Else, default to INFO
    numeric_level = getattr(logging, args.log[0].upper(), None)

    if not isinstance(numeric_level, int):
        numeric_level = 20

    log = logging.getLogger(__name__)
    log.setLevel(numeric_level)

    start_time = time.struct_time(time.localtime())
    log_file = "ssos_{:s}.log".format(time.strftime("%Y%m%d%H%M%S", start_time))

    file_log = logging.FileHandler(os.path.join(log_dir, log_file))
    file_log.setLevel(numeric_level)
    file_log.terminator = ""  # suppress newline characters
    log.addHandler(file_log)

    # If the quiet flag is set, logs to console are supressed
    if not args.quiet:
        stdout_log = logging.StreamHandler()
        stdout_log.setLevel(numeric_level)
        stdout_log.terminator = ""  # suppress newline characters
        log.addHandler(stdout_log)

    def log_excepthook(execption_type, exception_value, traceback, logger=log):
        """Function to pipe tracebacks into the logfile"""
        log.error(
            "Exception occurred during analysis: \n",
            exc_info=(execption_type, exception_value, traceback),
        )
        log.error("\n")

    sys.excepthook = log_excepthook
    return log, log_file, start_time


def unpack_header_kw(hdus, keywords, try_first=0):
    """Looks up keyword values in FITS headers. Iterates over
    extensions until keyword is found, else returns False

    Can look up single keyword or several at once

    hdus - HDUList, exposure headers
    keywords - string or list, keywords to look up
    try_first - header index to look up first

    returns
    ----
    vals - single parameter value or list of parameter values
    """

    vals = []

    kws = [keywords] if not isinstance(keywords, list) else keywords

    for kw in kws:
        for hdu_ind in set([try_first, *range(len(hdus))]):
            try:
                hdu = hdus[hdu_ind]
                vals.append(hdu.header[kw])
                break
            except (IndexError, KeyError):
                pass
        else:
            vals.append(False)

    if not isinstance(keywords, list):
        return vals[0]

    return vals


def check_wcs(images, sci_ext, log):
    """Checks input images for WCS header information
    and computes mean RA DEC of the observation field

    :images: list - absolute paths to images
    :sci_ext: int or str - index of SCI extension in images
    :log: logger instance - pipeline log
    :returns: bool -  presence of WCS coordinates,
              float, float - random pair of RA DEC coordinates
              str - OBJECT keyword value
    """

    for i, image in enumerate(images):

        kws = [
            "CTYPE1",
            "CTYPE2",
            "CRVAL1",
            "CRVAL2",
            "CD1_1",
            "CD1_2",
            "CD2_1",
            "CD2_2",
        ]
        # 'CUNIT1', 'CUNIT2']
        vals = {}

        with fits.open(image) as exp:

            for kw in kws:
                vals[kw] = unpack_header_kw(exp, kw, sci_ext)

                if not vals[kw]:  # keyword could not be found
                    log.info(
                        f"\n\nWCS header check failed for {image}, keyword"
                        f" {kw} not present.\n\n"
                    )
                    raise PipelineDependencyError(
                        f"The images do not contain valid WCS"
                        f" information in the header. Is the"
                        f" SCI_EXTENSION correct?"
                    )

            if i == 0:
                # Check for OBJECT keyword
                object_ = unpack_header_kw(exp, "OBJECT", sci_ext)

                if object_ is False:
                    object_ = ""

    return True, vals["CRVAL1"], vals["CRVAL2"], object_


def compute_image_center(header):
    """Uses WCS information to compute RA DEC of image center

    :header: astropy HDU.header - image header
    :returns: float, float - RA DEC of image center

    """
    # Projection keywords that throw off scamp
    patterns = [
        r"PV[0-9]+_[0-9]+",
        r"TR[0-9]+_[0-9]+",
        r"WAT[0-9]+_[0-9]+",
        r"LTM[0-9]+_[0-9]+",
    ]

    # remove bad keywords
    rm = []

    for kw in header:
        for pattern in patterns:
            if bool(re.match(pattern, kw)):
                rm.append(kw)
                break
    for kw in rm:
        del header[kw]

    warnings.simplefilter("ignore", category=wcs.FITSFixedWarning)
    warnings.simplefilter("ignore", category=IllegalSecondWarning)

    # Read in WCS
    w = wcs.WCS(header, fix=False)

    # Read in pixel limits
    naxis1 = header["naxis1"]
    naxis2 = header["naxis2"]

    ra, dec = w.wcs_pix2world(naxis1 / 2, naxis2 / 2, 0)
    return ra, dec


def create_clean_image(image, tmp_path, update=False):
    """Creates a copy of the input image with cleaned
    headers for SCAMP

    :image: str - absolute path to input image
    :tmp_path: str - aboslute path to temporary copy
    :update: bool - whether to change the image inplace
    """

    # Projection keywords that throw off scamp
    patterns = [
        r"PV[0-9]+_[0-9]+",
        r"TR[0-9]+_[0-9]+",
        r"WAT[0-9]+_[0-9]+",
        r"LTM[0-9]+_[0-9]+",
    ]
    # r'A_[0-9]_[0-9]', r'AP_[0-9]_[0-9]',
    # r'B_[0-9]_[0-9]', r'BP_[0-9]_[0-9]',
    # r'AP_ORDER', r'BP_ORDER']

    with fits.open(
        image, "update" if update else "readonly", output_verify=False
    ) as exp:
        for hdu in exp:
            # Iterate over headers and remove keywords
            # which mess up SCAMP
            rm = []

            for kw in hdu.header:
                for pattern in patterns:
                    if bool(re.match(pattern, kw)):
                        rm.append(kw)
                        break
            for kw in rm:
                del hdu.header[kw]

            # Change CTYPEX from ZPN to TAN
            # if 'CTYPE1' in exp[0].header:
            # exp[0].header['CTYPE1'] = 'RA---TAN'
            # exp[0].header['CTYPE2'] = 'DEC--TAN'
            if "EQUINOX" in exp[0].header:
                exp[0].header["EQUINOX"] = float(exp[0].header["EQUINOX"])

        if not update:
            exp.writeto(tmp_path, output_verify="ignore", overwrite=True)
        else:
            exp.flush()


def convert_to_mpc(infile):
    """Convert ssos CSV output to MPC format. Saves result to file.

    :infile: str - absolute path to ssos CSV output
    """
    # Outfile in same directory as infile
    outfile = infile.replace(".csv", "_mpc.txt")

    # Read in observations and sort
    data = pd.read_csv(infile, low_memory=False)
    data = data.sort_values(["SOURCE_NUMBER", "EPOCH"])

    # ------
    # Get some observation metadata
    add_details = input(f"Add template observation details header? [Y|n]:")
    if add_details in ["n", "N"]:
        add_details = False
    else:
        add_details = True

    MAG_BLANK = False
    if "MAG_CALIB" not in data.columns:
        mag_dec = input(
            f"No calibrated magnitude found in MAG_CALIB column. "
            f"Leave blank? [Y|n]:"
        )

        if mag_dec in ["n", "N"]:
            print("Add calibrated magnitudes as MAG_CALIB column.")
            sys.exit()
        else:
            MAG_BLANK = True

    if "OBS_CODE" not in data.columns:
        obs_code = input("No OBS_CODE column found. Enter observatory code: ")
    else:
        obs_code = data["OBS_CODE"].values[0]

    if len(obs_code) != 3:
        print("Observatory code has to consist of three characters.")
        sys.exit()

    # ------
    # Write output file line by line
    with open(outfile, "w") as output:

        if add_details:
            header = (
                f"COD {obs_code}\n"
                f"CON M. Name, Institute\n"
                f"CON [myemail@institute.fr]\n"
                f"MEA M. Mame, M. Collaborator\n"
                f"TEL 0.83-m reflector + CCD\n"
                f"COM Explanatory comment\n"
                f"NET Gaia-DR2\n"
                f"ACK Acknowledgment line\n"
                f"AC2 mycollaborator@institute.fr\n"
                f"\n"
            )

            output.write(header)

        # ------
        # Create text file in submission format
        desi = ""  # to keep track of discoveries
        for index, row in tqdm(
            data.iterrows(),
            total=data.shape[0],
            unit="obs",
            desc="Converting to MPC format",
        ):

            # Get number, MPC designation, discovery state
            number, desi, disco = _number_designation(row, desi)

            # Note 1, which handles errors, empty note
            note1 = " "
            # Append note 2, C for CCD
            note2 = "C"

            # Add date of observation
            date_obs = _date_of_observation(row)

            # Add coordinates of detection
            ra, dec = _coordinates(row)

            # Add magnitude and band
            if not MAG_BLANK:
                mag, band = _magnitude(row)
            else:
                mag, band = "", ""

            # Bring it all together
            # Columns     Format   Use
            #  1 -  5       A5     Minor planet number
            #  6 - 12       A7     Provisional or temporary designation
            # 13            A1     Discovery asterisk
            # 14            A1     Note 1
            # 15            A1     Note 2
            # 16 - 32              Date of observation
            # 33 - 44              Observed RA (J2000.0)
            # 45 - 56              Observed Decl. (J2000.0)
            # 57 - 65       9X     Must be blank
            # 66 - 71    F5.2,A1   Observed magnitude and band
            #                         (or nuclear/total flag for comets)
            # 72 - 77       X      Must be blank
            # 78 - 80       A3     Observatory code
            line = (
                f"{number:0>5}{desi:<7}{disco:<1}{note1}{note2}"
                f'{date_obs}{ra} {dec}{"":<9}{mag:<5}{band:<1}'
                f'{"":<7}{obs_code}\n'
            )
            if len(line) != 81:  # 80 + newline character
                print(
                    f"Error: This line has {len(line) - 1} instead of "
                    f"80 characters:"
                )
                print(line)
                sys.exit()
            output.write(line)


def _number_designation(row, prev_desi):

    # If object is unknown
    if not row.MATCHED:
        number = " " * 5
        desi = row["SOURCE_NUMBER"]
        disco = "*" if desi != prev_desi else " "
        return number, desi, disco

    # ------
    # Matched in SkyBoT
    # See if object is numbered
    if Names.asteroid_or_comet(row["SKYBOT_NAME"]) == "asteroid":

        if not pd.isna(row["SKYBOT_NUMBER"]):

            packed_number = Names.to_packed(str(int(row["SKYBOT_NUMBER"])))
            number = packed_number.rjust(5, "0")

            # If the asteroid is numbered, do not add name
            desi = ""
            disco = " "
            return number, desi, disco

        # Known but not numbered
        desi = Names.to_packed(row["SKYBOT_NAME"])
        return "", desi, ""

    else:  # it's a comet
        # Comet
        # Columns     Format   Use
        # 1 -  4       I4     Periodic comet number
        # 5            A1     Letter indicating type of orbit
        number = row["SKYBOT_NAME"][:-1].rjust(4, "0")  # Periodic comet number
        number += row["SKYBOT_NAME"][-1]  # Letter indicating type of orbit
        return number, "", ""


def _date_of_observation(row):
    # Append date of observation. Converting from MJD to YYYY MM DD.dddddd
    date_mjd = Time(row["MID_EXPOSURE_MJD"], format="mjd")

    # YYYY MM DD part
    date = date_mjd.iso.split(" ")[0].replace("-", " ")

    # Decimal day, taken from MJD
    date += "." + str(row["MID_EXPOSURE_MJD"]).split(".")[1][:6].ljust(6, "0")
    return date


def _coordinates(row):

    # Now append RA and DEC like ' 14 28 18.11 +21 34 01.7'
    ra = row["RA"]
    dec = row["DEC"]

    coords = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame="icrs")

    ra_seconds = "%.2f" % coords.ra.hms.s
    ra_seconds = ra_seconds.zfill(5)
    dec_seconds = "%.1f" % coords.dec.dms.s
    dec_seconds = dec_seconds.strip("-").zfill(4)

    coords = coords.to_string("hmsdms")
    ra, dec = coords.split()

    ra = ra.replace("h", " ").replace("m", " ").replace("s", "")
    dec = dec.replace("d", " ").replace("m", " ").replace("s", "")

    ra = " ".join(ra.split()[:-1] + [ra_seconds])  # cut to decimal places
    dec = " ".join(dec.split()[:-1] + [dec_seconds])  # cut to decimal places
    return ra, dec


def _magnitude(row):
    # Append observed magnitude and band
    mag = str(round(row["MAG_CALIB"], 1))  # only giving 0.1 accuracy
    band = row["FILTER"][0]

    return mag, band
