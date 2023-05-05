import os

from astropy.io import fits
from astropy.visualization import simple_norm, ZScaleInterval
from astropy.wcs import WCS
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import warnings

warnings.filterwarnings("ignore")


def press(event, ssos, index, source_number):

    ssos.loc[index, "INSPECTED"] = True

    if event.key in ["a", "right"]:
        ssos.loc[index, "ASTEROID"] = True
        print("%i: Asteroid" % source_number)

    if event.key in ["u", "up"]:
        ssos.loc[index, "UNKNOWN"] = True
        print("%i: UFO" % source_number)

    if event.key in ["d", "left"]:
        ssos.loc[index, "ASTEROID"] = False
        print("%i: Discarded" % source_number)

    plt.close()


def inspectCutouts(ana_dirs):

    # cheap hack for CMD API
    if "--force" in ana_dirs:
        force = True
    else:
        force = False

    for ana_dir in ana_dirs:

        if ana_dir == "--force":
            continue

        try:
            sso_path = os.path.join(
                sorted(
                    [
                        os.path.join(ana_dir, "cats/", csv)
                        for csv in os.listdir(os.path.join(ana_dir, "cats/"))
                        if csv.startswith("ssos") and csv.endswith("csv")
                    ]
                )[-1]
            )
        except IndexError:
            print("\nNo SSO candidates database found in %s" % ana_dir)
            continue

        # Open SSO candidates data
        ssos = pd.read_csv(sso_path)
        print(
            "Looking at %i candidates in %s" % (len(set(ssos.SOURCE_NUMBER)), ana_dir)
        )

        if "INSPECTED" not in ssos.columns:
            ssos["INSPECTED"] = False
            ssos["ASTEROID"] = False
            ssos["UNKNOWN"] = False
        elif not force:
            print(
                f"Already inspected {ana_dir}, skipping. Use --force to "
                f" override this behaviour."
            )
            continue

        for sn, sso in ssos.groupby("SOURCE_NUMBER"):

            # Get cutout filepaths
            cutouts = [
                os.path.join(
                    ana_dir, "cutouts", str(sn) + "_{:02d}".format(cat) + ".fits"
                )
                for cat in sso.CATALOG_NUMBER
            ]

            fig = plt.figure()
            fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
            fig.canvas.mpl_connect(
                "key_press_event", lambda event: press(event, ssos, sso.index, sn)
            )

            ims = []

            image = fits.open(cutouts[0])
            wcs = WCS(image[0].header, fobj=image)
            image = image[0]
            # ax = plt.subplot(projection=wcs)
            ax = plt.subplot()

            for i, cutout in enumerate(cutouts):

                try:
                    image = fits.open(cutout)[0]
                except OSError:
                    print(
                        "Could not find %s in cutouts directory. Was EXTRACT_CUTOUTS set to True?"
                        % cutout
                    )
                    continue

                data = image.data
                interval = ZScaleInterval()
                norm = simple_norm(
                    data,
                    "linear",
                    min_cut=min(interval.get_limits(data)),
                    max_cut=max(interval.get_limits(data)),
                )

                # When data is zero we get a large contrast,
                # set to background level instead
                try:
                    data[data == 0] = np.min(data[data != 0])
                except ValueError:  # complete frame zero
                    pass

                im = ax.imshow(
                    data,
                    cmap="gray",
                    norm=norm,
                    animated=True,
                )
                ax.text(
                    0.0,
                    1.05,
                    "Press Key - Left: Discard, Up: Unknown, Right: Asteroid",
                    verticalalignment="center",
                    transform=ax.transAxes,
                    zorder=100,
                )
                if type(sso.SKYBOT_NAME.values[0]) != float:
                    ax.text(
                        0,
                        1.1,
                        sso.SKYBOT_NAME.values[0],
                        verticalalignment="center",
                        transform=ax.transAxes,
                    )

                ims.append([im])

            ani = animation.ArtistAnimation(fig, ims, interval=150, blit=True)
            plt.show()

        print(
            "Found %i asteroids, %i unknown objects, and discarded %i candidates."
            % (
                len(set(ssos[ssos["ASTEROID"]]["SOURCE_NUMBER"])),
                len(set(ssos[ssos["UNKNOWN"]]["SOURCE_NUMBER"])),
                len(
                    set(ssos[(~ssos["ASTEROID"]) & (~ssos["UNKNOWN"])]["SOURCE_NUMBER"])
                ),
            )
        )
        ssos.to_csv(sso_path, index=False)
