#%%

import concurrent.futures
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypillometry as pp

print(f"\nSystem version\n{sys.version}")

#%% set up paths

PATH = {}
PATH["data_dir"] = "DataBlinksInterpolated_csv"
PATH["out_dir"] = "DataPreprocessed"
# PATH["outpdf_dir"] = "DataPreprocessed_pdf"

#%% function for deconvolution


def deconvolve_pupil(subject, subset=None):
    print(f"Reading data")

    dt = pd.read_csv(os.path.join(PATH["data_dir"], f"{subject}data.csv"))
    event = pd.read_csv(os.path.join(PATH["data_dir"], f"{subject}event.csv"))
    # add time column to event by merging with dt
    event = event.merge(dt[["sample", "time"]], "left", ["sample"])
    # select events
    event = event.query(
        "type in ['TRIALID', '91', '92', '93', '94', '51', '52', '53', '200', '210', '10', '20', '110']"
    )

    print(f"Creating PupilData object")
    d = pp.PupilData(
        dt["pupil"],
        time=dt["time"],
        name=f"subject_{subject}",
        event_onsets=event["time"],
        event_labels=event["type"],
    )

    # subset data (for testing)
    if subset:
        print(f"Subsetting")
        d = d.sub_slice(subset[0], subset[1], units="min")

    # downsample
    print(f"Downsampling")
    d = d.downsample(20)

    # estimate baseline
    # print(f"Estimating baseline")
    # d = d.estimate_baseline(method="envelope_iter_bspline_1", fsd=5, verbose=0)
    # d.sy = d.sy - d.baseline

    # z-score
    print(f"z-score")
    # z-score
    d = d.scale()

    # estimate response
    print(f"Estimating pupil response")
    d = d.estimate_response(npar=10.1, tmax=1300)

    # save data
    fname = os.path.join(PATH["out_dir"], f"{subject}.pd")
    print(f"Saving data {fname}")
    d.write_file(fname)

    # save pdf
    # d.plot_segments(pdffile=os.path.join(PATH["outpdf_dir"], f"{subject}.pdf"))

    return d


#%% test

# deconvolve_pupil(4, subset)

#%% subjects

subjects_list = [
    2,
    3,
    4,
    8,
    9,
    11,
    13,
    17,
    18,
    20,
    23,
    25,
    26,
    28,
    29,
    31,
    34,
    37,
    40,
    41,
    44,
    49,
    50,
    54,
    56,
    57,
    58,
    61,
    62,
    65,
    66,
    67,
    68,
    70,
    71,
    72,
    73,
]

#%% fit model with parallel loop

cores = 40
timeStart = time.time()
with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
    result = executor.map(deconvolve_pupil, subjects_list)
timeEnd = time.time()

print(f"\nTime to completion: {timeEnd - timeStart}")
