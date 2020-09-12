#%%

import concurrent.futures
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypillometry as pp
from scipy.stats import pearsonr

print(f"\nSystem version\n{sys.version}")

#%% set up paths

PATH = {}
PATH["data_dir"] = "DataPreprocessed"
PATH["event_dir"] = "DataBlinksInterpolated_csv"
PATH["out_dir"] = "Data_trial_PRF_coef_baseline"

#%% define function to export coefficients


def get_prf_coef(subject):
    try:
        print(f"Reading data {subject}")
        dd = pp.PupilData.from_file(os.path.join(PATH["data_dir"], f"{subject}.pd"))

        # events and onset times
        es = dd.event_labels
        ts = dd.event_onsets
        baseline = pp.stat_event_interval(dd.tx, dd.baseline, dd.event_onsets, [0, 0])
        pupil_response = dd.response_pars["coef"]

        # events
        event = pd.read_csv(os.path.join(PATH["event_dir"], f"{subject}event.csv"))
        event = event.query(
            "type in ['TRIALID', '91', '92', '93', '94', '51', '52', '53', '200', '210', '10', '20', '110']"
        )
        # split value column
        event_df = (
            event["value"].str.split(r"MSG\t\d+ ", expand=True).reset_index(drop=True)
        )
        event_df = event_df.rename(columns={0: "event", 1: "code"})
        event_df.loc[event_df["code"].str.contains("TRIAL"), "event"] = "fixation"
        event_df["ttl"] = np.float_(event_df["code"].str.replace("TRIALID ", ""))

        # convert ttl to condition
        new_trial_num = event_df.query("event == 'fixation'")["ttl"] / 1000
        event_df.loc[event_df["event"] == "fixation", "ttl"] = new_trial_num
        event_df.loc[event_df["event"] == "fixation", "trial"] = new_trial_num * 1000
        event_df.loc[event_df["ttl"].between(91, 94), "event"] = "cue"
        event_df.loc[event_df["ttl"].between(51, 53), "event"] = "stimulus"
        event_df.loc[event_df["ttl"].between(200, 210), "event"] = "response"
        event_df.loc[event_df["ttl"].isin([10, 20, 110]), "event"] = "feedback"
        event_df["trial"] = event_df["trial"].fillna(method="ffill")

        # add new columns
        event_df["event_onset"] = ts
        event_df["pupil_baseline"] = baseline
        event_df["pupil_response"] = pupil_response
        event_df["subject"] = subject

        event_df["npar"] = 10.1
        event_df["tmax"] = dd.response_pars["tmax"]
        event_df

        event_df["r"] = pearsonr(dd.sy, dd.response)[0]

        event_df.to_csv(os.path.join(PATH["out_dir"], f"{subject}.csv"), index=False)
    except Exception as e:
        print(f"error with subject {subject}")
        print(e)


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

timeStart = time.time()
with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
    result = executor.map(get_prf_coef, subjects_list)
timeEnd = time.time()

print(f"\nTime to completion: {timeEnd - timeStart}")
