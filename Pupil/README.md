# README

## Preprocessing

The blink-interpolated pupil data (csv files) are in the `DataBlinksInterpolated_csv` directory (download all data [here](https://osf.io/3buew/)). 

To fit the deconvolution model with the [pypillometry package](https://github.com/ihrke/pypillometry), use `deconvolve_pupil.py` (implemented with Python 3.8.5), which uses parallel processing to fit the model to each participant, and each participant's results/data will be saved in the `DataPreprocessed` directory (an example `6.pd` is provided). 

After fitting the model to all participants (each participant should have a `.pd` file), use `get_prf_coef.py` to get the pupil response coefficient for each event in the task (this step also requires the `*event.csv` files in the `DataBlinksInterpolated_csv` directory). The results of this step will be saved in the `Data_trial_PRF_coef_baseline` directory (1 csv file participant). 

Prepare the data (csv files in `Data_trial_PRF_coef_baseline`) for statistical modeling with `prep_pupil.R`. The results of this step is one csv file in `Data_analysis` (i.e., `pupil_coefs.csv`).

## Analysis

Run `pupil_coef_stats_figs.R` to replicate the results in the paper. This R script relies on `Data_analysis/pupil_coefs.csv`.
