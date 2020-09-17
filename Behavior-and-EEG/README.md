# README

## Preprocessing

The [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) MATLAB toolbox is required. Other custom MATLAB functions might be called too—get in touch if you want the continuous data or/and plan to segment (i.e., create epochs) the continuous data with `Epoching_FXCEEG.m`.

### Epoching continuous EEG data

To create epochs from the continuous EEG data that are locked to cue, stimulus, response, and feedback, use `Epoching_FXCEEG.m`. The resulting `.mat` files will be saved in a directory `Export`.

### Export single-trial EEG time-windows

To export single-trial time-windows from the epoched data from the step above, use `export_FXC_EEG.m`. The resulting `.mat` files (e.g., `P3b250550.mat`, `CNV10001500.mat`, `Baseline2000.mat`) will be saved in the same `Export` directory (see previous step). 

### Analysis

To run the analyses reported in the manuscript, use `Froemer_Lin_DeanWolf_Inzlicht_Shenhav_2020.R`. The behavioral data can be found  in the `Behavior` directory:

- Study 1: `FXC_102Data_010920.csv`
- Study 2: `FXCallSubDataTable.csv`

The EEG (`.mat`) files will be in the `Export` directory.

To plot the P3b and CNV ERP components, use `plotting_FXC_EEG.m`, which reads files from the `Export` directory and requires custom MATLAB functions ([export_fig](https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig)).
