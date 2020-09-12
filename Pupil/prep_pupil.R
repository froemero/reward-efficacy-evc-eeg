# code to prepare data for statistical modeling
rm(list = ls())
library(data.table); library(hausekeep)

# read data
datadir_prf <- "Data_trial_PRF_coef_baseline"  # data directory with pupil response coefficients
dt <- lapply(list.files(datadir_prf, full.names = TRUE), fread)
dt <- rbindlist(dt)
trials2use <- fread("./Data_analysis/trials2use.csv")  # trials to use (excludes trials with RTs < 200 and no responses; to match other analyses)

# select trials
dt <- dt[trials2use, on = c("subject", "trial")]
dt <- dt[!is.na(pupil_response)]

# select columns
dt <- select(dt, subject, npar, tmax, event:event_onset, pupil_response)

# recode
dt[ttl %in% c(91, 92), Efficacy := -0.5]
dt[ttl %in% c(93, 94), Efficacy := 0.5]
dt[ttl %in% c(91, 93), Reward := -0.5]
dt[ttl %in% c(92, 94), Reward := 0.5]

dt[ttl %in% c(10), isRewarded := -0.5]
dt[ttl %in% c(20, 110), isRewarded := 0.5]

dt[ttl %in% c(52), `:=` (congruencyC = -0.5, congruency = 'congruent')]
dt[ttl %in% c(53), `:=` (congruencyC = 0.0, congruency = 'neutral')]
dt[ttl %in% c(51), `:=` (congruencyC = 0.5, congruency = 'incongruent')]

# zscore pupil
dt[, pupil_responseZ := zscore(pupil_response), by = subject]
boxplot(dt$pupil_responseZ)

# remove outliers
dt[, pupil_responseZ_clean := outliersMAD(pupil_responseZ, 3, digits = 7), by = subject]
boxplot(dt$pupil_responseZ_clean)

# zscore trial
dt[, Trial := zscore(trial), by = .(subject, event)]

# save data
fwrite(dt, "./Data_analysis/pupil_coefs.csv")