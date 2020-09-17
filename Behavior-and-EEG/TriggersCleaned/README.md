# Matching behavioral data and EEG triggers

I should have fixed the triggers so that they match the behavioral data, which means you can also start analysing response-locked ERPs!

## _triggersFixed.csv

subject: subject ID

eventN: the trigger number recorded by EEG system (starts from 1)

eventCode: trigger code indicating cue/target/response/feedback onset

timeS: time (seconds) at which trigger was sent

samplingPoint: sampling point at which trigger was sent

trial: trial number (should match the trial number in behavioral data); certain subjects (e.g., 2, 4, 41) will not start counting from 1 because research assistants forgot to start EEG recording in time

trialExcluded: trials to be excluded in EEG analysis (based on behavioral data) (1: exclude trigger/trial; 2: include trigger/trial)

eventCodeNew: new version of eventCode that excludes trials/triggers based on behavioral data (-999 are triggers that should be rejected)

## behavioralData_trigger_match.csv

subject: subject ID

dataCuesN: number of valid trials in behavioral data

triggerCuesN: number of valid cue triggers in EEG

triggerTargetsN: number of valid target (Stroop) triggers in EEG

triggerResponseN: number of valid response triggers in EEG

triggerFeedbackN: number of valid response triggers in EEG

**If everything is corrrect, all the numbers above should match (i.e., should be the same for each participant)**





