% load source power
load('/home/zhibinz2/Documents/GitHub/EEG_pipeline/Celia/trouble_shooting20250225/Patient_power.mat')
sourcedelta=POWER.DELTA.Mean;
sourcebeta=POWER.BETA.Mean;

load('/home/zhibinz2/Documents/GitHub/EEG_pipeline/Celia/trouble_shooting20250225/scalp256/Patient_256_Power.mat')
eegdelta=POWER.DELTA.Mean;
eegbeta=POWER.BETA.Mean;