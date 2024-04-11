% This is for my own reference: my earlier unorganized code can be found at:
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/loop_source_data.m
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/Import_mne_headmodel.ipynb
% Created by Zhibin 4/11/2024
% #####################################################################################################################


% This script process the EEG data from 32 channel TMSi Device

% This script does the source localization (inverse solution) for the 
% cleaned EEG data from 32 channel TMSi Device
% And run PCA to aggreate source data in each ROI of the 448 cortical areas
clear
%% Construct brain model and create forward matrix in MNE


%% Source localization 


%% PCA to aggregate cortical source data