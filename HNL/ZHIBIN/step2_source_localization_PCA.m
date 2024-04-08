% This is for my own reference: my earlier unorganized code can be found at:
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/organize_62stroke/step1_clean_data.m
% Created by Zhibin 4/6/2024
% #####################################################################################################################

% This script does the source localization (inverse solution) for the 256 channel EGI raw EEG data
% And run PCA to aggreate source data in each ROI of the 448 cortical areas
clear
%% Construct brain model and create forward matrix in MNE
% Refer to MNE_construct_source.ipynb to create the leadfield matrix
% load forward matrix and source
load('base_files/MNE/EGI256/MNE_source_model.mat','leadfield','source_rr');
% leadfield: leadfield matrix from MNE (256 EEG channels X 5124 brain sources)
% source_rr: coordinates of the 5124 brain sources

%% Source localization 


% roiNames_250: labels for 463 ROIs (448 cortical + 15 subcortical ROIs)
% scale250_subcortROIs: indices for 15 subcortical ROIs in roiNames_250

% corti_ave_source_coor: coordinates for 448 cortical ROIs
% corti_ave_source_labl: indices for 448 cortical ROIs in roiNames_250
% corti_roiNames: names of the 448 cortical ROIs

% source_labels: indices for the 5124 brain sources in roiNames_250

% Brain: no need, delete from mat?
% parcels: contain labels in the brain, no need delete from mat?
clear brain parcels

% get inverse matrix for source localization
% You will to add the following repo to your path:
% https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh/tree/main/Simulations/util
addpath ../../../AdaptiveGraphicalLassoforParCoh/Simulations/util % this is my local path
[inversemat] = inversemodel(leadfield,'prctile',1);


subj_files=[0:1:60];
f=3; % select one patient to run source localization (f=1:61)

% Navigate to your cleaned preprocessed EEG data directory:
% cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_cleaned
% cd /ssd/zhibin/archive/EEG_stroke_62_cleaned
% load the preprocessed EEG file for one patient
load([num2str(subj_files(f)) '.mat'],'preprocessed_eeg','subject_ID', ...
    'chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');
display(['start processing subject file: ' num2str(subj_files(f)) '.mat']);

% Using inverse solution to covert preprocessed eeg to localized source data
source_data=inversemat*preprocessed_eeg;

%% PCA to aggregate cortical source data
% PCA to aggregate source data in each ROI
% fra_eigenvalues=zeros(1,max(unique(source_labels))); % no need to include
corti_source_data=[];
for sr=1:max(unique(source_labels))
    I=find(source_labels==sr); % need to load save and load this source_labels
    if ~isempty(I)
        [~, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
        % fra_eigenvalues(sr)=LATENT(1)/sum(LATENT); % weight of 1st eigenvalue
        corti_source_data=[corti_source_data SCORE(:,1)];
    end
end % 66s

% remove the non-zeros fraction of eigenvalues and subcortical rois
% corti_fra_eigenvalues=fra_eigenvalues(ave_source_label); % no need to include
% need to load save and load this ave_source_label
    
% remove subcortical roi
corti_source_data(:,ind)=[];
% need to load save and load this ind?

% save source data of the aggregated source data for cortical ROIs
% navegate the your destination directory, such as
% cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
save([num2str(subject_ID) '.mat'],'corti_fra_eigenvalues', ...
        'corti_source_data','corti_ave_source_coor','corti_ave_source_labl', ...
        'subject_ID','chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');

display(['Complete one file: ' num2str(subj_files(f)) '.mat ********************'])


