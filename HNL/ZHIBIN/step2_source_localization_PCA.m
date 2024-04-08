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
% load the brain model, ROI names and localation files for the 463 ROIs
load('base_files/brain_ROIs.mat','Brain','roiNames_250','scale250_subcortROIs','parcels')
% Brain: The brain model of Lausanne2008_fsaverageDSsurf
% parcels: contain indices of 463 ROIs in a 3d brain matrix (256x256x256)
% roiNames_250: labels for 463 ROIs (448 cortical + 15 subcortical ROIs) (scale 250)
% scale250_subcortROIs: indices for 15 subcortical ROIs in roiNames_250

% label the sources
Vertex=Brain.Vertex;
% Allign the 5124 sources with 463 ROIs and label them
x_shift=(max(Vertex(:,1))-max(source_rr(:,1))*1e3)/2+(min(Vertex(:,1))-min(source_rr(:,1))*1e3)/2;
y_shift=(max(Vertex(:,2))-max(source_rr(:,2))*1e3)/2+(min(Vertex(:,2))-min(source_rr(:,2))*1e3)/2;
z_shift=(max(Vertex(:,3))-max(source_rr(:,3))*1e3)/2+(min(Vertex(:,3))-min(source_rr(:,3))*1e3)/2;
source_x=source_rr(:,1) * 1e3 + x_shift;
source_y=source_rr(:,2) * 1e3 + y_shift;
source_z=source_rr(:,3) * 1e3 + z_shift;
source_xyz=[source_x source_y source_z];
num_source=size(source_xyz,1);
source_fsaverage = source_xyz+127.5; % 127.5 is based on the fsaverage volume being 256 x 256 x 256
source_labels=zeros(num_source,1);
for i = 1:length(source_fsaverage)
    vox = floor(source_fsaverage(i,:)); % change from ceil to floor,now we have 2 subcortical not mapped
    inds              = sub2ind([size(parcels)], vox(1), vox(2), vox(3));
    label             = parcels(inds); 
    source_labels(i) = label;
end
% find which ROIs not mapped
roiNames_250(setdiff(1:463,unique(source_labels))) % 227 457 459
setdiff(1:463,unique(source_labels))

% get coordinate and label for each ROIs to be aggregated
ave_source_coor=[];
ave_source_label=[];
for sr=1:max(unique(source_labels))
    I=find(source_labels==sr);
    if ~isempty(I)
        ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
        ave_source_label=[ave_source_label; sr];
    end
end 

% create a boolean of subcortical rois
bool_subcorti=zeros(1,length(ave_source_label));
for i=1:length(ave_source_label)
    clear tmp
    tmp=ave_source_label(i);
    if ismember(tmp, scale250_subcortROIs)
        bool_subcorti(i)=1;
    end
end
sum(bool_subcorti)
ind=find(bool_subcorti);
ind % use these subcortical indices to remove subcortical aggregated pca data

% remove subcortical ROIs
ave_source_coor(ind,:)=[];
ave_source_label(ind)=[];
% save the cortical ROI coordinates, indicies and labels 
corti_ave_source_coor=ave_source_coor;
corti_ave_source_labl=ave_source_label;
corti_roiNames=roiNames_250(corti_ave_source_labl);
save('./base_files/corti_ave_source.mat','corti_ave_source_coor','corti_ave_source_labl','corti_roiNames');
% corti_ave_source_coor: coordinates for 448 cortical ROIs
% corti_ave_source_labl: indices for 448 cortical ROIs in roiNames_250
% corti_roiNames: names of the 448 cortical ROIs

% get inverse matrix for source localization
% You will to add the following repo to your path:
% https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh/tree/main/Simulations/util
addpath ../../../AdaptiveGraphicalLassoforParCoh/Simulations/util % this is my local path
[inversemat] = inversemodel(leadfield,'prctile',1);

% select one patient to run source localization (f=1:61)
subj_files=[0:1:60]; % all patient files
f=3; % select one patient 

% load the cleaned EEG file
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


