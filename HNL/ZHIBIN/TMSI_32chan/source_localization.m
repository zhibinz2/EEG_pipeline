% This is for my own reference: my earlier unorganized code can be found at:
% open ../../../../MEG_EEG_Source_Localization/PCA_32chan_AGL/loop_source_data.m
% open ../../../../MEG_EEG_Source_Localization/PCA_32chan_AGL/Import_mne_headmodel.ipynb

% Created by Zhibin 4/11/2024
% #####################################################################################################################


% This script process the EEG data from 32 channel TMSi Device

% This script does the source localization (inverse solution) for the 
% cleaned EEG data from 32 channel TMSi Device
% And run PCA to aggreate source data in each ROI of the 448 cortical areas
clear
%% Construct brain model and create forward matrix in MNE


%% Source localization 

cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
load('source_rr.mat');
load('Lausanne2008_fsaverageDSsurf_60_125_250.mat')

% label the sources
Vertex=Brain.Vertex;
load('parcels.mat') % This is the labels
% Anni's labeling method
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

% find which label not mapped
roiNames_250(setdiff(1:463,unique(source_labels))) % 225 227

% load forward matrix
load('leadfield.mat');

% inverse model
addpath ../../AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat, stat, reconstructed] = inversemodel(leadfield,'prctile',1);

data_path= '../../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 

cd ../../Cleaned_data

%% PCA to aggregate cortical source data

cd 
data=load([data_path  'clean_' runid '.mat'],'dataL','dataR');
% select a subject's cleaned EEG data
preprocessed_eeg=data.dataL;
% select a trial
tr=1;
EEG_ori=preprocessed_eeg{tr}(:,1:32)';
source_data=inversemat*EEG_ori;

agr_source_data=[];
ave_source_coor=[];
ave_source_label=[];
for sr=1:max(unique(source_labels))
    I=find(source_labels==sr);
    if ~isempty(I)
        [COEFF, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
        agr_source_data=[agr_source_data SCORE(:,1)];
        ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
        ave_source_label=[ave_source_label; sr];
    end
end

% save the source data
save([data_path 'source_data/' num2str(runid) '/subj' num2str(subj) '_tr_' num2str(tr)  '.mat'], ...
                            'fra_eigenvalues','agr_source_data','ave_source_coor','ave_source_label');
