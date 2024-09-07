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
% Refer to ../TMSI_32chan/MNE_construct_source.ipynb to create the leadfield matrix
% load forward matrix and source
load('../base_files/MNE/TMSI32/leadfield.mat','leadfield');
load('../base_files/MNE/TMSI32/source_rr.mat','source_rr');
% leadfield: leadfield matrix from MNE (256 EEG channels X 5124 brain sources)
% source_rr: coordinates of the 5124 brain sources

%% Source localization 
load('../base_files/MNE/TMSI32/Lausanne2008_fsaverageDSsurf_60_125_250.mat','Brain','roiNames_250','scale250_subcortROIs')
% label the sources
Vertex=Brain.Vertex;

% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
load('../base_files/MNE/TMSI32/parcels.mat') % This is the labels

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
% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
load('../base_files/MNE/TMSI32/leadfield.mat');

% inverse model
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat, stat, reconstructed] = inversemodel(leadfield,'prctile',1);

data_path= '/home/zhibinz2/Documents/GitHub/Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 

cd /home/zhibinz2/Documents/GitHub/Cleaned_data

%% PCA to aggregate cortical source data
cd 
runid=num2str(20220713);
data=load([data_path  'clean_' runid '.mat'],'dataL','dataR');
% select a subject's cleaned EEG data
preprocessed_eeg=data.dataL;
% select a trial
tr=1;
EEG_ori=preprocessed_eeg{tr}(:,1:32)';
source_data=inversemat*EEG_ori;

fra_eigenvalues=zeros(1,max(unique(source_labels)));
agr_source_data=[];
ave_source_coor=[];
ave_source_label=[];
tic
for sr=1:max(unique(source_labels))
    I=find(source_labels==sr);
    if ~isempty(I)
        [COEFF, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
        fra_eigenvalues(sr)=LATENT(1)/sum(LATENT);
        agr_source_data=[agr_source_data SCORE(:,1)];
        ave_source_coor=[ave_source_coor; mean(source_fsaverage(I,:),1)];
        ave_source_label=[ave_source_label; sr];
    end
end
toc % 95s

% save the source data
cd /ssd/zhibin/1overf/Cleaned_sourcedata/source_data
% save([data_path 'source_data/' num2str(runid) '/subj' num2str(subj) '_tr_' num2str(tr)  '.mat'], ...
%                             'fra_eigenvalues','agr_source_data','ave_source_coor','ave_source_label');

%% Remove 16 subcortical and "zeros marked" sources not mapped and saved the aggreated source data
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
load('/home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/Lausanne2008_fsaverageDSsurf_60_125_250.mat')

marked_fra_eigenvalues=fra_eigenvalues(ave_source_label);
bool_temp=zeros(1,length(ave_source_label));
for i=1:length(ave_source_label)
    clear tmp
    tmp=ave_source_label(i);
    if ismember(tmp, scale250_subcortROIs)
        bool_temp(i)=1;
    end
end
clear ind
ind=find(bool_temp);

% remove subcortical ROIs
agr_source_data(:,ind)=[];
marked_fra_eigenvalues(ind)=[];
ave_source_coor(ind,:)=[];
ave_source_label(ind)=[];

% organize into the all sessions and subjects
ses=1;
corti_fra_eigenvalues{ses,subj,tr}=marked_fra_eigenvalues;
corti_ave_source_coor{ses,subj,tr}=ave_source_coor;
corti_ave_source_labl{ses,subj,tr}=ave_source_label;
cd /ssd/zhibin/1overf/Cleaned_sourcedata/cortical_source_data
% save('corti_fra_eigenvalues.mat','corti_fra_eigenvalues');
% save('corti_ave_source_coor.mat','corti_ave_source_coor');
% save('corti_ave_source_labl.mat','corti_ave_source_labl');
% save(['./' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat'],'agr_source_data')