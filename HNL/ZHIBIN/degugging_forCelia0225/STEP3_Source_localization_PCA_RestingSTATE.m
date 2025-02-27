% This is for my own reference: my earlier unorganized code can be found at:
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/organize_62stroke/step1_clean_data.m
% Created by Zhibin 4/6/2024
% #####################################################################################################################



% This script does the source localization (inverse solution) for the 256 channel EGI raw EEG data
% And run PCA to aggreate source data in each ROI of the 448 cortical areas
clear
DOSSIERS={'E:\DATASET\3_After_manual_rejection\Control','E:\DATASET\3_After_manual_rejection\patient'};
DOSSIERS_SAVE={'E:\DATASET\4_Source_reconstructed\Control','E:\DATASET\4_Source_reconstructed\patient'};
%% Construct brain model and create forward matrix in MNE
% Refer to MNE_construct_source.ipynb to create the leadfield matrix
% load forward matrix and source
load('/home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/base_files/MNE/EGI256/MNE_source_model.mat','leadfield','source_rr');
% leadfield: leadfield matrix from MNE (256 EEG channels X 5124 brain sources)
% source_rr: coordinates of the 5124 brain sources

%% Source localization 
% load the brain model, ROI names and localation files for the 463 ROIs
load('/home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/base_files/brain_ROIs.mat','Brain','roiNames_250','scale250_subcortROIs','parcels')
% Brain: The brain model of Lausanne2008_fsaverageDSsurf
% parcels: contain indices of 463 ROIs in a 3d brain matrix (256x256x256)
% roiNames_250: labels for 463 ROIs (448 cortical + 15 subcortical ROIs) (scale 250)
% scale250_subcortROIs: indices for 15 subcortical ROIs in roiNames_250

% Electrode for flip for left stroke 

LEFTIDX=[16;22;23;27;28;29;32;33;34;35;37;38;39;46;47;54;118;110;100;89;80;45;125;117;109;99;88;79;53;136;124;116;108;98;87;78;60;146;135;123;115;107;97;86;77;66;145;134;122;114;106;96;85;76;72;133;121;113;105;95;84;75;71;65;59;52;44;9;120;112;104;94;83;74;70;64;58;51;43;17;111;103;93;69;63;57;50;42;24;102;92;68;62;56;49;41;30;91;82;73;67;61;55;48;40;36;241;242;243;244;245;246;247;248;249;250;251;252;253;254;255;256];
RIGHTIDX=[7;14;6;20;13;5;25;19;12;4;18;11;3;10;2;1;127;128;129;130;131;132;138;139;140;141;142;143;144;148;149;150;151;152;153;154;155;156;157;158;159;160;161;162;163;164;165;166;167;168;169;170;171;172;173;174;175;176;177;178;179;180;181;182;183;184;185;186;187;188;189;190;191;192;193;194;195;196;197;198;199;200;201;202;203;204;205;206;207;208;209;210;211;212;213;214;215;216;217;218;219;220;221;222;223;224;238;239;240;234;235;236;237;230;231;232;233;226;225;227;228;229];

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

% Load the injured hemisphere

load('E:\DATASET\Script\Injured_hemi.mat')

% Boucle
for d=1:1:length(DOSSIERS)
save(DOSSIERS_SAVE{d},'corti_ave_source_coor','corti_ave_source_labl','corti_roiNames');
% corti_ave_source_coor: coordinates for 448 cortical ROIs
% corti_ave_source_labl: indices for 448 cortical ROIs in roiNames_250
% corti_roiNames: names of the 448 cortical ROIs

% get inverse matrix for source localization
% You will neded to add the following repo to your path:
% https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh/tree/main/Simulations/util
%addpath ../../../AdaptiveGraphicalLassoforParCoh/Simulations/util/ % this is my local path
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat] = inversemodel(leadfield,'prctile',1);

dossier1=DOSSIERS{1,d};
F=dir(fullfile(dossier1,'*.mat'));

% FILES=[62]

for f=1:1:length(F)

 fichier =  F(f).name;
 load(fullfile(F(f).folder,fichier));

% for f=1:1:length(FILES)
% 
% fichier =  F(FILES(f)).name;
% load(fullfile(F(FILES(f)).folder,fichier));

% Lateralisation EEG if left hemisphere is injured

if d==1
    D=1
elseif d==2
    D=3
end

index = find(cellfun(@(x) strcmp(x, subject_ID), Injured_hemi(:,D)));

if strcmp(Injured_hemi{index,D+1},'R')
    disp('no flip')
elseif strcmp(Injured_hemi{index,D+1},'L') %
    disp ('flip')

    for i = 1:length(LEFTIDX)
        % Sauvegarder temporairement la ligne de LEFTIDX
        tempRow = preprocessed_eeg(LEFTIDX(i), :);

        % Échanger les lignes
        preprocessed_eeg(LEFTIDX(i), :) = preprocessed_eeg(RIGHTIDX(i), :);
        preprocessed_eeg(RIGHTIDX(i), :) = tempRow;
    end
end

clear tempRow
% load the cleaned EEG file
% Navigate to your cleaned preprocessed EEG data directory, such as:
% cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_cleaned
% cd /ssd/zhibin/archive/EEG_stroke_62_cleaned
% load the preprocessed EEG file for one patient

display(['start processing subject file: ' fichier]);

% Lateralization


% Using inverse solution to covert preprocessed eeg to localized source data
cd /ssd/zhibin/archive/EEG_stroke_62_reorganized
load('2.mat')
source_data=inversemat*preprocessed_eeg;


%% PCA to aggregate cortical source data
% PCA to aggregate source data in each ROI
corti_source_data=[];
for sr=1:max(unique(source_labels))
    I=find(source_labels==sr); % need to load save and load this source_labels
    if ~isempty(I)
        [~, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
        corti_source_data=[corti_source_data SCORE(:,1)];
    end
end % 66s
    
% remove source data of the subcortical ROIs
% ind_rm=ind;
corti_source_data(:,ind_rm)=[];

% save the aggregated source data for the 448 cortical ROIs
% navegate the your destination directory, such as
% cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source

SaveFolder=DOSSIERS_SAVE{1,d};
cd   (SaveFolder)

save([num2str(subject_ID) '.mat'], ...
        'corti_source_data','corti_ave_source_coor','corti_ave_source_labl', ...
        'subject_ID','chanlocs','ch_labels','ch_dubious','ch_peripheral','Fs');

display(['Complete one file: ' fichier '********************'])

end
end

