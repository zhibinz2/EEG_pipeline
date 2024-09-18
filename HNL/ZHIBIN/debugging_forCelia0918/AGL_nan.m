% This is for my own reference: my earlier unorganized code can be found at:
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/organize_62stroke/step1_clean_data.m
% open /home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/EGI_256chan/step2_source_localization_PCA.m
% open /home/zhibinz2/Documents/GitHub/EEG_pipeline/Celia/Trouble_shooting/wetransfer_calt-mat_2024-09-08_1921/STEP4_Quantification_coherence.m
% open /home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/AGL_code/AGL_example.m
% Created by Zhibin 4/6/2024
% #####################################################################################################################

% This script does the source localization (inverse solution) for the 256 channel EGI raw EEG data
% And run PCA to aggreate source data in each ROI of the 448 cortical areas
clear
%% Construct brain model and create forward matrix in MNE
% Refer to ../EGI_256chan/MNE_construct_source.ipynb to create the leadfield matrix
% load forward matrix and source
load('../base_files/MNE/EGI256/MNE_source_model.mat','leadfield','source_rr');
% leadfield: leadfield matrix from MNE (256 EEG channels X 5124 brain sources)
% source_rr: coordinates of the 5124 brain sources

%% Source localization 
% load the brain model, ROI names and localation files for the 463 ROIs
load('../base_files/brain_ROIs.mat','Brain','roiNames_250','scale250_subcortROIs','parcels')
% Brain: The brain model of Lausanne2008_fsaverageDSsurf
% parcels: contain indices of 463 ROIs in a 3d brain matrix (256x256x256)
% roiNames_250: labels for 463 ROIs (448 cortical + 15 subcortical ROIs) (scale 250)
% scale250_subcortROIs: indices for 15 subcortical ROIs in roiNames_250

%% Source localization 
% load the brain model, ROI names and localation files for the 463 ROIs
load('../base_files/brain_ROIs.mat','Brain','roiNames_250','scale250_subcortROIs','parcels')
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

% corti_ave_source_coor: coordinates for 448 cortical ROIs
% corti_ave_source_labl: indices for 448 cortical ROIs in roiNames_250
% corti_roiNames: names of the 448 cortical ROIs

% get inverse matrix for source localization
% You will neded to add the following repo to your path:
% https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh/tree/main/Simulations/util
% addpath(genpath('../../../AdaptiveGraphicalLassoforParCoh/Simulations/util/')); % this is my local path
cd /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat] = inversemodel(leadfield,'prctile',1);

%% load Celia's source data
cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/Celia/Trouble_shooting/
load('ADLE.mat','corti_source_data');
AGL_example_data=corti_source_data;

% figure;
% plot(1:length(AGL_example_data),AGL_example_data);
% isnan(sum(AGL_example_data,'all'))

%% AGL
cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/AGL_code/
Fs=1000;
srnew = 200;
downsample = Fs/srnew;

% Celia's filter bands
passbands = [1 4; 4.5 8; 8.5 12; 12.5 30; 30.5 50];
bandlabels = {'Delta','Theta', 'Alpha', 'Beta', 'Gamma'};

% passbands = [1 3; 3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29];
% bandlabels = {'Delta','Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};

% make filter
attenuation=60;
filt_ds=cell(1,length(passbands));
for freqBand=1:length(passbands)
    % Select frequency
    passFreq1 = passbands(freqBand,1);
    passFreq2 = passbands(freqBand,2);
    filt_d = designfilt('bandpassiir','FilterOrder',20, ...
    'PassbandFrequency1',passFreq1,'PassbandFrequency2',passFreq2, ...
    'StopbandAttenuation1',attenuation,'PassbandRipple',0.2, ...
    'StopbandAttenuation2',attenuation,'SampleRate',srnew);
    filt_ds{freqBand}=filt_d;
end

% AGL variables
% https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh
addpath(genpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh'))
allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
n_Lambdas=length(allLambdas); % number of lambda values
min_LamdaIn = min(allLambdas);

% load a random individual connectome for AGL
load("AGL_example_data.mat",'Individual_Connectome')

stroke_pcoh5=nan(5,448,448);
% pick a frequency to process the data
for freq=1:5%; % beta2
    filtered_data = filter(filt_ds{freq},AGL_example_data);
    hilbertdata = hilbert(filtered_data');
    % combined real and imaginary part
    sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));
    sourceDataReal = [sourceDataReal*(1/mean(abs(sourceDataReal(:))))]'; % normalize data
    % split into 2 ensambles: 2 x #source x #(samples/2)
    n_split=2; n_sr=size(sourceDataReal,2);
    sam_len=size(sourceDataReal,1);
    sam_size=floor(sam_len/n_split); sam_range=1:sam_size;
    sourceDataReal = sourceDataReal(1:n_split*sam_size,:);
    datareshaped = reshape(sourceDataReal, sam_size, n_split, n_sr);
    datapermuted = permute(datareshaped,[2,3,1]); % split into 2 ensambles: 2x896x9300
    % compute covariance
    datapermuted_cov=nan(n_split,n_sr,n_sr);
    stroke_Cov = cov(sourceDataReal);
    for n=1:n_split
        datapermuted_cov(n,:,:) = cov([squeeze(datapermuted(n,:,:))]');
    end
    
    % AGL
    % penalty selection and fit precision
    % cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
    addpath(genpath('./AGL_util'));
    dataCovs_op=squeeze(datapermuted_cov);
    tic
    [penalizationIn_op,penalizationOut_op,minDev_op]=penaltyselection( ...
    Individual_Connectome,allLambdas,allLambdasOut,dataCovs_op);
    toc % take 701 seconds for this example
    
    tic
    [stroke_Pcov] = fitprecision( ...
    Individual_Connectome,penalizationIn_op,penalizationOut_op,min_LamdaIn,double(stroke_Cov));
    toc % take 277 seconds for this example
    
    % add this repo to your path (https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh)
    addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
    
    % convert real value covariance to complex matrix (448x448) 
    % then compute coherence
    stroke_coh=normalizeCSD(r2c(stroke_Cov)); % ordinary coherence
    stroke_Pcoh=normalizeCSD(r2c(logical(stroke_Pcov)));% partial coherence using boolean
    
    display(['freq' num2str(freq)])
    sum(isnan(stroke_Pcoh),'all')
    
    stroke_pcoh5(freq,:,:)=stroke_Pcoh;
end
figure
for freq=1:5
    subplot(2,5,freq)
    imagesc(squeeze(stroke_pcoh5(freq,:,:)))
    colorbar
    title([bandlabels{freq} '   nan#: ' num2str(sum(isnan(stroke_pcoh5(freq,:,:)),'all'))])
    subplot(2,5,5+freq)
    imagesc(logical(squeeze(stroke_pcoh5(freq,:,:))))
    colorbar
    title([bandlabels{freq} '   binary: edges # ' num2str(sum(triu(squeeze(logical(stroke_pcoh5(freq,:,:))),1),'all'))])
end
sgtitle('ADLE - partial coherence (EEG)')
save('output_debug.mat','stroke_pcoh5'); % reconstruire BUG pour patient file

figure;
imagesc(triu(squeeze(logical(stroke_Pcoh)),1));colorbar

