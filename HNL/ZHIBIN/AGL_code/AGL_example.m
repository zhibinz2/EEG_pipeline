% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/organize_62stroke
clear

Fs=1000;
srnew = 200;
downsample = Fs/srnew;

passbands = [1 3; 3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29];
bandlabels = {'Delta','Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};

% make filter
attenuation=60;
filt_ds=cell(1,length(passbands));
for freqBand=1:6
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
addpath(genpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh'))
allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
n_Lambdas=length(allLambdas); % number of lambda values
min_LamdaIn = min(allLambdas);

% load one individual Reduced connectome for AGL
cd /home/zhibinz2/Documents/GitHub/STROKE_P61/getLesionMaskConnectome_nonflip_p61
load('Connectome_p61.mat', 'redConnectome_p61')
p=1;% pick one patient
redC_mat=squeeze(redConnectome_p61(p,:,:));
cd /home/zhibinz2/Documents/GitHub/STROKE_P61
load("source_info.mat")
SC_p=logical(redC_mat(source_roi_index,source_roi_index));% select the cortical rois
% cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
% load('corti_ave_source_labl.mat')
% SC_p=logical(redC_mat(corti_ave_source_labl,corti_ave_source_labl));


% penalty selection and fit precision
cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_corti_source
load('0.mat');% load the EEG source time series from one subject
downsample_data=resample(double(corti_source_data),1,downsample,'Dimension',1); % downsample 

% pick a frequency to process the data
freq=1; % delta
filtered_data = filter(filt_ds{freq},downsample_data);
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
datapermuted_cov=nan(n_split,n_sr,n_sr);

% compute covariance
stroke_Cov = cov(sourceDataReal);
for n=1:n_split
    datapermuted_cov(n,:,:) = cov([squeeze(datapermuted(n,:,:))]');
end

% AGL
% cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL
addpath(genpath('./AGL_util'));
dataCovs_op=squeeze(datapermuted_cov);
tic
[penalizationIn_op,penalizationOut_op,minDev_op]=penaltyselection( ...
SC_p,allLambdas,allLambdasOut,dataCovs_op);
toc % take 701 seconds for this example

tic
[stroke_Pcov] = fitprecision( ...
SC_p,penalizationIn_op,penalizationOut_op,min_LamdaIn,stroke_Cov);
toc % take 277 seconds for this example

% add this repo to your path (https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh)
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util

% convert real value covariance to complex matrix (448x448) 
% then compute coherence
stroke_coh=normalizeCSD(r2c(stroke_Cov)); % coherence
stroke_Pcoh=normalizeCSD(r2c(logical(stroke_Pcov)));% partial coherence using boolean
