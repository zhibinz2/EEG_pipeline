% (my reference) cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/hilbert2cov.m

% This is a script for computing coherence and partial coherence using AGL with the hyperscanning cortical source data.
clear
cd ../cortical_source_data/
load('corti_ave_source_labl.mat')
load('corti_ave_source_coor.mat')
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
% every file has the same labels and coordinations
% every trial is the same, load one trial
source_labels=corti_ave_source_labl{1,1,1}; 
source_coor=corti_ave_source_coor{1,1,1};

Fs=2000;
srnew = 200;
downsample = Fs/srnew;

passbands = [3.5 6.5; 7 10; 10.5 13.5; 14 20; 21 29];
bandlabels = {'Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2'};

% make filter
attenuation=60;
filt_ds=cell(1,length(passbands));
for freqBand=1:5
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
% add this repo to your path (https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh)
addpath(genpath('../GitHub/AdaptiveGraphicalLassoforParCoh'))
allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
n_Lambdas=length(allLambdas); % number of lambda values
min_LamdaIn = min(allLambdas);

% load structure connectome for AGL
load('scale250_Connectome.mat','fc'); % Virtual-Tractography/ForZhibin/processed_data/scale250_Connectome.mat
SC=logical(fc(source_labels,source_labels));

%% select one file to process
ses=12; % sesssion 12
subj=2; % subject 2 (subject R)
tr=12; % trial 12
load(['./cortical_source_data/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat']);
% downsample to 200 Hz
downsample_data=resample(double(agr_source_data),1,downsample,'Dimension',1);

% compute covariance 
for ses=1:12
    tic
    for subj=1:2
        for tr=1:12
            load(['./cortical_source_data/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat']);

            hilbert_dataCov=cell(1,5);
            for freq=1:5
                
                % resample to 200 Hz
                downsample_data=resample(double(agr_source_data),1,downsample,'Dimension',1);
                filterd_data = filter(filt_ds{freq},downsample_data);
                hilbertdata = hilbert(filterd_data'); 

                % combine real and imaginary parts 
                sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));
                sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data

                hilbert_dataCov{freq} = cov(sourceDataReal');

                save(['./hilbert_datacov/' num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(tr) '.mat'],'hilbert_dataCov');
            end
        end
    end
    toc
end
% This will take one hour

% organize all covariance matricies together with resorted trial sequence
hilbert_dataCov_all=nan(12,2,12,5,896,896);
for ses=1:numSes
    tic
    clear conditions sortorders
    runid = num2str(seeds(ses,:));
    load(['../../Cleaned_data/clean_' runid '.mat'],'conditions'); % load the trial condition sequence
    % sort order
    [x,sortorder]=sort(conditions); % sort the trial condition sequence from 1 to 4
    for subj=1:2
        for tr=1:12
            % load the previous covariance
            load(['../../Cleaned_data/hilbert_datacov/' ...
                    num2str(seeds(ses,:)) '/subj' num2str(subj) '_tr_' num2str(sortorder(tr)) '.mat'] ...
                    ,'hilbert_dataCov');
            for freq=1:5
                hilbert_dataCov_all(ses,subj,tr,freq,:,:)=hilbert_dataCov{freq};
            end
        end
    end
    toc
end
% This will take one minute

% convert real value covariance to complex values for each trial
addpath(genpath('./AGL_util'));
Complex_dataCov_all=nan(12,2,12,5,448,448);
for ses=1:numSes
    for subj=1:2
        for tr=1:12
            for freq=1:5
                Complex_dataCov_all(ses,subj,tr,freq,:,:)=r2c(squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:)));
            end
        end
    end
end

% normalize to real value coherence (now they are back to 448 x 448) for each trial
coh_all=nan(12,2,12,5,448,448);
for ses=1:numSes
    for subj=1:2
        for tr=1:12
            for freq=1:5
                coh_all(ses,subj,tr,freq,:,:)=normalizeCSD(squeeze(Complex_dataCov_all(ses,subj,tr,freq,:,:)));
            end
        end
    end
end

%% AGL (method option 3 by subject, in order to speed up the process)
% For this option: Choose the penalties for each subject in each frequency band for all conditions.
% Each subject produces 24 trials in 2 session, 3 trials of each condition in each session.  
% So take 1 trial of each condition (8 trials from 2 session, 4 from each session) and average the covariance.  
% Now you have 3 average covariance matrices.  Then run the penalty selection function.
% Go back and fit each trial for that subject using that penalty.
% This way you would have to do 12 x 5 penality optimizations.  
% I think this option may be better, because the subject to subject differences may be greater than the conditon differences.

% organize the data into 3 ensambles for each subject
ave_hilcov_option3=nan(3,2,6,5,896,896); % 3 ensamble x 2 subjects x 6 double-sessions x 5 frequency x 894 x 894
for ensam = 1:3
    for subj =1:2
        for dl_ses=1:6
            for freq=1:5
                tic
                ave_hilcov_option3(ensam,subj,dl_ses,freq,:,:)=...
                    squeeze(mean(squeeze(mean(squeeze(hilbert_dataCov_all([1+2*(dl_ses-1):2+2*(dl_ses-1)], ...
                    subj,[1:3:12]+(ensam-1),freq,:,:)),2)),1));
                toc
            end
        end
    end
end

% penaltyselection for each subject
% add this repo to your path (https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh)
addpath(genpath('/../../AdaptiveGraphicalLassoforParCoh/Simulations/util'))
addpath(genpath('../../AdaptiveGraphicalLassoforParCoh/AGL'))
penalizationIn_op3=nan(2,6,5);
penalizationOut_op3=nan(2,6,5);
minDev_op3=nan(2,6,5); % We should also save the deviance values for the final fit, because they tell us how good a model it is.
parfor subj=1:2
    for dl_ses=1:6
        for freq=1:5
            tic
            dataCovs_op=squeeze(ave_hilcov_option3(:,subj,dl_ses,freq,:,:));
            [penalizationIn_op3(subj,dl_ses,freq),penalizationOut_op3(subj,dl_ses,freq),minDev_op3(subj,dl_ses,freq)]=...
            penaltyselection(SC,allLambdas,allLambdasOut,dataCovs_op);
            toc
        end
    end
end
% This step will take 3 days

% fitprecision for each trial
X_op3=nan(12,2,12,5,896,896);
for subj=1:2
    for dl_ses=1:6
        for freq=1:5
            penalizationIn=penalizationIn_op3(subj,dl_ses,freq);
            penalizationOut=penalizationOut_op3(subj,dl_ses,freq);
            for tr=1:12
                tic
                ses=1+2*(dl_ses-1); % synchronization session
                dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
                [X_op3(ses,subj,tr,freq,:,:)] = fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);
                ses=2+2*(dl_ses-1); % syncopation session
                dataCov=squeeze(hilbert_dataCov_all(ses,subj,tr,freq,:,:));
                [X_op3(ses,subj,tr,freq,:,:)] = fitprecision(SC,penalizationIn,penalizationOut,min_LamdaIn,dataCov);
                toc
            end
        end
    end
end

% convert the X_op3 (896x896) to complex values (448x448) using r2c function
Complex_X_op3_all=nan(12,2,12,5,448,448);
tic
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                X=squeeze(X_op3(ses,subj,tr,freq,:,:));
                Complex_X_op3_all(ses,subj,tr,freq,:,:)=r2c(X);
            end
        end
    end
end
toc 
% This step takes 33 s 

% normalize Complex_X_op3_all to real value partial coherence
Pcoh_all=nan(12,2,12,5,448,448);
tic
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                cov_c=squeeze(Complex_X_op3_all(ses,subj,tr,freq,:,:));
                newG=normalizeCSD(cov_c);
                Pcoh_all(ses,subj,tr,freq,:,:)=newG;
            end
        end
    end
end
toc % This will take 12 seconds

