%% Load data
cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/Celia
load('P_RODJ20240426_After_Reject.mat')
Grip=Data_cut.Grip;
GripL=Grip.L;
GripR=Grip.R;

MVC=Data_cut.MVC;
Flex=MVC.Flex;
FlexL=Flex.L;FlexL_Raw=FlexL.Raw;FlexL_Norm=FlexL.Norm;
FlexR=Flex.R;FlexR_Raw=FlexR.Raw;FlexR_Norm=FlexR.Norm;

Info=Data_cut.Info;
Labels=Info.Labels;
Paretic=Info.Paretic;
Reject=Info.Reject; RejectL=Reject.L;RejectR=Reject.R;

% randomly take a 5-seconds contraction period 
% and concatenate the EEG and EMG into a matrix of "31x10000"

% In Grip, the data corresponding to the required grip task
% Performed with right and left limbs
% Data is clean and ready to use. 
% It's segmented 2 sec before the stimulus, 
% during the 5secondes of the displaying of the stimulus 
% and 3 seconds after the end of the displaying 
% (=10sec with a sampling frequency of 500Hz)
Fs=500;

% The labels are in Info.Labels and the paretic limb in Info.Paretic, 
% but the data is already flipped to correspond to a lesion in the right hemisphere and a left paretic limb.
% So for exemple: for the left movement 
% Channel 1 to 29 is EEG
% Channel 33 is agonist muscle
% Channel 34 is Antagonist muscle
EEGchs=1:29;
L_Flex=33;
L_Exte=34;


data=GripL.Data; % pick data from L (Lesion side)
tr=10; % pick a trial

figure;
subplot(211);
plot(1:5000,squeeze(data(EEGchs,:,tr))');
subplot(212);
plot(1:5000,squeeze(data([L_Flex L_Exte],:,tr))');

%% WCA_WMSC Method (over time)
% pick 10 trials
trs=[1:10];
trs_data=data([EEGchs L_Flex L_Exte],:,trs);

rate=500;

% pick a frequency to investigate
freq=2;
% select the time points to compute network measurements
ts=(rate:rate:size(s1,2));

% loop through all 31 channels in your adjacency matrix of 31x31
WMSC_mat=nan(10,31,31); % 10 time proints x 31 channels  x 31 channels
tic
for i=1:31 
    for j=1:31
        s1=squeeze(trs_data(i,:,trs))';
        s2=squeeze(trs_data(j,:,trs))';

        cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/cortical_muscular_coherece_forCelia
        [WMSC,t,freqs] = WCA_WMSC(s1,s2,rate); % You have to fix the freq range in the function that use your WCA code. I have no idea.
      
        WMSC_tmp=WMSC(freq,ts);

        WMSC_mat(:,i,j)=WMSC_tmp';
    end
end
toc % it takes about 30s

%% HNL Method
Fs=500;
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

% pick a frequency to process the data
freq=1; % delta
% select one trial
tr=1;
example_data=data([EEGchs L_Flex L_Exte],:,tr)';
filtered_data = filter(filt_ds{freq},example_data);
hilbertdata = hilbert(filtered_data');
% combined real and imaginary part
sourceDataReal = cat(1,real(hilbertdata),imag(hilbertdata));
sourceDataReal = [sourceDataReal*(1/mean(abs(sourceDataReal(:))))]'; % normalize data
% compute covariance
stroke_Cov = cov(sourceDataReal);
% https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh
addpath(genpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh'))
% convert real value covariance to complex matrix (448x448) 
% then compute coherence
stroke_coh=normalizeCSD(r2c(stroke_Cov)); % ordinary coherence

