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

%% WCA_WMSC Method


%% AGL Method
