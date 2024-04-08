% This is for my own reference: my earlier unorganized code can be found at:
% open /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/organize_62stroke/step1_clean_data.m
% Created by Zhibin 4/3/2024
% #####################################################################################################################


% This script process the 256 channel EGI raw EEG data
%% some basic settings
clear
% These are 49 channels on the pheripheral and Cz chanel (#257)
ch_peripheral_cz=[241 242 243 238 239 240 ...
    244 245 246 247 251 256 91 102 111 120 133 145 165 174 187 199 208 216 229 233 237 236 235 234 ...
    232 228 217 209 200 188 175 166 156 146 134 121 112 103 92 82 255 250 ...
    257];
% We will remove these pherical channels, as they contain lots EMG artifect
% Cz channel is the default reference on the EGI device, which has values of zeros on the raw data
% After removal of 49 peripheral channels and the Cz, there are 208 channels left (257-49=208)

% remove Cz channel from ch_peripheral_cz 
ch_peripheral=ch_peripheral_cz;
ch_peripheral(end)=[]; 

% original sampling frequency
Fs=1000;

% channel location information (coordinates)
chanlocs = load('./base_files/MNE/EGI256/chanlocs.mat');
chanlocs = chanlocs.chanlocs;
% we have to keep only 256 channel to meet MNE's data format for the building of brain model 
chanlocs(257)=[]; % Channel 257 is Cz, which is default device reference (values are all zeros)

% Add the original channel sequence index into channel information
for ch=1:256 
    chanlocs(ch).urchan=ch;
end

% Create channel text labels
ch_labels = cell(256,1);
for c = 1:256
    ch_labels{c}=num2str(c);
end
clear c ch
% label some landmark channels
ch_labels{18}='Fp2'; ch_labels{37}='Fp1';
ch_labels{36}='F3'; ch_labels{224}='F4';
ch_labels{59}='C3'; ch_labels{183}='C4';
ch_labels{69}='T7'; ch_labels{202}='T8';
ch_labels{87}='P3'; ch_labels{153}='P4';
ch_labels{96}='P7'; ch_labels{170}='P8';
ch_labels{94}='LM'; ch_labels{190}='RM';
ch_labels{116}='O1'; ch_labels{150}='O2';
ch_labels{31}='NAS'; 
ch_labels{21}='Fz'; ch_labels{101}='Pz'; ch_labels{126}='Oz';
% ch_labels{257}='Cz';

% Create another set of channel location information with only the central 208 channels
% by removing 49 perpheral channels (then to check channel index for eyechans later)
central_chanlocs=chanlocs;
central_chanlocs(ch_peripheral)=[];
% use "urchan" to obtain an array of the orginal index for these central_chanlocs
urchan=extractfield(central_chanlocs,'urchan');

% setting the eye channels for ICA
eyechans = [248 252 253 67 249 254 73 230 226 225 219 231 227 218]; % channels around the eyes in the original 256 indices
eyechans=sort(eyechans);
% the above eyechans  67    73   218   219   225   226   227   230   231   248   249   252   253   254
% corresponds to      67    73   192   193   199   200   201   202   203   204   205   206   207   208 
% in the updated chanel index of central_chanlocs (eyechans_ind below)

% Identify the updated indices in central_chanlocs for these eyechan channels 
eyechans_ind=nan(1,length(eyechans));
for c=1:length(eyechans)
    eyechans_ind(c)=find(urchan==eyechans(c));
end

%% filter settings (band pass 0.25 - 50 Hz)
% make a high pass filter
Hp = makefilter(Fs,0.25,0.01,6,20,0); 
% make a low pass filter
Lp = makefilter(Fs,50,51,6,20,0);  

%% setting the thresholds for artifact detection and removal
% You can adjust these parameters as needed.
% Identify bad channel and epochs by using standard deviation 
var_threshold = 2.5;  % normalized variance threshold to reject trials.
chan_threshold = 2.5;  % to mark bad channels (smaller values are stricter)
% Identify eye movement component using correlation coefficients
corr_threshold = 0.4;  % threshold for identifying whether an ICA component contains eye movement. (smaller values are stricter)


%% Select the raw data of one patient: filter, clean and run ICA
subj_files=[0:1:60];
f=3; % select one patient (f=1:61)

tic
% Navigate the data directory such as
% cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_reorganized
load([num2str(subj_files(f)) '.mat']);
subject_ID=subj_files(f);
display(['start processing subject file: ' num2str(subj_files(f)) '.mat']);

% average re-reference
rData=Data-ones(size(Data,1),1)*mean(Data,1);
% remove 49 peripheral channels and cz reference
rData(ch_peripheral_cz,:)=[];
rData=rData';
% detrend the EEG data 
detrend_data=detrend(rData,1);
% add paddings to allow some buffer zone of edge effect
padding=zeros(round(size(detrend_data,1)/10), size(detrend_data,2));
detrend_pad=cat(1,padding,detrend_data,padding);
% high pass filter
pad_hp=filtfilthd(Hp,detrend_pad);
% Low pass filter (take a few seconds)
pad_lp=filtfilthd(Lp,pad_hp);
% remove the paddings that contain the edge artifacts
filtered_data=pad_lp((size(padding,1)+1):(size(padding,1)+size(detrend_data,1)),:);

% remove some variables to clean up some space on the memory 
clearvars Data rData detrend_data detrend_pad pad_hp pad_lp padding
% keep the output filtered_data

% organize into one second epochs
nepochs=floor(length(filtered_data)/Fs);
nchans=size(filtered_data,2);
epochdata=zeros(nepochs,Fs,nchans);
for e = 1: nepochs
    epochdata(e,:,:)=filtered_data((1+(e-1)*Fs):(Fs+(e-1)*Fs),:);
end
clearvars filtered_data e 
% keep the output epochdata

% identify bad chan and epochs by standard deviation
% compute standard deviations
eegstd=squeeze(std(epochdata,[],1));
chanstd=sum(eegstd,2);
epochstd=sum(eegstd,1);
    
% Threshhold epochs by standard deviation criteria and remove them.
% In principle you could threshold channels this way too.  
% But, I think with 32 channels you need to avoid that.  
% With 128 or more chanels you could.
badchan = find(epochstd / median(epochstd) > chan_threshold);
goodchan = setdiff(1:nchans,badchan);
    
badepoch = find(chanstd / median(chanstd) > var_threshold);
goodepoch = setdiff(1:nepochs, badepoch);
    
% remove bad epochs and connect back to 2 dimenstion matrix (time x chan)
good_epochdata=squeeze(epochdata(goodepoch,:,:));
good_filtered_data=zeros(length(goodepoch)*Fs,nchans);
for n = 1:length(goodepoch)
    good_filtered_data(((1+(n-1)*Fs):(Fs+(n-1)*Fs)),:)=squeeze(good_epochdata(n,:,:));
end

clearvars epochdata epochstd eegstd cahnstd good_epochdata n
% keep the output good_filtered_data

% run ICA
[icasig, A, W] = fastica(good_filtered_data');  % a few minutes

% Combine badchan and eyechans into dubious chans 
% to compute the correlation of each component with these dubious contaminated channels
dubious_chans =unique([badchan eyechans_ind]);

% Compute correlation with dubious chans (It take a few seconds)
corrs=zeros(length(dubious_chans),size(icasig,1));
for du = 1:length(dubious_chans)
    corrs(du, :)=corr(icasig',good_filtered_data(:,dubious_chans(du)));
end
corrsmax = max(corrs);

% Detect which components are not too correlated with dubious channels
badcomponents = abs(corrsmax) >= corr_threshold; 
display(['bad components: ' num2str(find(badcomponents))] )% display the bad components
goodcomponents = abs(corrsmax) < corr_threshold;

% Further detect bad ones in the goodcomponents by examing the weight
% proportion in A (mixing matrix)
proportion_threshold=0.6; % if any cahnnel weighted higher than 0.6 in A (This value can be adjusted as needed)
chancomponents = zeros(size(icasig,1),1);
B = zeros(nchans,size(icasig,1)); % 208 channel x 208 components
for n = 1:size(icasig,1) % loop through all 208 components
    B(:,n)=A(:,n).^2 / sum(A(:,n).^2);
    chancomponents(n)=max(B(:,n));
    if chancomponents(n) > proportion_threshold
        display(['detect bad component - ' num2str(n)])
        goodcomponents(n)=0; % remove it from good components
    end
end
    
% Restore the data without the bad components.
mixedsig=A(:,goodcomponents)*icasig(goodcomponents,:);  

% Inserting the periphereal channels back in place with zeros
% Because we need to use full 256 channel in oredr to import into MNE
preprocessed_eeg = zeros(256, size(mixedsig, 2));
i = 1;
for c = 1:256
    if ~ismember(c, ch_peripheral)
        preprocessed_eeg(c, :) = mixedsig(i, :);
        i = i + 1;
    end
end
clearvars icasig A W corrsmax du  n  mixedsig good_filtered_data
% Keep the output preprocessed_eeg

% Get the original index of the dubious chans in the 256 index
ch_dubious=nan(1,length(dubious_chans));
for i = 1:length(dubious_chans)
    ch_dubious(i) = central_chanlocs(dubious_chans(i)).urchan;
end

% Save the cleaned data at your destination directory such as 
% cd /home/zhibinz2/Documents/GitHub/archive/EEG_stroke_62_cleaned/
save([num2str(subject_ID) '.mat'],'preprocessed_eeg','Fs','ch_dubious','ch_peripheral','ch_labels','chanlocs','subject_ID')

display(['Complete one file: ' num2str(subj_files(f)) '.mat ********************'])
toc

