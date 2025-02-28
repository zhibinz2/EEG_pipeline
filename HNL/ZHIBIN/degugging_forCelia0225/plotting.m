addpath(genpath('/home/zhibinz2/Documents/GitHub/BrainNetViewer_20191031'))
BrainNet
% load source power
load('/home/zhibinz2/Documents/GitHub/EEG_pipeline/Celia/trouble_shooting20250225/Patient_power.mat')
sourcedelta=POWER.DELTA.Sum;
sourcebeta=POWER.BETA.Sum;
sourcenames=POWER.INFO.Name;
load('/home/zhibinz2/Documents/GitHub/EEG_pipeline/Celia/trouble_shooting20250225/scalp256/Patient_256_Power.mat')
eegdelta=POWER.DELTA.Sum;
eegbeta=POWER.BETA.Sum;
eegnames=POWER.INFO.Name;
eegnames=eegnames(1,:);

for freq=1:5
    subplot(1,5,freq)
    


[isMatch, idx2] = ismember(sourcenames, eegnames);
idx1 = find(isMatch);  % Indices in cellArray1
idx2 = idx2(isMatch); 


sourcepowersum=sourcedelta;
eegpowersum=eegdelta;

sourcepowersum=sourcebeta;
eegpowersum=eegbeta;
for pt=1:100
    sourcedeltasum=sourcepowersum(:,idx1(pt));
    eegdeltasum=eegpowersum(:,idx2(pt));
    clf
    subplot(121);
    topoplot(eegdeltasum, chanlocs, 'electrodes', 'on');colorbar; colormap('jet');clim([min(eegdeltasum) max(eegdeltasum)])
    subplot(122);
    addpath /home/zhibinz2/Documents/GitHub/STROKE_P61/lesion_mask_on_L_p61_20240427/networkx/new_simu
    % value_display=std(corti_source_data);
    value_display=sourcedeltasum;
    roi3dbrain(value_display, X448_fs, Y448_fs, Z448_fs, ...
        pointsize*2,source_roi_index,ROI448_IC,ROI463_IC,triang,'jet');
    text(3,100,0,'I');text(170,100,0,'C');text(70,8,0,'Posterior')
    xlim([0 175]);ylim([0 205]);
    view(0,90); 
    [caz,cel] = view;
    view(caz,cel)
    clim([min(value_display) max(value_display)])
    % title(titles{nx})
    % subtitle([bandlabels{freq}])
    axis off
    sgtitle(num2str(pt))
    pause(1)
end

%delta 
pt=2;

%beta
pt=10 
%% plot in my own way
addpath(genpath('/home/zhibinz2/Documents/GitHub/eeglab'))
% load corti_source_data
cd /ssd/zhibin/archive/EEG_stroke_62_cleaned20240423
% figure;
% for p=1:60; % 12
% p=52; % remove bad channel from leadfield matrix
% p=12
% p=53 
% p= 5 
% p= 35
% p = 48 
% p =1
% p=23
p=34
load([num2str(p) '.mat']);
% reduced the size
wind=[91000:99000];
preprocessed_eeg=preprocessed_eeg(:,wind);
% preprocessed_eeg = zscore(preprocessed_eeg, 0, 2); % Normalize across each channel
% plot(wind,preprocessed_eeg);
% preprocessed_eeg=preprocessed_eeg;

% Define delta band range
% freq_band = [0.5 3]; % Hz
freq_band = [3 7]; % Hz
freq_band = [8 12]; % Hz
freq_band = [13 20]; % Hz
freq_band = [20 30]; % Hz
% Sampling rate
fs = 1000; % Hz

% % Compute power spectral density (PSD) using Welch's method
% [pxx, f] = pwelch(preprocessed_eeg', [], [], [], fs); % Transpose for pwelch
% % Find indices corresponding to delta band
% freq_idx = (f >= freq_band(1)) & (f <= freq_band(2));
% % Sum power in delta band for each channel
% freq_power = sum(pxx(freq_idx, :), 1); % Sum across frequency range
% % Plot topography using EEGLAB's topoplot

% [b, a] = butter(4, freq_band / (fs / 2), 'bandpass');
% filteredEEG = filtfilt(b, a, preprocessed_eeg')';
% Variance = var(filteredEEG, 0, 2);

N = size(preprocessed_eeg, 2);  % Number of samples
EEG_fft = fft(preprocessed_eeg, [], 2);
freqs = (0:N-1) * (fs / N);
freq_idx = find(freqs >= freq_band(1) & freqs <= freq_band(2));
powerSpectrum = abs(EEG_fft).^2 / N;
freq_power = sum(powerSpectrum(:, freq_idx), 2);

% figure;
% subplot(6,10,p)
% clf
topoplot(freq_power, chanlocs,'electrodes', 'on');colorbar;colormap('jet');
% topoplot(Variance, chanlocs,'electrodes', 'on');colorbar;colormap('jet');
% clim([min(freq_power) max(freq_power)]);title('power')
% topoplot(freq_power, chanlocs, 'electrodes', 'labelpoint');
% topoplot(std(preprocessed_eeg'), chanlocs, 'electrodes', 'labelpoint');
topoplot(std(preprocessed_eeg'), chanlocs, 'electrodes', 'off');
% colorbar;
title(num2str(p))
% % title('Delta Band Power (0.5-4 Hz)');
 title('Theta');
% pause(0.25);
end
sgtitle('std 91000-99000')

figure;
clf
% topoplot(freq_power, chanlocs,'electrodes', 'on');
topoplot(std(preprocessed_eeg'), chanlocs, 'electrodes', 'labelpoint');
colorbar;
title(num2str(p))
% title('Delta Band Power (0.5-4 Hz)');
% end


% interpolate zero channels
stdVals = std(preprocessed_eeg, 0, 2);
zeroStdIdx = find(stdVals == 0);
validChannels = find(stdVals > 0);
% Interpolate each zero-variance channel using nearest valid neighbors
for i = 1:length(zeroStdIdx)
    ch = zeroStdIdx(i);
    % Find the nearest valid neighboring channels
    lower = max(validChannels(validChannels < ch));
    upper = min(validChannels(validChannels > ch));
    if isempty(lower)
        % If no lower channel, use only the upper channel
        preprocessed_eeg(ch, :) = preprocessed_eeg(upper, :);
    elseif isempty(upper)
        % If no upper channel, use only the lower channel
        preprocessed_eeg(ch, :) = preprocessed_eeg(lower, :);
    else
        % If both neighbors exist, use average for interpolation
        preprocessed_eeg(ch, :) = (preprocessed_eeg(lower, :) + preprocessed_eeg(upper, :)) / 2;
    end
end



%leadfield remove badchans
badchans=[];badchansplus=[];
badchans=unique([badchans badchansplus'])
% badchans=ch_dubious;
%load leadfield
load('/home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/base_files/MNE/EGI256/MNE_source_model.mat')
leadfield(badchans,:)=zeros;

% source
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat] = inversemodel(leadfield,'prctile',50);
source_data=inversemat*preprocessed_eeg;
% figure;plotx(std(source_data'))

% correlation
reconEEG=leadfield*source_data;
corrmat=corr(reconEEG',preprocessed_eeg');
% imagesc(corrmat);colorbar
% figure;clf;topoplot(diag(corrmat), chanlocs,'electrodes', 'on');colorbar
figure; topoplot(diag(corrmat), chanlocs,'electrodes', 'numbers');colorbar
plot(1:size(preprocessed_eeg,1),diag(corrmat))
badchansplus=find(diag(corrmat)<0.9)
plot(1:length(wind),preprocessed_eeg);

% Compute power spectral density (PSD) using Welch's method
[pxx, f] = pwelch(reconEEG', [], [], [], fs); % Transpose for pwelch
freq_idx = (f >= freq_band(1)) & (f <= freq_band(2));
freq_power = sum(pxx(freq_idx, :), 1); % Sum across frequency range
% Plot topography using EEGLAB's topoplot
figure;
clf
topoplot(freq_power, chanlocs,'electrodes', 'on');colorbar;colormap('jet');clim([min(freq_power) max(freq_power)]);title('power')
% topoplot(std(reconEEG'), chanlocs,'electrodes', 'labelpoint');
colorbar;
title(num2str(p))
% title('Delta Band Power (0.5-4 Hz)');
% pause(0.5);


% % Compute power spectral density (PSD) using Welch's method
% [pxx, f] = pwelch(source_data', [], [], [], fs); % Transpose for pwelch
% % Find indices corresponding to delta band
% freq_idx = (f >= freq_band(1)) & (f <= freq_band(2));
% % Sum power in delta band for each channel
% freq_power = sum(pxx(freq_idx, :), 1); % Sum across frequency range


N = size(source_data, 2);  % Number of samples
EEG_fft = fft(source_data, [], 2);
freqs = (0:N-1) * (fs / N);
freq_idx = find(freqs >= freq_band(1) & freqs <= freq_band(2));
powerSpectrum = abs(EEG_fft).^2 / N;
freq_power = sum(powerSpectrum(:, freq_idx), 2);


% cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/base_files/MNE/EGI256
% source_rrX=source_rr(:,1);source_rrY=source_rr(:,2);source_rrZ=source_rr(:,3);
addpath /home/zhibinz2/Documents/GitHub/STROKE_P61/lesion_mask_on_L_p61_20240427/networkx/new_simu
figure
values_display=freq_power';
% values_display=std(source_data');
s=scatter3(source_x, source_y, source_z, pointsize, values_display, ...
        "filled", 'MarkerFaceAlpha',0.9);
xlabel('x');ylabel('y');zlabel('z');
view(0,90); 
[caz,cel] = view;
view(caz,cel)
colormap("cool")
colorbar; 
% clim([0 20])



% PCA to aggregate source data in each ROI
corti_source_data=[]; LATENT_DATA=[];
for sr=1:max(unique(source_labels))
    I=find(source_labels==sr); % need to load save and load this source_labels
    if ~isempty(I)
        [coeff, SCORE, LATENT] = pca(source_data(I,:)','Centered',false);
        LATENT_DATA=[LATENT_DATA LATENT(1)];
        corti_source_data=[corti_source_data SCORE(:,1)];
    end
end % 66s
% remove source data of the subcortical ROIs
corti_source_data(:,ind_rm)=[];


sum(coeff.^2)
    
figure;plot(var(source_data'))
figure;plot(var(corti_source_data))
figure;plot(LATENT_DATA)

% load source aggreagate
% cd /ssd/zhibin/archive/flip_EEG_stroke_62_corti_source20240424
% clear corti_source_data
% load([num2str(p) '.mat']);

% corti_source_data=corti_source_data(wind',:);

% Compute power spectral density (PSD) using Welch's method
[pxx, f] = pwelch(corti_source_data, [], [], [], fs); % Transpose for pwelch
% Find indices corresponding to delta band
freq_idx = (f >= freq_band(1)) & (f <= freq_band(2));
% Sum power in delta band for each channel
freq_power = sum(pxx(freq_idx, :), 1); % Sum across frequency range

N = size(corti_source_data', 2);  % Number of samples
EEG_fft = fft(corti_source_data', [], 2);
freqs = (0:N-1) * (fs / N);
freq_idx = find(freqs >= freq_band(1) & freqs <= freq_band(2));
powerSpectrum = abs(EEG_fft).^2 / N;
freq_power = sum(powerSpectrum(:, freq_idx), 2);



addpath /home/zhibinz2/Documents/GitHub/STROKE_P61/lesion_mask_on_L_p61_20240427/networkx/new_simu
figure
% value_display=std(corti_source_data);
value_display=freq_power;
roi3dbrain(value_display, X448_fs, Y448_fs, Z448_fs, ...
    pointsize*2,source_roi_index,ROI448_IC,ROI463_IC,triang,'sky');
text(3,100,0,'I');text(170,100,0,'C');text(70,8,0,'Posterior')
xlim([0 175]);ylim([0 205]);
view(0,90); 
[caz,cel] = view;
view(caz,cel)
clim([0 max(value_display)])
% title(titles{nx})
% subtitle([bandlabels{freq}])
axis off

%
