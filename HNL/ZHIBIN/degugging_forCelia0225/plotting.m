addpath(genpath('/home/zhibinz2/Documents/GitHub/BrainNetViewer_20191031'))
BrainNet


%% plot in my own way
addpath(genpath('/home/zhibinz2/Documents/GitHub/eeglab'))
% load corti_source_data
cd /ssd/zhibin/archive/EEG_stroke_62_cleaned20240423
figure;
% for p=1:60; % 12
% p=52; % remove bad channel from leadfield matrix
p=12
% p=50 
% p= 5 
% p= 35
% p = 48 
% p =1
load([num2str(p) '.mat']);
% reduced the size
wind=[91000:96000];
preprocessed_eeg=preprocessed_eeg(:,wind);
preprocessed_eeg = zscore(preprocessed_eeg, 0, 2); % Normalize across each channel
% plot(wind,preprocessed_eeg);
% preprocessed_eeg=preprocessed_eeg;

% Define delta band range
freq_band = [30 50]; % Hz

% Sampling rate
fs = 1000; % Hz

% Compute power spectral density (PSD) using Welch's method
[pxx, f] = pwelch(preprocessed_eeg', [], [], [], fs); % Transpose for pwelch
% Find indices corresponding to delta band
freq_idx = (f >= freq_band(1)) & (f <= freq_band(2));
% Sum power in delta band for each channel
freq_power = sum(pxx(freq_idx, :), 1); % Sum across frequency range
% Plot topography using EEGLAB's topoplot

figure;
% subplot(6,10,p)
% clf
topoplot(freq_power, chanlocs,'electrodes', 'on');colorbar;colormap('jet');clim([min(freq_power) max(freq_power)]);title('power')
% topoplot(freq_power, chanlocs, 'electrodes', 'labelpoint');
% topoplot(std(preprocessed_eeg'), chanlocs, 'electrodes', 'labelpoint');
% colorbar;
% title(num2str(p))
% % title('Delta Band Power (0.5-4 Hz)');
% pause(0.25);
% end
% sgtitle('zscore power freq 30-50 wind 91000-96000')

figure;
clf
% topoplot(freq_power, chanlocs,'electrodes', 'on');
topoplot(std(preprocessed_eeg'), chanlocs, 'electrodes', 'labelpoint');
colorbar;
title(num2str(p))
% title('Delta Band Power (0.5-4 Hz)');
% end

%leadfield remove badchans
badchans=[];badchansplus=[];
badchans=unique([badchans badchansplus'])
cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/base_files/MNE/EGI256
load('MNE_source_model.mat')
leadfield(badchans,:)=zeros;

% source
addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
[inversemat] = inversemodel(leadfield,'prctile',75);
source_data=inversemat*preprocessed_eeg;
% figure;plotx(std(source_data'))

% correlation
reconEEG=leadfield*source_data;
corrmat=corr(reconEEG',preprocessed_eeg');
% imagesc(corrmat);colorbar
% figure;clf;topoplot(diag(corrmat), chanlocs,'electrodes', 'on');colorbar
figure; topoplot(diag(corrmat), chanlocs,'electrodes', 'labelpoint');colorbar
plot(1:size(preprocessed_eeg,1),diag(corrmat))
badchansplus=find(diag(corrmat)<0.9)
plot(1:length(wind),preprocessed_eeg);

% Compute power spectral density (PSD) using Welch's method
[pxx, f] = pwelch(reconEEG', [], [], [], fs); % Transpose for pwelch
% Find indices corresponding to delta band
freq_idx = (f >= freq_band(1)) & (f <= freq_band(2));
% Sum power in delta band for each channel
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


% Compute power spectral density (PSD) using Welch's method
[pxx, f] = pwelch(source_data', [], [], [], fs); % Transpose for pwelch
% Find indices corresponding to delta band
freq_idx = (f >= freq_band(1)) & (f <= freq_band(2));
% Sum power in delta band for each channel
freq_power = sum(pxx(freq_idx, :), 1); % Sum across frequency range

% cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/base_files/MNE/EGI256
% source_rrX=source_rr(:,1);source_rrY=source_rr(:,2);source_rrZ=source_rr(:,3);
cd /home/zhibinz2/Documents/GitHub/STROKE_P61/lesion_mask_on_L_p61_20240427/networkx/new_simu
figure
values_display=freq_power;
% values_display=std(source_data');
s=scatter3(source_rrX, source_rrY, source_rrZ, pointsize, values_display, ...
        "filled", 'MarkerFaceAlpha',0.9);
xlabel('x');ylabel('y');zlabel('z');
view(0,90); 
[caz,cel] = view;
view(caz,cel)
colormap("cool")
colorbar; 



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


cd /home/zhibinz2/Documents/GitHub/STROKE_P61/lesion_mask_on_L_p61_20240427/networkx/new_simu
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
% clim([0 2])
% title(titles{nx})
% subtitle([bandlabels{freq}])
axis off

%
