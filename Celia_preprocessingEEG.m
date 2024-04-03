%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Routine 1 : prétraitement

% Required :
% EEGLAB v2023.0 with plugins "Biosig" v3.8.1, "ICLabel" v1.4,
% "clean_rawdata" v2.8 and "zapline-plus" v1.2.1
% Other functions : load_acq, butter_filtfilt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization

clear
close all
clc
addpath(genpath('D:\Recherche\eVAICI'));

dirpath = uigetdir('','Select folder with all subject data') ;
list_sub = dir([dirpath '\S*']) ;
fileSel = listdlg('PromptString','Please select subject:','ListString',{list_sub.name},'ListSize',[200 200], 'Name', 'Subject selection');
dirpathsave = uigetdir('','Select folder to save preprocessed data EEG');

EEG_sel_fields = {'nbchan' 'trials' 'pnts' 'srate' 'xmin' 'xmax' 'times' 'data' 'icachansind' 'chanlocs' 'chaninfo' 'ref' 'event' 'etc'};

% ICLabel threshold for Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other
IC_threshold = [NaN NaN;0.9 1;0.9 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN];

%% Subject loop
for s = 1 : numel(fileSel)
    fprintf('-----------------------------< PREPROCESSING OF SUBJECT %s >-----------------------------\n',...
        list_sub(fileSel(s)).name(2:end));

    F_EEG = dir([dirpath '\' list_sub(fileSel(s)).name '\bdf\*.bdf']);


    for file = 1 : numel(F_EEG) % For each EEG file (f) of each subject (s)

        filename =  F_EEG(file).name;
        num = filename(1:3);
        cond = filename (5:7);
        vit = filename(9:11);
        repet = filename(13);

        disp(['%%%%%%%%%%%%%% - CONDITION => ' cond ' // SPEED => ' vit ' // REPETITION => ' repet ' - %%%%%%%%%%%%%%' ])

        %% Load EEG data

        EEG = pop_biosig([F_EEG(file).folder '\' num '_' cond '_' vit '_' repet '.bdf']);


        %% Preprocess EEG data
        EEG.data = double(EEG.data);
        EEG.data(1:64,:) = butter_filtfilt(EEG.data(1:64,:),EEG.srate,2,1,'high');
        %EEG = pop_resample(EEG, rateint); % resampling  between 250 Hz and 500 Hz is recommended to use Zapline-plus

        [EEG.data(1:64,:), EEG.etc.zapline.config, EEG.etc.zapline.results, ~] = clean_data_with_zapline_plus(EEG.data(1:64,:),...
            EEG.srate,'noisefreqs',[50 100],'minfreq',48,'maxfreq',102, 'plotResults',0); % remove line noise

        EEG_remove = pop_clean_rawdata(EEG_remove, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off',...
            'BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');

        % A visualization step could be necessary here to reject potential time window, in particular at the begening and the ending of the recording

        %%% Run ICA, label and reject components with ICLabel %%%
        EEG_remove = pop_runica(EEG_remove, 'icatype', 'runica', 'extended',1,'interrupt','on');
        EEG_remove = pop_iclabel(EEG_remove, 'default');

        EEG_remove = pop_icflag(EEG_remove, IC_threshold);
        EEG_clean = pop_subcomp(EEG_remove, [],0);
        EEG_clean.etc.ic_classification.ICLabel = EEG_remove.etc.ic_classification.ICLabel;
        EEG_clean.etc.ic_classification.ICLabel.ic_removed = EEG_clean.etc.ic_classification.ICLabel.classifications  >= IC_threshold(:,1)';
        EEG_clean.etc.ic_classification.ICLabel.ic_removed_nb = sum(EEG_clean.etc.ic_classification.ICLabel.ic_removed);

        % Interpolate bad channels
        EEG_clean = eeg_interp (EEG_clean, EEG.chanlocs(1:64), 'spherical');

        % Average re-referencing
        EEG_clean = pop_reref(EEG_clean, []); %There is no consensus in the literature to do the re-referencing before or after the ICA...

        
        % Making epoch if necessary, accordingly to the triggers
        % Storage in a "DATA" structure


        clearvars -except rateint dirpath dirpathsave list_sub fileSel EEG_sel_fields IC_threshold s F_* num file DATA
    end

    % Save preprocessed data
    save(fullfile(dirpathsave,strcat(num, '_PREPROCESSED.mat')),'DATA','-v7.3');

    clear DATA F_* num
end
