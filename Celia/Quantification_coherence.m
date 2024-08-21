<<<<<<< HEAD
clear
DOSSIERS={'E:\DATASET\4_Source_reconstructed\Control','E:\DATASET\4_Source_reconstructed\patient'};
dossiers={'Control','Patient'};
% HNL Method
Fs=500;
srnew = 200;
downsample = Fs/srnew;

passbands = [1 4;4.5 8;8.5 12; 12.5 30; 30.5 50];
bandlabels = {'Delta','Theta', 'Alpha', 'Beta', 'Gamma'};

load('E:\DATASET\Script\Individual_Connectome.mat');

% make filter
attenuation=60;
filt_ds=cell(1,length(passbands));
for freqBand=1:length(bandlabels)
    % Select frequency
    passFreq1 = passbands(freqBand,1);
    passFreq2 = passbands(freqBand,2);
    filt_d = designfilt('bandpassiir','FilterOrder',20, ...
        'PassbandFrequency1',passFreq1,'PassbandFrequency2',passFreq2, ...
        'StopbandAttenuation1',attenuation,'PassbandRipple',0.2, ...
        'StopbandAttenuation2',attenuation,'SampleRate',srnew);
    filt_ds{freqBand}=filt_d;
end

for d=1:1:length(DOSSIERS)

    dossier1=DOSSIERS{1,d};
    F=dir(fullfile(dossier1,'*.mat'));

    for f=1:1:length(F)

        fichier =  F(f).name;

        disp(fichier);
        load(fullfile(F(f).folder,fichier));

        % select one trial
        example_data=corti_source_data;

        for Hz=1:1:length(bandlabels)
            % pick a frequency to process the data

            filtered_data = filter(filt_ds{Hz},example_data);
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

            % set a threshold of 0.4 and convert to binary logical values
            stroke_coh_binary=zeros(size(stroke_coh));
            stroke_coh_binary(stroke_coh>0.4)=1;

            figure;
            subplot(121)
            imagesc(stroke_coh);colorbar
            title('coherence')
            subplot(122)
            imagesc(stroke_coh_binary);colorbar
            sgtitle('Binary matricies for network analysis')

            Coh.(dossiers{d}).(bandlabels{Hz}).Data(:,:,f)=stroke_coh;


            % signficativité de la cohérence

            % AGL variables
            % https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh
            addpath(genpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh'))
            allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
            allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
            n_Lambdas=length(allLambdas); % number of lambda values
            min_LamdaIn = min(allLambdas);

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
            toc


            tic
            [stroke_Pcov] = fitprecision( ...
                Individual_Connectome,penalizationIn_op,penalizationOut_op,min_LamdaIn,stroke_Cov);
            toc

            % add this repo to your path (https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh)
            addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util

            % convert real value covariance to complex matrix (448x448)
            % then compute coherence
            stroke_coh=normalizeCSD(r2c(stroke_Cov)); % ordinary coherence
            stroke_Pcoh=normalizeCSD(r2c(logical(stroke_Pcov)));% partial coherence using boolean

            %%%%%%%%%%%%%%%%%%%%SAVE structure and clear
        end
    end
end

save('Coherence.mat','Coh')



=======
clear
DOSSIERS={'E:\DATASET\4_Source_reconstructed\Control','E:\DATASET\4_Source_reconstructed\patient'};
dossiers={'Control','Patient'};
% HNL Method
Fs=500;
srnew = 200;
downsample = Fs/srnew;

passbands = [1 4;4.5 8;8.5 12; 12.5 30; 30.5 50];
bandlabels = {'Delta','Theta', 'Alpha', 'Beta', 'Gamma'};

load('E:\DATASET\Script\Individual_Connectome.mat');

% make filter
attenuation=60;
filt_ds=cell(1,length(passbands));
for freqBand=1:length(bandlabels)
    % Select frequency
    passFreq1 = passbands(freqBand,1);
    passFreq2 = passbands(freqBand,2);
    filt_d = designfilt('bandpassiir','FilterOrder',20, ...
        'PassbandFrequency1',passFreq1,'PassbandFrequency2',passFreq2, ...
        'StopbandAttenuation1',attenuation,'PassbandRipple',0.2, ...
        'StopbandAttenuation2',attenuation,'SampleRate',srnew);
    filt_ds{freqBand}=filt_d;
end

for d=1:1:length(DOSSIERS)

    dossier1=DOSSIERS{1,d};
    F=dir(fullfile(dossier1,'*.mat'));

    for f=1:1:length(F)

        fichier =  F(f).name;

        disp(fichier);
        load(fullfile(F(f).folder,fichier));

        % select one trial
        example_data=corti_source_data;

        for Hz=1:1:length(bandlabels)
            % pick a frequency to process the data

            filtered_data = filter(filt_ds{Hz},example_data);
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

            % set a threshold of 0.4 and convert to binary logical values
            stroke_coh_binary=zeros(size(stroke_coh));
            stroke_coh_binary(stroke_coh>0.4)=1;

            figure;
            subplot(121)
            imagesc(stroke_coh);colorbar
            title('coherence')
            subplot(122)
            imagesc(stroke_coh_binary);colorbar
            sgtitle('Binary matricies for network analysis')

            Coh.(dossiers{d}).(bandlabels{Hz}).Data(:,:,f)=stroke_coh;


            % signficativité de la cohérence

            % AGL variables
            % https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh
            addpath(genpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh'))
            allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
            allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);
            n_Lambdas=length(allLambdas); % number of lambda values
            min_LamdaIn = min(allLambdas);

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
            toc


            tic
            [stroke_Pcov] = fitprecision( ...
                Individual_Connectome,penalizationIn_op,penalizationOut_op,min_LamdaIn,stroke_Cov);
            toc

            % add this repo to your path (https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh)
            addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util

            % convert real value covariance to complex matrix (448x448)
            % then compute coherence
            stroke_coh=normalizeCSD(r2c(stroke_Cov)); % ordinary coherence
            stroke_Pcoh=normalizeCSD(r2c(logical(stroke_Pcov)));% partial coherence using boolean

            %%%%%%%%%%%%%%%%%%%%SAVE structure and clear
        end
    end
end

save('Coherence.mat','Coh')



>>>>>>> ec17b3a49212a2b46099454f24f2ebf03c929aeb
