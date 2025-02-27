clear

eeglab
close all

DOSSIERS={'E:\DATASET\4_Source_reconstructed\Control','E:\DATASET\4_Source_reconstructed\patient'};
DOSSIERS_Save={'D:\DATASET\4_Power_Calculation'};

for d=1:1:length(DOSSIERS)
    dossier1=DOSSIERS{1,d};
F=dir(fullfile(dossier1,'*.mat'));

for f=1:1:length(F)
close all
fichier =  F(f).name;
load(fullfile(F(f).folder,fichier));

% Fourier calcul
ch_peripheral={'241', '242', '243', '238', '239', '240', ...
    '244', '245', '246', '247', '251', '256', '91', '102', '111', '120', '133', '145', '165', '174', '187', '199', '208', '216', '229', '233', '237', '236', '235', '234', ...
    '232', '228', '217', '209', '200', '188', '175', '166', '156', '146', '134', '121', '112', '103', '92', '82', '255', '250'};


    [spectra,freqs,speccomp,contrib,specstd]=spectopo(preprocessed_eeg, 0,1000,'freqrange',[0.5 30]);

    ch_peripheral_numeric = str2double(ch_peripheral);
% Ajouter les NaN

% Trier les indices pour éviter les problèmes lors de l'insertion
ch_peripheral_sorted = sort(ch_peripheral_numeric);

% Créer une matrice temporaire avec une taille augmentée
% Create a temporary matrix with an increased size
newNumRows = size(spectra, 1) + length(ch_peripheral_sorted);
newSpectra = NaN(newNumRows, size(spectra, 2));
X=1
a=0
for x=1:1:size(spectra, 1) + length(ch_peripheral_sorted)
    if ismember(x, ch_peripheral_sorted);
    else
        a=a+1;
        newSpectra(x,:)=spectra(a,:);
    end
end

spectra=newSpectra;

disp(strcat('************************ Traitement*****',fichier));

IdxDelta = [1:3]';
IdxTheta = [4:7]';
IdxAlpha = [8:12]';
IdxBeta = [13:30]';
IdxLowBeta=[13:19]';
IdxHighBeta=[20:30]';
IdxGamma = [31:45]';
IdxAll= [1:45]';

% Mean in each frequency band
Power_delta=mean(spectra(:,IdxDelta),2);
Power_theta=mean(spectra(:,IdxTheta),2);
Power_alpha=mean(spectra(:,IdxAlpha),2);
Power_low_beta=mean(spectra(:,IdxLowBeta),2);
Power_high_beta=mean(spectra(:,IdxHighBeta),2);
Power_beta=mean(spectra(:,IdxBeta),2);
Power_gamma=mean(spectra(:,IdxGamma),2);
Power_all=mean(spectra(:,IdxAll),2);

% Sum in each frequency band
spectra=10.^(spectra/10); %power brut

Power_delta_sum=sum(spectra(:,IdxDelta),2);
Power_theta_sum=sum(spectra(:,IdxTheta),2);
Power_alpha_sum=sum(spectra(:,IdxAlpha),2);
Power_beta_sum=sum(spectra(:,IdxBeta),2);
Power_low_beta_sum=sum(spectra(:,IdxLowBeta),2);
Power_high_beta_sum=sum(spectra(:,IdxHighBeta),2);
Power_gamma_sum=sum(spectra(:,IdxGamma),2);
Power_all_sum=sum(spectra(:,IdxAll),2);


% Storage
%Mean
clear POWER
POWER.DELTA.Mean(:,f)=Power_delta;
POWER.THETA.Mean(:,f)=Power_theta;
POWER.ALPHA.Mean(:,f)=Power_alpha;
POWER.BETA.Mean(:,f)=Power_beta;
POWER.LOWBETA.Mean(:,f)=Power_low_beta;
POWER.HIGHBETA.Mean(:,f)=Power_high_beta;
POWER.GAMMA.Mean(:,f)=Power_gamma;
POWER.ALL.Mean(:,f)=Power_all;

%MeanBrut
Power_delta=mean(spectra(:,IdxDelta),2);
Power_theta=mean(spectra(:,IdxTheta),2);
Power_alpha=mean(spectra(:,IdxAlpha),2);
Power_beta=mean(spectra(:,IdxBeta),2);
Power_low_beta=mean(spectra(:,IdxLowBeta),2);
Power_high_beta=mean(spectra(:,IdxHighBeta),2);
Power_gamma=mean(spectra(:,IdxGamma),2);
Power_all=mean(spectra(:,IdxAll),2);

POWER.DELTA.MeanB(:,f)=Power_delta;
POWER.THETA.MeanB(:,f)=Power_theta;
POWER.ALPHA.MeanB(:,f)=Power_alpha;
POWER.BETA.MeanB(:,f)=Power_beta;
POWER.LOWBETA.MeanB(:,f)=Power_low_beta;
POWER.HIGHBETA.MeanB(:,f)=Power_high_beta;
POWER.GAMMA.MeanB(:,f)=Power_gamma;
POWER.ALL.MeanB(:,f)=Power_all;

%Normalisation of sum within frequency band by total sum
POWER.DELTA.Sum(:,f)=Power_delta_sum*100./Power_all_sum;
POWER.THETA.Sum(:,f)=Power_theta_sum*100./Power_all_sum;
POWER.ALPHA.Sum(:,f)=Power_alpha_sum*100./Power_all_sum;
POWER.BETA.Sum(:,f)=Power_beta_sum*100./Power_all_sum;
POWER.LOWBETA.Sum(:,f)=Power_low_beta_sum*100./Power_all_sum;
POWER.HIGHBETA.Sum(:,f)=Power_high_beta_sum*100./Power_all_sum;
POWER.GAMMA.Sum(:,f)=Power_gamma_sum*100./Power_all_sum;
POWER.ALL.Sum(:,f)=Power_all_sum;


POWER.SPECTRA(:,:,f)=spectra;
POWER.INFO.Name{1,f}=fichier;
end

% Save  Data
SaveFolder=DOSSIERS_Save{1,1};
cd   (SaveFolder)


% Retirer ces éléments de X
POWER.INFO.Label = ch_labels;

 if d==1
    save(['Control.mat'],'POWER') 
 elseif d==2
    save(['Patient.mat'],'POWER') 
 end

clearvars -except d f DOSSIERS_Save DOSSIERS
end




%% Representation

load('E:\DATASET\4_Power_Calculation\Without_source_reconstruction\Patient_256_Power.mat')

cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/base_files/MNE/EGI256
load('chanlocs.mat');


A=mean(POWER.BETA.Sum,2)
topoplot(A,chanlocs)%,'maplimits',[-45 20])
colorbar

