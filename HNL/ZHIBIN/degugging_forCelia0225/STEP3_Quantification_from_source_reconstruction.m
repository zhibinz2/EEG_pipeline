clear

eeglab
close all

DOSSIERS={'E:\DATASET\4_Source_reconstructed\Control','E:\DATASET\4_Source_reconstructed\patient'};
DOSSIERS_Save={'E:\DATASET\4_Power_Calculation'};

for d=1:1:length(DOSSIERS)
    dossier1=DOSSIERS{1,d};
F=dir(fullfile(dossier1,'*.mat'));

for f=1:1:length(F)
close all
fichier =  F(f).name;
load(fullfile(F(f).folder,fichier));

% Fourier calcul

    [spectra,freqs,speccomp,contrib,specstd]=spectopo(corti_source_data', 0,1000,'freqrange',[0.5 30]);



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
POWER.INFO.Label = corti_ave_source_labl;

 if d==1
    save(['Control.mat'],'POWER') 
 elseif d==2
    save(['Patient.mat'],'POWER') 
 end

clearvars -except d f DOSSIERS_Save DOSSIERS
end




