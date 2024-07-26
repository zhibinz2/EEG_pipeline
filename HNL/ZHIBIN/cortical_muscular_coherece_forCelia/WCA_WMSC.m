function [WMSC,t,freq] = WCA_WMSC(s1,s2,rate)
% WCA_WCS_WMSC Summary of this function goes here

% Inputs:
% - s1, s2: two signals ( N trials x n observations )
% - rate: sampling frequency

% Outputs:
% - WMSC: Wavelet magnitude-squared coherence

% ----- Echantillonnage des signaux ------------------------------------------------------

n = size(s1,2) ; % Number of observation
N = size(s1,1) ; % Number of trials

dt = 1/rate ; % Période d'échantillonnage (sampling period)
T = n / rate ; % Durée totale en seconde (total duration in seconds)
t = linspace(0,T,rate*T) ; % Time vector

% ----- Paramètres WaveCrossSpec (WaveCrossSpec parameters)
nvoice = 5;
J1 = 5;
wavenumber = 10;

% alpha = 0.05 ;

% ----------------------------------------------------------------------------------------

% ----- Paramétrisation WaveCrossSpec (Bigot et al., 2011) -------------------------------
% ----- WaveCrossSpec parameterization
Args=struct('Pad',1,...
    'DJ',nvoice, ...
    'S0',2*dt,...
    'J1',J1,...
    'Mother','Morlet',...
    'Cycles',wavenumber) ;

% ----------------------------------------------------------------------------------------

% ----- Données -------------------------------------------------------
% ----- Data 
S1 = mean(s1,1) ; % Signal moyen de l'electrode "1"
S2 = mean(s2,1) ; % Signal moyen de l'electrode "2"

% ----- Time-Frequency Coherence Analysis via Wavelet Transform --------------------------
% ----------
% Wavelet transforms (need wavelet function)
% Found one in https://github.com/grinsted/wavelet-coherence
% https://github.com/grinsted/wavelet-coherence/blob/master/private/wavelet.m
% cd /home/zhibinz2/Documents/GitHub/wavelet-coherence/private
cd C:\Users\zhouz\GitHub\wavelet-coherence\private
% cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/Celia
for i = 1:N
    [WT_S1{i},period,scale,coi] = wavelet(s1(i,:)',dt,Args.Pad,Args.DJ,Args.S0,Args.J1,Args.Mother,Args.Cycles) ;
    [WT_S2{i},period,scale,coi] = wavelet(s2(i,:)',dt,Args.Pad,Args.DJ,Args.S0,Args.J1,Args.Mother,Args.Cycles) ;
end

freq = 1./period ;
nb_scale = length(freq) ;

% ----------
% Wavelet auto spectrum / wavelet power spectrum / Scalogram
WS_S1 = zeros(nb_scale,n) ;
WS_S2 = zeros(nb_scale,n) ;
for i = 1:N ;
    WS_S1 = WS_S1+abs(WT_S1{i}).^2 ; % Scalogram
    WS_S2 = WS_S2+abs(WT_S2{i}).^2 ; % Scalogram
end
WPS_S1 = WS_S1/N ; % Wavelet power spectrum
WPS_S2 = WS_S2/N ; % Wavelet power spectrum


% ----------
% Wavelet cross spectrum
WCS = zeros(nb_scale,n) ;
for i = 1:N ;
    WCS = WCS + (WT_S1{i}).*conj(WT_S2{i}) ;
end
% ----------
% Wavelet magnitude-squared coherence
WMSC = zeros(nb_scale,n) ;
WMSC = (abs(WCS).^2)./(WS_S1.*WS_S2) ;  % Wavelet magnitude-squared coherence

% Plot
% imagesc(WMSC) ;
end


