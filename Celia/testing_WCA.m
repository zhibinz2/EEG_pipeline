% ----------------------------------------------------------------------------------------
% Time-Frequency Coherence Analysis via Wavelet Transform
% ----------
% David AMARANTINI, PhD
% david.amarantini@inserm.fr / david.amarantini@univ-tlse3.fr
% ToNIC, Toulouse NeuroImaging Center, Universit� de Toulouse, Inserm, UPS, France
% ----------
% Outputs
% - S1: Mean signal "1"
% - S2: Mean signal "2"
% - WPS_S1: Wavelet power spectrum (scalogram) of signal "1"
% - WPS_S2: Wavelet power spectrum (scalogram) of signal "2"
% - WCPS: Wavelet cross power spectrum
% - WPD: Wavelet phase difference
% - WMSC: Wavelet magnitude-squared coherence
% - SRoWCS: Significant regions of the wavelet cross-spectrum
% - SRoWC: Significant regions of the wavelet coherence
% - freq: Vector of frequencies
% ----------------------------------------------------------------------------------------


% load my data example
cd /home/zhibinz2/Documents/GitHub/EEG_pipeline/HNL/ZHIBIN/GraphicalLasso_example
load('trialdata.mat')


% ----- Echantillonnage des signaux ------------------------------------------------------
s1 = trialdata(:,:,1); % Signal 1 centr� : N trials x n observations (use channel 1)
s2 = trialdata(:,:,3);; % Signal 2 centr�: N trials x n observations (use channel 1)
rate =  2000; % Frequency (my sampling frequency)

n = size(s1,2) ; % Number of observation
N = size(s1,1) ; % Number of trials

dt = 1/rate ; % P�riode d'�chantillonnage (sampling period)
T = n / rate ; % Dur�e totale en seconde (total duration in seconds)
t = linspace(0,T,rate*T) ; % Time vector

% ----- Param�tres WaveCrossSpec (WaveCrossSpec parameters)
nvoice = 5;
J1 = 5;
wavenumber = 10;

alpha = 0.05 ;

% ----------------------------------------------------------------------------------------

% ----- Param�trisation WaveCrossSpec (Bigot et al., 2011) -------------------------------
% ----- WaveCrossSpec parameterization
Args=struct('Pad',1,...
    'DJ',nvoice, ...
    'S0',2*dt,...
    'J1',J1,...
    'Mother','Morlet',...
    'Cycles',wavenumber) ;

% ----------------------------------------------------------------------------------------

% ----- Donn�es -------------------------------------------------------
% ----- Data 
S1 = mean(s1,1) ; % Signal moyen de l'electrode "1"
S2 = mean(s2,1) ; % Signal moyen de l'electrode "2"

% ----------------------------------------------------------------------------------------

% ----- Time-Frequency Coherence Analysis via Wavelet Transform --------------------------
% ----------
% Wavelet transforms (need wavelet function)
% Found one in https://github.com/grinsted/wavelet-coherence
% https://github.com/grinsted/wavelet-coherence/blob/master/private/wavelet.m
cd /home/zhibinz2/Documents/GitHub/wavelet-coherence/private
for i = 1:N
    waitbar(i/N) ; % Progress waitbar
    [WT_S1{i},period,scale,coi] = wavelet(s1(i,:)',dt,Args.Pad,Args.DJ,Args.S0,Args.J1,Args.Mother,Args.Cycles) ;
    [WT_S2{i},period,scale,coi] = wavelet(s2(i,:)',dt,Args.Pad,Args.DJ,Args.S0,Args.J1,Args.Mother,Args.Cycles) ;
    freq = 1./period ;
    nb_scale = length(freq) ;
end


% ----------
% Wavelet auto spectrum / wavelet power spectrum / Scalogram
WS_S1 = zeros(nb_scale,n) ;
WS_S2 = zeros(nb_scale,n) ;
for i = 1:N ;
    waitbar(i/N) ; % Progress waitbar
    WS_S1 = WS_S1+abs(WT_S1{i}).^2 ; % Scalogram
    WS_S2 = WS_S2+abs(WT_S2{i}).^2 ; % Scalogram
end
WPS_S1 = WS_S1/N ; % Wavelet power spectrum
WPS_S2 = WS_S2/N ; % Wavelet power spectrum


% ----------
% Wavelet cross spectrum
WCS = zeros(nb_scale,n) ;
for i = 1:N ;
    waitbar(i/(N+6)) ; % Progress waitbar
    WCS = WCS + (WT_S1{i}).*conj(WT_S2{i}) ;
end
% ----------
% Wavelet cross power spectrum
waitbar((i+1)/(N+6)) ; % Progress waitbar
WCPS = abs(WCS/N) ;
% ----------
% Wavelet squared cross power spectrum
waitbar((i+2)/(N+6)) ; % Progress waitbar
WSCPS = WCPS.^2 ;
% ----------
% Wavelet phase difference in radians
waitbar((i+3)/(N+6)) ; % Progress waitbar
WPD = atan(imag(WCS)./real(WCS)) ;
% ----------
% Wavelet magnitude-squared coherence
waitbar((i+4)/(N+6)) ; % Progress waitbar
WMSC = zeros(nb_scale,n) ;
WMSC = (abs(WCS).^2)./(WS_S1.*WS_S2) ;  % Wavelet magnitude-squared coherence

% ----------
% Test of significativity of the cross spectrum
% Covariance matrices of the data and computation of their largest eigenvalues
waitbar((i+5)/(N+6)) ; % Progress waitbar
beta = 1 ;
indice = 1:n ;
X1 = zeros(N,length(indice)) ;
X2 = zeros(N,length(indice)) ;
for i = 1:N ;
    X1(i,:) = s1(i,indice) ;
    X2(i,:) = s2(i,indice) ;
end
eps1 = sqrt(max(eig(transpose(X1)*X1/N))/((1+sqrt(length(indice)/N)+sqrt(-2*log(beta)/N))^2)) ;
eps2 = sqrt(max(eig(transpose(X2)*X2/N))/((1+sqrt(length(indice)/N)+sqrt(-2*log(beta)/N))^2)) ;
% Significant regions of the wavelet cross power spectrum
SRoWCS = zeros(nb_scale,n) ;
thr = eps1*eps2*(-log(alpha/2)/N+sqrt(-2*log(alpha/2)/N)) ;
for s = 1:nb_scale ;
    SRoWCS(s,:) = WSCPS(s,:) > thr^2 ;
end
% ----------
% Test of significativity of the coherence
% Significant regions of the wavelet magnitude-squared coherence
waitbar((i+6)/(N+6)) ; % Progress waitbar
SRoWC = zeros(nb_scale,n) ;
thr = 1-alpha^(1/(N-1)) ;
for s = 1:nb_scale
    SRoWC(s,:) = WMSC(s,:) > thr ;
end
% ----------------------------------------------------------------------------------------

% ----- Display --------------------------------------------------------------------------
scrsz = get(0,'ScreenSize') ;
% ----------
% Plot summary figure
indiceX = floor(linspace(1,length(t),5)) ;
indiceY = floor(linspace(1,length(freq),5)) ;
h = figure('NumberTitle','off','Name','EEG-EEG Time-frequency Coherence Analysis','Position',[0 0 scrsz(3) scrsz(4)]) ;
% ----------
subplot(7,2,1) ;
plot(t,mean(s1,1)) ;
title('S1 - Mean signal') ;
xlabel('Time (s)')
ylabel('Amplitude')
% ----------
subplot(7,2,2) ;
plot(t,mean(s2,1)) ;
title('S2 - Mean signal') ;
xlabel('Time (s)')
ylabel('Amplitude')
% ----------
subplot(7,2,3) ;
imagesc(WPS_S1) ;
title('S1 - Power spectrum') ; % Plot Wavelet power spectrum
set(gca,'YTickLabel',[]) ;
for p =1:1:length(indiceY) ;
    text(-20,indiceY(p),sprintf('%0.0f',freq(indiceY(p))),'HorizontalAlignment','right') ;
end
set(gca,'XTickLabel',[]) ;
for p =1:1:length(indiceX) ;
    text(indiceX(p),300,sprintf('%0.2f',t(indiceX(p))),'HorizontalAlignment','center') ;
end
xlabel('Time (s)','Position',[length(t)/2,(1.25)*length(freq),0])
ylabel('Frequency (Hz)','Position',[-60,length(freq)/2,0])
% ----------
subplot(7,2,4) ;
imagesc(WPS_S2) ;
title('S2 - Power spectrum') ; % Plot Wavelet power spectrum
set(gca,'YTickLabel',[]) ;
for p =1:1:length(indiceY) ;
    text(-20,indiceY(p),sprintf('%0.0f',freq(indiceY(p))),'HorizontalAlignment','right') ;
end
set(gca,'XTickLabel',[]) ;
for p =1:1:length(indiceX) ;
    text(indiceX(p),300,sprintf('%0.2f',t(indiceX(p))),'HorizontalAlignment','center') ;
end
xlabel('Time (s)','Position',[length(t)/2,(1.25)*length(freq),0])
ylabel('Frequency (Hz)','Position',[-60,length(freq)/2,0])
% ----------
subplot(7,2,[5 6]) ;
imagesc(WCPS) ;
title('S1/S2 - Cross power spectrum') ; % Plot wavelet cross power spectrum
set(gca,'YTickLabel',[]) ;
for p =1:1:length(indiceY) ;
    text(-7,indiceY(p),sprintf('%0.0f',freq(indiceY(p))),'HorizontalAlignment','right') ;
end
set(gca,'XTickLabel',[]) ;
for p =1:1:length(indiceX) ;
    text(indiceX(p),300,sprintf('%0.2f',t(indiceX(p))),'HorizontalAlignment','center') ;
end
xlabel('Time (s)','Position',[length(t)/2,(1.25)*length(freq),0])
ylabel('Frequency (Hz)','Position',[-25,length(freq)/2,0])
% ----------
subplot(7,2,[7 8]) ;
imagesc(SRoWCS) ;
title('S1/S2 - Significant regions of the cross spectrum') ; % Plot significant regions of the cross-spectrum
set(gca,'YTickLabel',[]) ;
for p =1:1:length(indiceY) ;
    text(-7,indiceY(p),sprintf('%0.0f',freq(indiceY(p))),'HorizontalAlignment','right') ;
end
set(gca,'XTickLabel',[]) ;
for p =1:1:length(indiceX) ;
    text(indiceX(p),300,sprintf('%0.2f',t(indiceX(p))),'HorizontalAlignment','center') ;
end
xlabel('Time (s)','Position',[length(t)/2,(1.25)*length(freq),0])
ylabel('Frequency (Hz)','Position',[-25,length(freq)/2,0])
% ----------
subplot(7,2,[9 10]) ;
imagesc(WMSC) ;
title('S1/S2 - Magnitude-squared coherence') ; % Plot wavelet cross power spectrum
set(gca,'YTickLabel',[]) ;
for p =1:1:length(indiceY) ;
    text(-7,indiceY(p),sprintf('%0.0f',freq(indiceY(p))),'HorizontalAlignment','right') ;
end
set(gca,'XTickLabel',[]) ;
for p =1:1:length(indiceX) ;
    text(indiceX(p),300,sprintf('%0.2f',t(indiceX(p))),'HorizontalAlignment','center') ;
end
xlabel('Time (s)','Position',[length(t)/2,(1.25)*length(freq),0])
ylabel('Frequency (Hz)','Position',[-25,length(freq)/2,0])
% ----------
subplot(7,2,[11 12]) ;
imagesc(SRoWC) ;
title('S1/S2 - Significant regions of the coherence') ; % Plot significant regions of the coherence
set(gca,'YTickLabel',[]) ;
for p =1:1:length(indiceY) ;
    text(-7,indiceY(p),sprintf('%0.0f',freq(indiceY(p))),'HorizontalAlignment','right') ;
end
set(gca,'XTickLabel',[]) ;
for p =1:1:length(indiceX) ;
    text(indiceX(p),300,sprintf('%0.2f',t(indiceX(p))),'HorizontalAlignment','center') ;
end
xlabel('Time (s)','Position',[length(t)/2,(1.25)*length(freq),0])
ylabel('Frequency (Hz)','Position',[-25,length(freq)/2,0])
% ----------
subplot(7,2,[13 14]) ;
imagesc(WPD) ;
title('S1/S2 - Phase difference (rad)') ; % Plot wavelet phase difference
set(gca,'YTickLabel',[]) ;
for p =1:1:length(indiceY) ;
    text(-7,indiceY(p),sprintf('%0.0f',freq(indiceY(p))),'HorizontalAlignment','right') ;
end
set(gca,'XTickLabel',[]) ;
for p =1:1:length(indiceX) ;
    text(indiceX(p),300,sprintf('%0.2f',t(indiceX(p))),'HorizontalAlignment','center') ;
end
xlabel('Time (s)','Position',[length(t)/2,(1.25)*length(freq),0])
ylabel('Frequency (Hz)','Position',[-25,length(freq)/2,0])
% ----------------------------------------------------------------------------------------

disp('Done!')
% ----------------------------------------------------------------------------------------