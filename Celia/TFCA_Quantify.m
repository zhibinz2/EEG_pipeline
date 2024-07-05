function [N,SURFACE,MEANincl0,MEANexcl0,VOLUME] = TFCA_Quantify(srate,freq,SRoWCS,WMSC,woi_bnds,f_bnds)
% -------------------------------------------------------------------------
% TFCA_Quantify, last modified 26-03-2019
% QUANTIFICATION of time-frequency wavelet-based magnitude-squared coherence
% -----
% David AMARANTINI, PhD    <(-_-)>
%	david.amarantini@inserm.fr / david.amarantini@univ-tlse3.fr
%	Paul Sabatier University (UPS)
%	TOulouse NeuroImaging Center (ToNIC, UMR 1214 Inserm / UPS)
% -------------------------------------------------------------------------
% Syntax:
% [N,SURFACE,MEANincl0,MEANexcl0,VOLUME] = TFCA_Quantify(srate,freq,SRoWCS,WMSC,woi_bnds,f_bnds) ;
% -----
% Required input argument list:
% - srate: Sampling frequency
% - freq: Frequency vector (Hz)
% - SRoWCS: Significant areas of correlation in the time-frequency plane
% - WMSC: Time-frequency wavelet-based magnitude-squared coherence
% -----
% Optional input argument list:
% - woi_bnds (2x1 matrix): Lower and upper bounds of the time window of interest in second (s)
% - f_bnds (2x1 matrix): Lower and upper bounds of the frequency band of interest in Hertz (Hz)
% -----
% Output argument list:
% - N: Number of time-frequency points where correlation is significant in the time-frequency plane
% - SURFACE: Surface of the areas where correlation is significant in the time-frequency plane
% - MEANincl0: Mean including 0 of magnitude-squared coherence values where correlation is significant inside the limits of interest
% - MEANexcl0: Mean excluding 0 of magnitude-squared coherence values where correlation is significant inside the limits of interest
% - VOLUME: Volume under magnitude-squared coherence values where correlation is significant inside the limits of interest
% -------------------------------------------------------------------------

disp('<(-_-)> Let''s go for wavelet magnitude-squared coherence quantification!') ;

%% ------------------------------------------------------------------------
%% QUANTIFICATION: Quantify time-frequency magnitude-squared coherence on a time window of interest and in a frequency band
% -------------------------------------------------------------------------
% freq bin size = 1/sampling range
dfreq = mean(diff(freq)) ;
% -------------------------------------------------------------------------
% period
dt = 1/srate ;
% -------------------------------------------------------------------------
% Time window of interest
switch nargin ;
    case {5 6} ;
        if ~isempty(woi_bnds) ;
            woi_lb = woi_bnds(1) ;
            woi_ub = woi_bnds(2) ;
        else woi_lb = 1./srate ;
             woi_ub = size(WMSC,2)./srate ;
        end
end
woi_band = woi_lb:dt:woi_ub ;
% -----
n_woi_lb = round(woi_lb.*srate) ; % woi_lb expressed in number of observations
% n_woi_ub = round(woi_ub.*srate) ; % woi_ub expressed in number of observations
% n_band = n_woi_lb:1:n_woi_ub ; % Observation band of interest
n_band = n_woi_lb:1:n_woi_lb+(length(woi_band)-1) ; % Observation band of interest
% -------------------------------------------------------------------------
% Frequency band of interest
switch nargin ;
    case {6} ;
        if ~isempty(f_bnds) ;
            f_lb = f_bnds(1) ;
            f_ub = f_bnds(2) ;
        else f_lb = freq(1) ;         
             f_ub = freq(end) ;
        end
end
% -----
f_band = find(freq<=f_ub & freq>=f_lb) ; % Observation numbers of the frequency band of interest
% -------------------------------------------------------------------------
% Wavelet magnitude-squared coherence values where correlation is significant in the time-frequency plane
WMSCwSWCPS = zeros(size(WMSC)) ;
WMSCwSWCPS(find(SRoWCS==1)) = WMSC(find(SRoWCS==1)) ;
% -----
ind_SRoWCS = find(WMSCwSWCPS(f_band,n_band)) ; % Matrix indices where correlation is significant in the time-frequency plane inside the limits of interest (i.e, the time window of interest within the frequency bounds of interest)
% -------------------------------------------------------------------------
% Dependent variables
N = length(ind_SRoWCS) ; % Number of time-frequency points where correlation is significant in the time-frequency plane
switch N ;
    case 0 ;
        SURFACE = 0 ; % Surface of the areas where correlation is significant in the time-frequency plane
        MEANincl0 = 0 ; % Mean including 0 of magnitude-squared coherence values where correlation is significant inside the limits of interest
        MEANexcl0 = 0 ; % Mean excluding 0 of magnitude-squared coherence values where correlation is significant inside the limits of interest
        VOLUME = 0 ; % Volume under magnitude-squared coherence values where correlation is significant inside the limits of interest
    otherwise
        % Wavelet magnitude-squared coherence values where correlation is significant in the time-frequency plane inside the limits of interest (i.e, the time window of interest within the frequency bounds of interest)
        iloi_WMSCwSWCPS = WMSCwSWCPS(f_band,n_band) ;
        % -------------------------------------------------------------------------
        SURFACE = -N*dfreq*dt ; % Surface of the areas where correlation is significant in the time-frequency plane
        MEANincl0 = mean(mean(iloi_WMSCwSWCPS)) ; % Mean including 0 of magnitude-squared coherence values where correlation is significant inside the limits of interest
        MEANexcl0 = mean(iloi_WMSCwSWCPS(ind_SRoWCS)) ; % Mean excluding 0 of magnitude-squared coherence values where correlation is significant inside the limits of interest
        VOLUME = -trapz(woi_band,trapz(freq(f_band),iloi_WMSCwSWCPS,1),2) ; % Volume under magnitude-squared coherence values where correlation is significant inside the limits of interest
end
%% ------------------------------------------------------------------------
%% END MESSAGE
fprintf('=> Wavelet magnitude-squared coherence quantification completed! <(-_-)>\n\n') ;
%% ------------------------------------------------------------------------