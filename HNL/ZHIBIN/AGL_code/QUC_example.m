clear
SC=SC_p;
allLambdas;
allLambdasOut;
dataCovs_op;

min_LamdaIn = min(allLambdas);
n_Lambdas=length(allLambdas); % number of lambda values

n_ensam=size(dataCovs_op,1); % number of ensambles
GforFit =[double(SC),double(SC) ; double(SC), double(SC)]; % boolean


% initiate the deviance
allDevs = zeros(length(allLambdas),length(allLambdasOut),n_ensam,n_ensam);

mins = 1; % select one ensemble   % loop through ensembles for cross validation
dataCov = squeeze(dataCovs_op(mins,:,:));
scalingVal = max(triu(abs(dataCov),1),[],'all'); % maximum value of the cov


lambda = 5; % select one lambdas In

tmpDev = zeros(length(allLambdasOut),n_ensam); % 13x4 deviance

lambdaOut = 1 % select one lambdas Out.

current_lamdaOut=allLambdasOut(lambdaOut);
current_lamdaIn=allLambdas(lambda);

% if current_lamdaOut < current_lamdaIn % if Out < In, bypass this loop. We want Out > In
%     continue
% end


In_network = current_lamdaIn*scalingVal*GforFit;

Out_network = current_lamdaOut*scalingVal*(double(~(GforFit))-eye(length(GforFit)));

diag_network = min_LamdaIn*scalingVal*eye(length(GforFit));

addpath(genpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh'))


mode='default';
S=dataCov;
L=In_network + Out_network + diag_network;
tol=1e-4;
msg=0;
maxIter=200;

tic
X = QUIC(mode, S, L, tol, msg, maxIter);
toc



clear
addpath(genpath('..GitHub/AdaptiveGraphicalLassoforParCoh')); % https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh
load('QUC_example_data.mat');
tic
[Xoutput, W, opt, time, iter, dGap] = QUIC('default', S, L, 1e-4, 0, 200);
toc % it takes about 138 seconds
% examining the inputs and output
figure;
subplot(141);imagesc(logical(X));title('my output X');colorbar;colormap('jet');
subplot(142);imagesc(S);title('S');colorbar;colormap('jet');clim([-0.05 0.05])
subplot(143);imagesc(L);title('L');colorbar;colormap('jet');clim([0 20])
subplot(144);imagesc(logical(outputX));title('your output X');colorbar;colormap('jet');

imagesc(logical(Xoutput));title('your output X');colorbar;colormap('jet');
imagesc(W);title('your output W');colorbar;colormap('jet');

save('QUC_example_data.mat','X','mode','S','L');

save('QUC_all_outputs.mat','Xoutput', 'W', 'opt', 'time', 'iter', 'dGap')