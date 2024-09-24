% fullConnectome463=logical(squeeze(fullConnectome_p61(1,:,:)));
% fullConnectome448=fullConnectome463(corti_ave_source_labl,corti_ave_source_labl);
% roi448names;
% roi463names=ROI463_texts;

% save("extract_fullconnectome.mat",'fullConnectome463','fullConnectome448',"roi448names",'roi463names');

%% load files
load('extract_fullconnectome.mat')

figure;
subplot(121)
imagesc(fullConnectome463);title('463 ROI full connectome');
subplot(122)
imagesc(fullConnectome448);title('448 ROI full connectome');


% In 'extract_fullconnectome.mat', you will find:
% fullConnectome463 - full intact binary structure connectome of all 463 ROIs
% fullConnectome448 - full intact binary structure connectome of 448 cortical ROIs
% roi463names - names for all 463 ROIs
% roi448names - names for 448 cortical ROIs
