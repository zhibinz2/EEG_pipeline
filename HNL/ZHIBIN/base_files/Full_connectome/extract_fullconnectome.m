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

