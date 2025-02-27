clear
% M = readmatrix('E:\DATASET\Script\BrainNet\Data\SurfTemplate\BrainMesh_ICBM152.txt');
M = readmatrix('/home/zhibinz2/Documents/GitHub/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152.txt');
M(81925:end,:)=[];
Vertices = M;
load('E:\DATASET\4_Power_Calculation\Control.mat')

load('E:\DATASET\4_Source_reconstructed\Control\ARAT.mat')

X=corti_ave_source_coor(:,1);
Y=corti_ave_source_coor(:,2);
Z=corti_ave_source_coor(:,3);

puissance=mean(POWER.DELTA.Sum,2);


puissance=value_display;
puissance=double(puissance);
interp_values = griddata(X, Y, Z, puissance, Vertices(:,1), Vertices(:,2), Vertices(:,3), 'nearest'); % wrong way to interperlate

figure
% values_display=freq_power;
values_display=interp_values;
s=scatter3(Vertices(:,1), Vertices(:,2), Vertices(:,3), pointsize, values_display, ...
        "filled", 'MarkerFaceAlpha',0.9);
xlabel('x');ylabel('y');zlabel('z');
view(0,90); 
[caz,cel] = view;
view(caz,cel)
colormap("cool")
colorbar; 

% right way to do it
% align
lowDensityCoords = corti_ave_source_coor; 
highDensityCoords = Vertices; 
lowCentroid = mean(lowDensityCoords, 1);
highCentroid = mean(highDensityCoords, 1);
translatedLowDensityCoords = lowDensityCoords - lowCentroid + highCentroid;
scaleFactor = mean(vecnorm(highDensityCoords - highCentroid, 2, 2)) / ...
              mean(vecnorm(lowDensityCoords - lowCentroid, 2, 2));
alignedLowDensityCoords = highCentroid + scaleFactor * (lowDensityCoords - lowCentroid);
figure;
scatter3(lowDensityCoords(:,1), lowDensityCoords(:,2), lowDensityCoords(:,3), 50, 'r', 'filled'); hold on;
scatter3(highDensityCoords(:,1), highDensityCoords(:,2), highDensityCoords(:,3), 10, 'b'); 
legend('Low Density (Original)', 'High Density');
title('Before Alignment');
figure;
scatter3(alignedLowDensityCoords(:,1), alignedLowDensityCoords(:,2), alignedLowDensityCoords(:,3), 50, 'r', 'filled'); hold on;
scatter3(highDensityCoords(:,1), highDensityCoords(:,2), highDensityCoords(:,3), 10, 'b'); 
legend('Aligned Low Density', 'High Density');
title('After Alignment');
%Find the nearest aligned low-density point for each high-density point
[idx, ~] = knnsearch(alignedLowDensityCoords, highDensityCoords);
lowDensityValues=value_display';
highDensityValues = lowDensityValues(idx);


figure;
subplot(121)
scatter3(lowDensityCoords(:,1), lowDensityCoords(:,2), lowDensityCoords(:,3), 10, lowDensityValues, ...
        "filled", 'MarkerFaceAlpha',0.9);
view(0,90); colorbar;
subplot(122)
scatter3(highDensityCoords(:,1), highDensityCoords(:,2), highDensityCoords(:,3), 10, highDensityValues, ...
        "filled", 'MarkerFaceAlpha',0.9);
view(0,90); colorbar;


%
% F = scatteredInterpolant(X, Y, Z, puissance, 'natural'); % Ou 'linear', 'cubic'
% interp_values = F(Vertices(:,1), Vertices(:,2), Vertices(:,3));
sum(isnan(highDensityValues))
interp_values(isnan(highDensityValues)) = 0;

% interp_values(isnan(highDensityValues)) = 0;

writematrix(highDensityValues, 'scalp_surface.txt');


BrainNet