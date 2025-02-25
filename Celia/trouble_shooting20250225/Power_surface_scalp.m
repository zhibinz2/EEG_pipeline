clear
M = readmatrix('E:\DATASET\Script\BrainNet\Data\SurfTemplate\BrainMesh_ICBM152.txt');
M(81925:end,:)=[];
Vertices = M;
load('E:\DATASET\4_Power_Calculation\Control.mat')

load('E:\DATASET\4_Source_reconstructed\Control\ARAT.mat')

X=corti_ave_source_coor(:,1);
Y=corti_ave_source_coor(:,2);
Z=corti_ave_source_coor(:,3);

puissance=mean(POWER.DELTA.Sum,2);

interp_values = griddata(X, Y, Z, puissance, Vertices(:,1), Vertices(:,2), Vertices(:,3), 'nearest');
%
% F = scatteredInterpolant(X, Y, Z, puissance, 'natural'); % Ou 'linear', 'cubic'
% interp_values = F(Vertices(:,1), Vertices(:,2), Vertices(:,3));

interp_values(isnan(interp_values)) = 0;

writematrix(interp_values, 'scalp_surface.txt');
