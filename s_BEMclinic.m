clc
clear all
close all

% Load Geometries
patientNumber = 'CHU1088'; % Specify the patient number
load(fullfile('C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical\', patientNumber, [patientNumber, '.mat']));
%%

% 
% triH = TriRep(triHeart.Triangulation,triHeart.X);
% triT = TriRep(triTorso.Triangulation,triTorso.X);

triHH = checkAndFixMesh(triH);

heartX = triHH.X;
heartTri = triHH.Triangulation;


Dessin(triHH,[],1:size(triHH,1));
pause;
triTT =checkAndFixMesh(triT);

torsoX = triTT.X;
torsoTri = triTT.Triangulation;
%%
heartModel = struct('node', heartX', 'face', heartTri', 'sigma', [0.22 0]);
torsoModel = struct('node', torsoX', 'face', torsoTri', 'sigma', [0 0.22]); 
completeModel.surface = {torsoModel, heartModel};

% Run BEM code to compute transfer matrix
transferA = bemMatrixPP(completeModel); 
% Save the transfer matrix to the loaded subfolder
save(fullfile('C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical\', patientNumber, 'transferA.mat'), 'transferA');