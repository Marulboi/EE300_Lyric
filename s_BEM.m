clc
clear all
close all

% Load Geometries

% Define torso nodes and elements
[torsoX,torsoTri] = loadtri('BSMUMC050Plain_thoraxnew.tri');
% torsoX = sim.itorso.verts;
% torsoTri = sim.itorso.faces;
tri = TriRep(torsoTri,torsoX)

% Define heart nodes and elements
[ heartX,heartTri] = loadtri('BSMUMC050Plain_ventricularpericardium.tri');
triH = TriRep(heartTri,heartX);
triNew = removeNodesMesh(triH,[251 252]);

vtkWriteBinary(triNew,[],'triHeart');
triH = checkAndFixMesh(triNew);
heartX = triH.X;
heartTri = triH.Triangulation;
% heartX = sim.peri3small.verts;
% heartTri = sim.peri3small.faces;
% 
% Define structures of heart and torso in format require for BEM code
% I think you can add other elements (lungs etc) into this and the
% conductivity is defined by sigma. This is homogeneous
heartModel = struct('node', heartX', 'face', heartTri', 'sigma', [0.22 0]);
torsoModel = struct('node', torsoX', 'face', torsoTri', 'sigma', [0 0.22]); 
completeModel.surface = {torsoModel, heartModel};

% Run BEM code to compute transfer matrix
transferA = bemMatrixPP(completeModel); 

% % Use only nodes corresponding to electrodes
% Electrodes = dlmread('BSMUMC050Plain.el');
% Electrodes = sort(Electrodes(2:end,2));
% transferA = transferA(Electrodes,:);
% [matU,vectS,matV] = csvd(transferA);

% In standard formulation, Vector b contains torso potentials only
vectB = TORSOPOTENTIALS;
nSamp = size(vectB,2);
% L-curve method : draw average l-curve to select reg param
hWait = waitbar(0, 'Calculation regularization parameter (L-curve)');
rho =[]; eta = [];
regvalC = zeros(nSamp, 1);
for i=1:nSamp
    waitbar(i/nSamp, hWait);
    [regvalC(i),rho(i,:), eta(i,:), regval(i,:)]= l_curve(matU,vectS,vectB(:, i), 'tikh');
   %plot(regval(1,:))
    %regvalC(i) = l_cornerJD(rho(i,:),eta(i,:),regval(i,:), 10^4);
   % regvalC(i) = l_cornerJD(rho(i,:),eta(i,:),regval(i,:))
end
close(hWait);
avrho = mean(rho, 1); aveta = mean(eta, 1); avreg = mean(regval, 1);

% Regularization parameter is the median over the whole time series (median
% to avoid outliers)
lc_lambda = median(regvalC);

% Ask user input for reg parameter lambda
lambda = lc_lambda * ones(nSamp, 1);

% Compute heart potentials
vectA = zeros(size(matV, 1), nSamp);
for i=1:nSamp, vectA(:, i) = tikhonov(matU, vectS, matV, vectB(:, i), lambda(i)); end
        