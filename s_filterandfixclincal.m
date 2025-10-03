clc
clear all
close all

load("C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical\CHU1013\Patient.mat")
load("C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical\CHU1013\ci.mat")

freq = 1000;

a = {[16 17 20 ], [89 90 91], [123 124 125 126], ...
            49, 73, 77, 101, 104, 115};

[~, indices] = cellfun(@(x) ismember(x,ChannelList), a, 'UniformOutput', false);

%%
% triHeart.Triangulation(size(triHeart.Triangulation,1)+1,:) = [939,938,895];

triHe = TriRep(triHeart.Triangulation,triHeart.X);
triTo = TriRep(triTorso.Triangulation,triTorso.X);

% Dessin(triHe,[],1:size(triHe,1));

% colormap winter



triH = checkAndFixMesh(triHe);
triT = checkAndFixMesh(triTo);

heartModel = struct('node', triH.X, 'face', triH.Triangulation, 'sigma', [0.22 0]);
torsoModel = struct('node', triT.X, 'face', triT.Triangulation, 'sigma', [0 0.22]); 

% completeModel.surface = {torsoModel, heartModel};
% 
% transferA = bemMatrixPP(completeModel); 
% 
% [U, S, V] = svd(transferA);



%%

    ECG = signalObject(sig,'samp',freq);           % Create ECG object with sampling rate

    ECG.applyLowPassFilter(150,2,'FIR');                     % Apply low-pass filter

    ECG.applyNotchFilter(50,3,'Fourier');

    ECG.removeBaseline('savitzkygolay', [3 2500]);
 
    plot(ECG,1);
    
    %%
    temp = ECG.processedSignal;    
%%

badchannels = [4,5,7,12,13,16:23,36,37,39,47,61,88,89,93:94,113:114,116:122,131:135,212,213,244,249];
    
    temp(badchannels,:) = nan;

    % 
    %  inter = signalObject(temp,'samp',freq);    
    % plot(inter,1);
%%
    
    recordLength = 1:size(sig,1);
    % beforeint = NaN(252,recordLength);
    % tempchannels = ismember(1:252,ChannelList);
    % beforeint(tempchannels,:) = temp;
    
    goodleads = recordLength(any(~isnan(temp),2));
    badleads = recordLength(any(isnan(temp),2));
    temp(badleads,:) = 0;


    %% 
    [ interpSig ] = interpolateInvFwd( triT, badleads, goodleads,temp, 0.8 );
    inter = signalObject(interpSig,'samp',freq);    
    plot(inter,1);


  %%
  % 
  save('C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical\CHU1013\CHU1013.mat', 'sig', 'interpSig', 'triH', 'triT','badleads');
  % writematrix(interpSig, 'C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical\CHU1088\interpSig.csv');
  % 
  % 
