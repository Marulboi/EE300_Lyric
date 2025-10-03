clc
clear all 
close all
%%
foldername = 'C:\Users\emire\OneDrive\Desktop\resp\Pig\Pig5\Signal_select_auto\';
matFiles = dir(fullfile(foldername, '*.mat'));

for k = 1:length(matFiles)
   fileName = matFiles(k).name;
   fullPath = fullfile(foldername, fileName);
   beatdata = fullPath;
pacingSiteName = extractAfter(beatdata, 'Signal_select_auto\');
pacingSiteName = extractBefore(pacingSiteName, '_Signal_select_Beat.mat');

load(['C:\Users\emire\OneDrive\Desktop\resp\Pig\Pig5\Signals\',pacingSiteName,'.mat'])


load(beatdata)
load('C:\Users\emire\OneDrive\Desktop\resp\Pig\Pig5\pig5_correctGeo_corr.mat')
%%
% triHeart = TriRep(meshes.sock_capped.faces,meshes.sock_capped.verts);
 triHeart = TriRep(meshes.sock.faces,meshes.sock.verts);



heartpot = rec.sock.Ve_filtered;
freq = 2048;
offset=10;
offset2 = +20;
%%

%%
for i = 1:length(Signal_select_Beat.QRS_Onset)
    QRS = Signal_select_Beat.QRS_Onset(i)+offset2:Signal_select_Beat.QRS_Onset(i)+Signal_select_Beat.QRS_range;
    
    mdVdT = getTimesElectrogram(heartpot,QRS,[]); % filter window up to 9 max
    mdVdTreal = getTimesElectrogram(rec.sock.Ve_filtered,QRS,[]);
    

    AT(i,:) = mdVdT - QRS(i);
    AT_ms(i,:) = (mdVdT - QRS(i))/(freq/1000);
    
    AT_JD(i,:) = depolMapJD(triHeart,heartpot(:,QRS));

    AT_real(i,:) = mdVdTreal - QRS(i);
    AT_ms_real(i,:) = (mdVdTreal - QRS(i))/(freq/1000);
   

    % Signal_select_Beat.AT(i,:) = mdVdT - QRS(i);
    % Signal_select_Beat.AT_ms(i,:) = (mdVdT - QRS(i))/(freq/1000);

 
    mdVdT(rec.sock.BadChannels) = nan;


    % Signal_select_Beat.badATLeads = badbyeeye;
    badATLeads = (AT(i,:)==1);

    interpAt = medianSpatialFilter(triHeart.X(isnan(mdVdT),:), triHeart.X(~isnan(mdVdT),:),mdVdT(~isnan(mdVdT))',10); % 30 max 
    inter_mdVdT =mdVdT;
    inter_mdVdT(isnan(mdVdT)) = interpAt;

    % data{1} = heartpot(:,QRST);
    % data{2} = inter_mdVdT+offset;
    % 
    % close all
    % Dessin(triHeart,[],inter_mdVdT); hold on
    % DessinAT_Pts(triHeart.X,data,inter_mdVdT);   
    % % % % figure;
    % % % % plot(gradient(heartpot(136,QRST)));
    % pause;

 end


end