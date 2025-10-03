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
offset2 = 0;
%%

%%
for i = 1:length(Signal_select_Beat.QRS_Onset_auto)
    QRS = Signal_select_Beat.QRS_Onset_auto(i)+offset2:Signal_select_Beat.QRS_Onset_auto(i)+Signal_select_Beat.QRS_range;
    QRST = Signal_select_Beat.QRS_Onset_auto(i)+offset2-offset:Signal_select_Beat.QRS_Onset_auto(i)+Signal_select_Beat.Beat_range;
    
    mdVdT = getTimesElectrogram(heartpot,QRS,[]); % filter window up to 9 max
    mdVdTreal = getTimesElectrogram(rec.sock.Ve_filtered,QRS,[]);
    

    AT(i,:) = mdVdT - QRS(i);
    AT_ms(i,:) = (mdVdT - QRS(i))/(freq/1000);
    
    AT_JD(i,:) = depolMapJD(triHeart,heartpot(:,QRS));

    AT_real(i,:) = mdVdTreal - QRS(i);
    AT_ms_real(i,:) = (mdVdTreal - QRS(i))/(freq/1000);
   

    % Signal_select_Beat.AT(i,:) = mdVdT - QRS(i);
    % Signal_select_Beat.AT_ms(i,:) = (mdVdT - QRS(i))/(freq/1000);

 
    % mdVdT(rec.sock.BadChannels) = nan;
    badbyeeye = [];
    mdVdT(badbyeeye) = nan; % Changes for eah pacing site probably

    % Signal_select_Beat.badATLeads = badbyeeye;
   badATLeads = badbyeeye;

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

%save(beatdata, 'Signal_select_Beat')

%% Calculate the peak variance for the heart potential

range = maxamp - minamp;
stop = size(range,1)-1;
denominator = max(range(2:stop,:), [], 1);
denominator(denominator == 0) = 1;  % Replace zeros with 1 to avoid division by zero
peakVar = std(range(2:stop,:) ./ denominator);
peakVar = min(peakVar,0.1);
data{1} =(range(2:stop,:)./ denominator)';
data{2} = (range(2:stop,:) ./ denominator)';
 data{1} = heartpot;

close all
Dessin(triHeart,[],peakVar); hold on
DessinAT_Pts(triHeart.X,data,peakVar)
savefig([pacingSiteName,'.fig']);
end

