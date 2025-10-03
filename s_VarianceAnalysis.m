clc
clear all 
close all
%%
foldername = 'C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\Signal_select_auto\';
matFiles = dir(fullfile(foldername, '*.mat'));

for k = 1:length(matFiles)
   fileName = matFiles(k).name;
   fullPath = fullfile(foldername, fileName);
   beatdata = fullPath;
pacingSiteName = extractAfter(beatdata, 'Signal_select_auto\');
pacingSiteName = extractBefore(pacingSiteName, '_Signal_select_Beat.mat');

load(['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\Signals\',pacingSiteName,'.mat'])
load("C:\Users\emir.ege-nemutlu\Desktop\resp\inverse solutions.mat")


load(beatdata)
load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\pig5_correctGeo_corr.mat')
%%
% triHeart = TriRep(meshes.sock_capped.faces,meshes.sock_capped.verts);
 triHeart = TriRep(meshes.sock.faces,meshes.sock.verts);



heartpot = rec.sock.Ve_filtered;
freq = 2048;
offset=10;
offset2 = 0;

%%
for i = 1:length(Signal_select_Beat.QRS_Onset_auto)
    QRS = Signal_select_Beat.QRS_Onset_auto(i)+offset2:Signal_select_Beat.QRS_Onset_auto(i)+Signal_select_Beat.QRS_range;
    QRST = Signal_select_Beat.QRS_Onset_auto(i)+offset2-offset:Signal_select_Beat.QRS_Onset_auto(i)+Signal_select_Beat.Beat_range;
    
    maxamp(i,:) = max(heartpot(:,QRS),[],2);
    minamp(i,:) = min(heartpot(:,QRS),[],2);  

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
