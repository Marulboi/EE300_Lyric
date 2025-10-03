clc;
%close all;
clear;

load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\Signals\EndoPacingSite-LV-AnteriorBase.mat');
load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\ForwardM\ForwMat_Sock_HInTOut.mat','Trf_HT_coarse');
load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\Signals_Beats\EndoPacingSite-LV-AnteriorBase_Signal_select_Beat.mat')


nLeadsBody = size(rec.vest.Ve_filtered,1);
brokenLeads = zeros(nLeadsBody,1);
brokenLeads(rec.vest.BadChannels) = 1;
lambda_range = logspace(-5, 1, 500);
Asub = Trf_HT_coarse((brokenLeads == 0),:); 
[U, S, V] = svd(Asub);
S_diag = diag(S); 
solsize = size(S,2);
start = Signal_select_Beat.QRS_Onset(1) + Signal_select_Beat.BeatOffset;
stop = Signal_select_Beat.QRS_Onset(end) + Signal_select_Beat.BeatOffset+Signal_select_Beat.Beat_range;

timeInterval = start:stop;

Sig_sub = rec.vest.Ve_filtered((brokenLeads == 0),timeInterval);
[sol_norm_avg, res_norm_avg, eta, rho] = avargeLcurve(U,S_diag,lambda_range,Sig_sub);
    

figure;
hold on;
for i = 1:50:length(timeInterval)
    t = timeInterval(i);                      % time value

    plot3( t*ones(size(lambda_range)),rho(:,i), eta(:,i));  % plot line at Y = time
end
ylabel('eta (Sol norm) (log)');
xlabel('Time');
zlabel('rho (res) (log)');
title('3D Lines with Log X and Z');
grid on;
view(3);


set(gca, 'YScale', 'log', 'ZScale', 'log');