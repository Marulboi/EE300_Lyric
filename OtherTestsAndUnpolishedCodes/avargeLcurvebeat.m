clc;
%close all;
clear;

load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\Signals\EndoPacingSite-LV-AnteriorBase.mat');
load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\pig5_correctGeo_corr.mat');
load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\Signal_select_auto\EndoPacingSite-LV-AnteriorBase_Signal_select_Beat.mat')


start = Signal_select_Beat.QRS_Onset_auto(1);
stop = Signal_select_Beat.QRS_Onset_auto(end) + Signal_select_Beat.QRS_range;
timeInterval = start:stop;

nLeadsBody = size(rec.vest.Ve_filtered,1);
brokenLeads = zeros(nLeadsBody,1);
brokenLeads(rec.vest.BadChannels) = 1;
lambda_range = logspace(-5, 1, 500);
Asub = Trf_HT_coarse((brokenLeads == 0),:); 
[U, S, V] = svd(Asub);
S_diag = diag(S); 
solsize = size(S,2);
figure;
    n = length(Signal_select_Beat.QRS_Onset_auto);
    cols = ceil(sqrt(n));
    rows = ceil(n / cols);

tiledlayout(cols,rows);
averageres = zeros(n,2);
averagesol = zeros(n,2);
k=1;
for i = Signal_select_Beat.QRS_Onset_auto' 
    nexttile;
    % start = (i+Signal_select_Beat.BeatOffset);
    % stop = (i+Signal_select_Beat.BeatOffset+Signal_select_Beat.Beat_range);
    start = (i);
    stop = (i+Signal_select_Beat.QRS_range);

    Sig_sub = rec.vest.Ve_filtered((brokenLeads == 0),(start:stop));
    [sol_norm_avg, res_norm_avg, eta, rho] = avargeLcurve(U,S_diag,lambda_range,Sig_sub);
    hold on
    for j = 1:size(eta,2)
        loglog(rho(:,j), eta(:,j), 'Color', [0.5, 0.5, 0.5, 0.2]); % RGBA with alpha = 0.2
    end
    averageres(k,:) = [i,(res_norm_avg(1))];
    averagesol(k,:) = [i, (sol_norm_avg(1))];
    k = k +1;
    loglog(res_norm_avg,sol_norm_avg,'Color','red');

    margin = 0.1;  % 10% padding
    y_min = min(sol_norm_avg);
    y_max = max(sol_norm_avg);
    y_range = y_max - y_min;
    ylim([y_min - margin*y_range, y_max + margin*y_range]);
    x_min = min(res_norm_avg);
    x_max = max(res_norm_avg);
    x_range = x_max - x_min;
    xlim([x_min - margin * x_range, x_max + margin * x_range]);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    title(sprintf('Plot %d', i));
end

%%
figure;
% === Subplot 1: Residual and Solution Norms ===
subplot(2,1,1)
yyaxis left;
semilogy(averageres(:,1), averageres(:,2), 'b', 'DisplayName', 'Residual Norm');
ylabel('Residual Norm');
xlabel('Time Index (QRS Onset)');
hold on;

yyaxis right;
semilogy(averagesol(:,1), averagesol(:,2), 'r', 'DisplayName', 'Solution Norm');
ylabel('Solution Norm');

title('Residual and Solution Norm Evolution Over Time');
legend('Location','best');
xlim([min(timeInterval), max(timeInterval)]);  % Align x-axes

% === Subplot 2: Raw Signal Segment with Envelope ===
subplot(2,1,2)
signal = rec.vest.Ve_filtered(33, timeInterval);
plot(timeInterval, signal, 'Color', [0.5, 0.5, 0.5, 0.2])
hold on;

% Compute the envelope
[upper_env, ~] = envelope(signal, 500, 'peak');  % 50: window size, can be tuned
plot(timeInterval, upper_env, 'r', 'LineWidth', 1.5)

xlabel('Time Index (QRS Onset)');
ylabel('Amplitude (Lead 2)');
title('Filtered Signal on Lead 2 with Peak Envelope');
xlim([min(timeInterval), max(timeInterval)]);
legend('Raw Signal', 'Envelope (Peaks)');