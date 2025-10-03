clc;
clear;

%% === INITIALIZATION ===
paceCHU = { 'CHU992', 'CHU1013', 'CHU1020', ...
           'CHU1027', 'CHU1056', 'CHU1058', 'CHU1062', ...
           'CHU1069', 'CHU1077', 'CHU1088'};
log_file = 'C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical\interpcsv\pngNN\infolog.log';
[pace, top_vals, bot_vals,rest_vals] = parse_resp_log(log_file);

% Output matrices for each quantity
lambda_L_data = [];
lambda_GCV_data = [];
lambda_RGCV_data = [];
lambda_CRESO_data = [];
lambda_U_data = [];
lambda_L_fix_data = [];
rho_data = [];
eta_data = [];
column_labels = {};



%% === MAIN LOOP ===
for i = 1:length(paceCHU)
    disp(paceCHU{i})
    idx = find(ismember(pace, ['interpSig',paceCHU{i}]));
    pathbase = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical\',paceCHU{i},'\'];
    path1 = [pathbase,paceCHU{i},'.mat'];
    path2 = fullfile(pathbase, 'transferA.mat');
    path3 = fullfile('C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical',paceCHU{i},'\Signal_select_auto\',['interpSig',paceCHU{i},'.mat']);
    if exist(path1, 'file') && exist(path2, 'file') && exist(path3, 'file')
        load(path1);
        load(path2);
        load(path3);
    else
        % One or both files do not exist — skip or continue
        fprintf('Skipping: %s or %s or %s not found.\n', path1, path2,path3);
        continue;  % or use 'continue' if inside a loop
    end

    if length(top_vals{idx}) > length(bot_vals{idx})
        exhalationpeaks = top_vals{idx};
        inhlationpeaks = bot_vals{idx};
    else
        inhlationpeaks = top_vals{idx};
        exhalationpeaks = bot_vals{idx};
    end
    rest_peaks = rest_vals{idx};

    [timeInterval, ~, respcycle] = timeintervalgenerator2(Signal_select_Beat, inhlationpeaks, exhalationpeaks,rest_peaks);
    [~, meanhelper] = timeintervalgenerator(Signal_select_Beat);
    N_time = length(timeInterval);
    nLeadsBody = size(interpSig, 1);
    disp([respcycle;Signal_select_Beat.QRS_Onset_auto'])
end
