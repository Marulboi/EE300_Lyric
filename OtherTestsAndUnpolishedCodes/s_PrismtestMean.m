clc;
clear;

%% === INITIALIZATION ===
pig_num = num2str(1);
load(['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig',pig_num,'\ForwardM\ForwMat_Sock_HInTOut.mat'],'Trf_HT_coarse');
log_file = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig',pig_num,'\SignalsCSV\pngN\infolog.log'];
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


%%
% Output directory
outdir = ['ExportedCSVmeansNEWmeth',pig_num];

mkdir(outdir);


%% === MAIN LOOP ===
for i = 1:length(pace)
    path1 = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig',pig_num,'\Signals\', pace{i}, '.mat'];
    if length(top_vals{i}) < 2 || length(bot_vals{i}) < 2
        disp(['Skipped1 :', pace{i}]);
        continue
    end
    disp(['Computing :', pace{i}]);

    [folder, name, ~] = fileparts(path1);
    folder2 = strrep(folder, 'Signals', 'Signal_select_auto');
    filename2 = [name '_Signal_select_Beat.mat'];
    path2 = fullfile(folder2, filename2);
    if exist(path1, 'file') && exist(path2, 'file')
        load(path1);
        load(path2);
    else
        % One or both files do not exist — skip or continue
        fprintf('Skipping: %s or %s not found.\n', path1, path2);
        continue;  % or use 'continue' if inside a loop
    end
    load(path1); load(path2);

    if length(top_vals{i}) > length(bot_vals{i})
        exhalationpeaks = top_vals{i};
        inhlationpeaks = bot_vals{i};
    else
        inhlationpeaks = top_vals{i};
        exhalationpeaks = bot_vals{i};
    end
    rest_peaks = rest_vals{i};

    [timeInterval, ~, respcycle] = timeintervalgenerator2(Signal_select_Beat, inhlationpeaks, exhalationpeaks,rest_peaks);
    [~, meanhelper] = timeintervalgenerator(Signal_select_Beat);
    N_time = length(timeInterval);
    nLeadsBody = size(rec.vest.Ve_filtered, 1);
    brokenLeads = zeros(nLeadsBody, 1);
    brokenLeads(rec.vest.BadChannels) = 1;

    lambda_range = logspace(-6, 0, 500);
    robustness = 0.5;
    Asub = Trf_HT_coarse((brokenLeads == 0), :);
    Sig_sub = rec.vest.Ve_filtered((brokenLeads == 0), timeInterval);

    [U, S, ~] = svd(Asub);
    S_diag = diag(S);

    lambda_L = zeros(N_time, 1);
    lambda_GCV = zeros(N_time, 1);
    lambda_RGCV = zeros(N_time, 1);
    lambda_CRESO = zeros(N_time, 1);
    lambda_U = zeros(N_time, 1);
    lambda_L_fix = zeros(N_time, 1);
    rho = zeros(1, N_time);
    eta = zeros(1, N_time);

    parfor t1 = 1:N_time
        b_t = Sig_sub(:, t1);
        [lambda_L(t1), ~, ~] = l_curve(U, S_diag, b_t);
        norm_svd_func = @(lambda) norms_svd(U, S_diag, b_t, lambda);
        [r_tmp, e_tmp] = norms_svd(U, S_diag, b_t, lambda_range);
        rho(t1) = r_tmp(end);
        eta(t1) = e_tmp(end);
        lambda_L_fix(t1) = my_l_curve(norm_svd_func, lambda_range);
        lambda_GCV(t1) = gcv(U, S_diag, b_t);
        lambda_RGCV(t1) = rgcv_svd(U, S, b_t, lambda_range, robustness);
        [lambda_CRESO(t1), ~] = creso(U, S_diag, b_t, lambda_range);
        lambda_U(t1) = newUcurve(U, S_diag, b_t, lambda_range);
    end

    lambda_L_piecewise     = zeros(1, 3);
    lambda_GCV_piecewise   = zeros(1, 3);
    lambda_RGCV_piecewise  = zeros(1, 3);
    lambda_CRESO_piecewise = zeros(1, 3);
    lambda_U_piecewise     = zeros(1, 3);
    lambda_L_fix_piecewise = zeros(1, 3);
    rho_piecewise          = zeros(1, 3);
    eta_piecewise          = zeros(1, 3);
    
    if isempty(meanhelper)
        disp(['skipping:',pace{i}])
        continue
    end
    for b = 1:meanhelper(end)
        lambda_L_piecewise(meanhelper == b)     = mean(lambda_L(meanhelper == b));
        lambda_GCV_piecewise(meanhelper == b)   = mean(lambda_GCV(meanhelper == b));
        lambda_RGCV_piecewise(meanhelper == b)  = mean(lambda_RGCV(meanhelper == b));
        lambda_CRESO_piecewise(meanhelper == b) = mean(lambda_CRESO(meanhelper == b));
        lambda_U_piecewise(meanhelper == b)     = mean(lambda_U(meanhelper == b));
        lambda_L_fix_piecewise(meanhelper == b) = mean(lambda_L_fix(meanhelper == b));
        rho_piecewise(meanhelper == b) = mean(rho(meanhelper == b));
        eta_piecewise(meanhelper == b) = mean(eta(meanhelper == b));
    end

    current_lambda_L = [];
    current_lambda_GCV = [];
    current_lambda_RGCV = [];
    current_lambda_CRESO = [];
    current_lambda_U = [];
    current_lambda_L_fix = [];
    current_rho = [];
    current_eta = [];

    for k = Signal_select_Beat.QRS_Onset_auto'
        idx = find(timeInterval == k + 2,1);
        current_lambda_L(end+1)     = lambda_L_piecewise(idx);
        current_lambda_GCV(end+1)   = lambda_GCV_piecewise(idx);
        current_lambda_RGCV(end+1)  = lambda_RGCV_piecewise(idx);
        current_lambda_CRESO(end+1) = lambda_CRESO_piecewise(idx);
        current_lambda_U(end+1)     = lambda_U_piecewise(idx);
        current_lambda_L_fix(end+1) = lambda_L_fix_piecewise(idx);
        current_rho(end+1)          = rho_piecewise(idx);
        current_eta(end+1)          = eta_piecewise(idx);
    end

       lambda_L_data = [lambda_L_data;[mean(current_lambda_L(respcycle==1)), mean(current_lambda_L(respcycle==2)), mean(current_lambda_L(respcycle==3))]];
       lambda_GCV_data = [lambda_GCV_data;[mean(current_lambda_GCV(respcycle==1)), mean(current_lambda_GCV(respcycle==2)), mean(current_lambda_GCV(respcycle==3))]];
       lambda_RGCV_data = [lambda_RGCV_data;[mean(current_lambda_RGCV(respcycle==1)), mean(current_lambda_RGCV(respcycle==2)), mean(current_lambda_RGCV(respcycle==3))]];
       lambda_CRESO_data = [lambda_CRESO_data;[mean(current_lambda_CRESO(respcycle==1)), mean(current_lambda_CRESO(respcycle==2)), mean(current_lambda_CRESO(respcycle==3))]];
       lambda_U_data = [lambda_U_data;[mean(current_lambda_U(respcycle==1)), mean(current_lambda_U(respcycle==2)), mean(current_lambda_U(respcycle==3))]];
       lambda_L_fix_data = [lambda_L_fix_data;[mean(current_lambda_L_fix(respcycle==1)), mean(current_lambda_L_fix(respcycle==2)), mean(current_lambda_L_fix(respcycle==3))]];
       rho_data = [rho_data;[mean(current_rho(respcycle==1)), mean(current_rho(respcycle==2)), mean(current_rho(respcycle==3))]];
       eta_data = [eta_data;[mean(current_eta(respcycle==1)), mean(current_eta(respcycle==2)), mean(current_eta(respcycle==3))]];

    disp('==================')
    disp(eta_data)
    disp('==================')
    
    column_labels = [column_labels, pace{i}];




end

%% Save outputs to folder
lambda_L_data_T =  lambda_L_data;
lambda_GCV_data_T =  lambda_GCV_data;
lambda_RGCV_data_T =  lambda_RGCV_data;
lambda_CRESO_data_T =  lambda_CRESO_data;
lambda_U_data_T = lambda_U_data;
lambda_L_fix_data_T =  lambda_L_fix_data;
rho_data_T =  rho_data;
eta_data_T =  eta_data;

toplabel = {"inhalation","exhalation","rest"};
leftlabel = [" ",column_labels]';
leftlabel = cellstr(leftlabel);
lambda_L_data_T = [toplabel; num2cell(lambda_L_data_T)];
lambda_L_data_T = [leftlabel, lambda_L_data_T];

lambda_GCV_data_T = [toplabel; num2cell(lambda_GCV_data_T)];
lambda_GCV_data_T = [leftlabel, lambda_GCV_data_T];

lambda_RGCV_data_T = [toplabel; num2cell(lambda_RGCV_data_T)];
lambda_RGCV_data_T = [leftlabel, lambda_RGCV_data_T];

lambda_CRESO_data_T = [toplabel; num2cell(lambda_CRESO_data_T)];
lambda_CRESO_data_T = [leftlabel, lambda_CRESO_data_T];

lambda_U_data_T = [toplabel; num2cell(lambda_U_data_T)];
lambda_U_data_T = [leftlabel, lambda_U_data_T];

lambda_L_fix_data_T = [toplabel; num2cell(lambda_L_fix_data_T)];
lambda_L_fix_data_T = [leftlabel, lambda_L_fix_data_T];

rho_data_T = [toplabel; num2cell(rho_data_T)];
rho_data_T = [leftlabel, rho_data_T];

eta_data_T = [toplabel; num2cell(eta_data_T)];
eta_data_T = [leftlabel, eta_data_T];
%%

writecell(lambda_L_data_T,     fullfile(outdir, 'Lambda_L.csv'));
writecell(lambda_GCV_data_T,   fullfile(outdir, 'Lambda_GCV.csv'));
writecell(lambda_RGCV_data_T,  fullfile(outdir, 'Lambda_RGCV.csv'));
writecell(lambda_CRESO_data_T, fullfile(outdir, 'Lambda_CRESO.csv'));
writecell(lambda_U_data_T,     fullfile(outdir, 'Lambda_U.csv'));
writecell(lambda_L_fix_data_T, fullfile(outdir, 'Lambda_L_fix.csv'));
writecell(rho_data_T,          fullfile(outdir, 'Rho.csv'));
writecell(eta_data_T,          fullfile(outdir, 'Eta.csv'));

disp('All files saved successfully to ExportedCSVs/');
