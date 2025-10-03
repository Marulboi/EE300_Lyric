clc;
clear;

%% === INITIALIZATION ===
paceCHU = { 'CHU976','CHU992', 'CHU1013', 'CHU1020', ...
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


%%
% Output directory
outdir = 'ExportedCSVmeansClinicalNEW';

mkdir(outdir);


%% === MAIN LOOP ===
for i = 1:length(paceCHU)

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


    lambda_range = logspace(-6, 0, 500);
    robustness = 0.5;
    Asub = transferA;
    Sig_sub = interpSig(:, timeInterval);

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
        disp(['skipping:',paceCHU{i}])
        continue
    end
    for b = 1:3
        lambda_L_piecewise(b)     = mean(lambda_L(respcycle == b));
        lambda_GCV_piecewise(b)   = mean(lambda_GCV(respcycle == b));
        lambda_RGCV_piecewise(b)  = mean(lambda_RGCV(respcycle == b));
        lambda_CRESO_piecewise(b) = mean(lambda_CRESO(respcycle == b));
        lambda_U_piecewise(b)     = mean(lambda_U(respcycle == b));
        lambda_L_fix_piecewise(b) = mean(lambda_L_fix(respcycle == b));
        rho_piecewise(b) = mean(rho(respcycle == b));
        eta_piecewise(b) = mean(eta(respcycle == b));
    end

       lambda_L_data = [lambda_L_data;[lambda_L_piecewise(1),lambda_L_piecewise(2),lambda_L_piecewise(3)]];
       lambda_GCV_data = [lambda_GCV_data;[lambda_GCV_piecewise(1),lambda_GCV_piecewise(2),lambda_GCV_piecewise(3)]];
       lambda_RGCV_data = [lambda_RGCV_data;[lambda_RGCV_piecewise(1),lambda_RGCV_piecewise(2),lambda_RGCV_piecewise(3)]];
       lambda_CRESO_data = [lambda_CRESO_data;[lambda_CRESO_piecewise(1),lambda_CRESO_piecewise(2),lambda_CRESO_piecewise(3)]];
       lambda_U_data = [lambda_U_data;[lambda_U_piecewise(1),lambda_U_piecewise(2),lambda_U_piecewise(3)]];
       lambda_L_fix_data = [lambda_L_fix_data;[lambda_L_fix_piecewise(1),lambda_L_fix_piecewise(2),lambda_L_fix_piecewise(3)]];
       rho_data = [rho_data;[rho_piecewise(1),rho_piecewise(2),rho_piecewise(3)]];
       eta_data = [eta_data;[eta_piecewise(1),eta_piecewise(2),eta_piecewise(3)]];


    disp('==================')
    disp(eta_data)
    disp('==================')
    
    column_labels = [column_labels, paceCHU{i}];




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
