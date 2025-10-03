clc;
clear;

%% === INITIALIZATION ===
paceCHU = {'CHU992', 'CHU1013', 'CHU1020', ...
           'CHU1027', 'CHU1056', 'CHU1058', 'CHU1062', ...
           'CHU1069', 'CHU1077', 'CHU1088','CHU976'};
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


leftlabel = {' ';'Inhalation';'Exhalation';'rest'};
%%
% Output directory
outdir = 'ExportedCSVClinicalN';

mkdir(outdir);


%% === MAIN LOOP ===
for i = 1:length(paceCHU)
    disp(paceCHU{i})
    id = find(ismember(pace, ['interpSig',paceCHU{i}]));
    pathbase = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical\',paceCHU{i},'\'];
    path1 = [pathbase,paceCHU{i},'.mat'];
    path2 = fullfile(pathbase, 'transferA.mat');
    path3 = fullfile('C:\Users\emir.ege-nemutlu\Desktop\resp\Clinical',paceCHU{i},'\Signal_select_auto\',['interpSig',paceCHU{i},'.mat']);
    
    if exist(path1, 'file')
        load(path1);
    else
        fprintf('Skipping: %s not found.\n', path1);
        continue;
    end

    if exist(path2, 'file')
        load(path2);
    else
        fprintf('Skipping: %s not found.\n', path2);
        continue;
    end

    if exist(path3, 'file')
        load(path3);
    else
        fprintf('Skipping: %s not found.\n', path3);
        continue;
    end


    if length(top_vals{id}) > length(bot_vals{id})
        exhalationpeaks = top_vals{id};
        inhlationpeaks = bot_vals{id};
    else
        inhlationpeaks = top_vals{id};
        exhalationpeaks = bot_vals{id};
    end
    rest_peaks = rest_vals{id};

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
    disp('loading complete parfor everything')
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

    lambda_L_piecewise     = zeros(size(lambda_L));
    lambda_GCV_piecewise   = zeros(size(lambda_GCV));
    lambda_RGCV_piecewise  = zeros(size(lambda_RGCV));
    lambda_CRESO_piecewise = zeros(size(lambda_CRESO));
    lambda_U_piecewise     = zeros(size(lambda_U));
    lambda_L_fix_piecewise = zeros(size(lambda_L_fix));
    rho_piecewise          = zeros(size(rho));
    eta_piecewise          = zeros(size(eta));
    
    if isempty(meanhelper)
        disp(['skipping:',paceCHU{i}])
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
    
    QRS_Onset_auto = Signal_select_Beat.QRS_Onset_auto';
    QRS_Onset_auto(Signal_select_Beat.QRS_Onset_auto'<0)=1;

    for k = QRS_Onset_auto
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

    lambda_L_data     = padcat(lambda_L_data,     padcat(current_lambda_L(respcycle==1)', current_lambda_L(respcycle==2)', current_lambda_L(respcycle==3)'));
    lambda_GCV_data   = padcat(lambda_GCV_data,   padcat(current_lambda_GCV(respcycle==1)', current_lambda_GCV(respcycle==2)', current_lambda_GCV(respcycle==3)'));
    lambda_RGCV_data  = padcat(lambda_RGCV_data,  padcat(current_lambda_RGCV(respcycle==1)', current_lambda_RGCV(respcycle==2)', current_lambda_RGCV(respcycle==3)'));
    lambda_CRESO_data = padcat(lambda_CRESO_data, padcat(current_lambda_CRESO(respcycle==1)', current_lambda_CRESO(respcycle==2)', current_lambda_CRESO(respcycle==3)'));
    lambda_U_data     = padcat(lambda_U_data,     padcat(current_lambda_U(respcycle==1)', current_lambda_U(respcycle==2)', current_lambda_U(respcycle==3)'));
    lambda_L_fix_data = padcat(lambda_L_fix_data, padcat(current_lambda_L_fix(respcycle==1)', current_lambda_L_fix(respcycle==2)', current_lambda_L_fix(respcycle==3)'));
    rho_data          = padcat(rho_data,          padcat(current_rho(respcycle==1)', current_rho(respcycle==2)', current_rho(respcycle==3)'));
    eta_data          = padcat(eta_data,          padcat(current_eta(respcycle==1)', current_eta(respcycle==2)', current_eta(respcycle==3)'));
    
    % disp('==================')
    % disp(eta_data)
    % disp('==================')
    
    column_labels = [column_labels, paceCHU{i}];
end

%% Save outputs to folder
lambda_L_data_T = concat_blocks_by_3(lambda_L_data);
lambda_GCV_data_T = concat_blocks_by_3(lambda_GCV_data);
lambda_RGCV_data_T = concat_blocks_by_3(lambda_RGCV_data);
lambda_CRESO_data_T = concat_blocks_by_3(lambda_CRESO_data);
lambda_U_data_T = concat_blocks_by_3(lambda_U_data);
lambda_L_fix_data_T = concat_blocks_by_3(lambda_L_fix_data);
rho_data_T = concat_blocks_by_3(rho_data);
eta_data_T = concat_blocks_by_3(eta_data);

toplabel = repelem(column_labels, size(lambda_L_data, 1));

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
%%
function C = padcat(varargin)
% PADCAT - Pads and concatenates numeric vectors/matrices column-wise
% C = padcat(A, B, ...) pads with NaNs to the tallest input and concatenates horizontally.

    narginchk(1, Inf);
    maxRows = max(cellfun(@(x) size(x, 1), varargin));
    C = [];

    for k = 1:nargin
        A = varargin{k};
        if isvector(A)
            A = A(:); % make column
        end
        padSize = maxRows - size(A, 1);
        if padSize > 0
            A = [A; NaN(padSize, size(A, 2))];
        end
        C = [C, A];
    end
end
function out = concat_blocks_by_3(cellmat)
% CONCAT_BLOCKS_BY_3 - Concatenates Nx(3K) cell array into 3xM matrix
% by splitting horizontally into blocks of 3 columns, transposing each,
% and concatenating the results horizontally.
%
% INPUT:
%   cellmat - N×(3*K) cell array, each cell must contain a scalar
%
% OUTPUT:
%   out - 3×(K*N) numeric matrix

    [nRows, nCols] = size(cellmat);
    if mod(nCols, 3) ~= 0
        error('Number of columns must be a multiple of 3.');
    end

    numBlocks = nCols / 3;
    out = [];

    for i = 1:numBlocks
        block = cellmat(:, 3*(i-1)+1 : 3*i);   % N×3 cell block 
        out = [out, block'];                   % Transpose → 3×N, then concatenate
    end
end
