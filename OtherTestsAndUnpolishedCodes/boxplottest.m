clc;
clear;

%% === INITIALIZATION ===
PigNum = 5;
load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\ForwardM\ForwMat_Sock_HInTOut.mat','Trf_HT_coarse');
log_file = 'C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\SignalsCSV\png\infolog.log';
[pace, top_vals, bot_vals] = parse_resp_log(log_file);

% Output matrices for each quantity
lambda_L_data = [NaN];
lambda_GCV_data = [NaN];
lambda_RGCV_data = [NaN];
lambda_CRESO_data = [NaN];
lambda_U_data = [NaN];
lambda_L_fix_data = [NaN];
rho_data = [NaN];
eta_data = [NaN];
column_labels = {'test'};

% Output directory
outdir = 'ExportedCSVs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% === MAIN LOOP ===
for i = 1:length(pace)
    path1 = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\Signals\', pace{i}, '.mat'];
    if length(top_vals{1}) < 2 || length(bot_vals{1}) < 2
        disp(['Skipped :', pace{i}]);
        continue
    end
    disp(['Computing :', pace{i}]);

    [folder, name, ~] = fileparts(path1);
    folder2 = strrep(folder, 'Signals', 'Signal_select_auto');
    filename2 = [name '_Signal_select_Beat.mat'];
    path2 = fullfile(folder2, filename2);
    load(path1); load(path2);

    if length(top_vals) > length(bot_vals)
        exhalationpeaks = top_vals{i};
        inhlationpeaks = bot_vals{i};
    else
        inhlationpeaks = top_vals{i};
        exhalationpeaks = bot_vals{i};
    end

    [timeInterval, ~, respcycle] = timeintervalgenerator2(Signal_select_Beat, inhlationpeaks, exhalationpeaks);
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
        rho(t1) = r_tmp(1);
        eta(t1) = e_tmp(1);
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

    j = 1;
    current_lambda_L = [];
    current_lambda_GCV = [];
    current_lambda_RGCV = [];
    current_lambda_CRESO = [];
    current_lambda_U = [];
    current_lambda_L_fix = [];
    current_rho = [];
    current_eta = [];

    for k = Signal_select_Beat.QRS_Onset_auto'
        current_lambda_L(end+1)     = lambda_L_piecewise(find(timeInterval == k + 2,1));
        current_lambda_GCV(end+1)   = lambda_GCV_piecewise(find(timeInterval == k + 2,1));
        current_lambda_RGCV(end+1)  = lambda_RGCV_piecewise(find(timeInterval == k + 2,1));
        current_lambda_CRESO(end+1) = lambda_CRESO_piecewise(find(timeInterval == k + 2,1));
        current_lambda_U(end+1)     = lambda_U_piecewise(find(timeInterval == k + 2,1));
        current_lambda_L_fix(end+1) = lambda_L_fix_piecewise(find(timeInterval == k + 2,1));
        current_rho(end+1)          = rho_piecewise(find(timeInterval == k + 2,1));
        current_eta(end+1)          = eta_piecewise(find(timeInterval == k + 2,1));
        j = j + 1;
    end
lambda_L_data     = padcat(lambda_L_data,     padcat(num2cell(current_lambda_L(respcycle==1)'), num2cell(current_lambda_L(respcycle==2)'), num2cell(current_lambda_L(respcycle==3)')));
lambda_GCV_data   = padcat(lambda_GCV_data,   padcat(num2cell(current_lambda_GCV(respcycle==1)'), num2cell(current_lambda_GCV(respcycle==2)'), num2cell(current_lambda_GCV(respcycle==3)')));
lambda_RGCV_data  = padcat(lambda_RGCV_data,  padcat(num2cell(current_lambda_RGCV(respcycle==1)'), num2cell(current_lambda_RGCV(respcycle==2)'), num2cell(current_lambda_RGCV(respcycle==3)')));
lambda_CRESO_data = padcat(lambda_CRESO_data, padcat(num2cell(current_lambda_CRESO(respcycle==1)'), num2cell(current_lambda_CRESO(respcycle==2)'), num2cell(current_lambda_CRESO(respcycle==3)')));
lambda_U_data     = padcat(lambda_U_data,     padcat(num2cell(current_lambda_U(respcycle==1)'), num2cell(current_lambda_U(respcycle==2)'), num2cell(current_lambda_U(respcycle==3)')));
lambda_L_fix_data = padcat(lambda_L_fix_data, padcat(num2cell(current_lambda_L_fix(respcycle==1)'), num2cell(current_lambda_L_fix(respcycle==2)'), num2cell(current_lambda_L_fix(respcycle==3)')));
rho_data          = padcat(rho_data,          padcat(num2cell(current_rho(respcycle==1)'), num2cell(current_rho(respcycle==2)'), num2cell(current_rho(respcycle==3)')));
eta_data          = padcat(eta_data,          padcat(num2cell(current_eta(respcycle==1)'), num2cell(current_eta(respcycle==2)'), num2cell(current_eta(respcycle==3)')));


    column_labels = [column_labels, {[pace{i} ":1"], [pace{i} ":2"], [pace{i} ":3"]}];
end

% Save outputs to folder
writecell([column_labels; lambda_L_data],     fullfile(outdir, 'Lambda_L.csv'));
writecell([column_labels; lambda_GCV_data],   fullfile(outdir, 'Lambda_GCV.csv'));
writecell([column_labels; lambda_RGCV_data],  fullfile(outdir, 'Lambda_RGCV.csv'));
writecell([column_labels; lambda_CRESO_data], fullfile(outdir, 'Lambda_CRESO.csv'));
writecell([column_labels; lambda_U_data],     fullfile(outdir, 'Lambda_U.csv'));
writecell([column_labels; lambda_L_fix_data], fullfile(outdir, 'Lambda_L_fix.csv'));
writecell([column_labels; rho_data],          fullfile(outdir, 'Rho.csv'));
writecell([column_labels; eta_data],          fullfile(outdir, 'Eta.csv'));


disp('All files saved successfully to ExportedCSVs/');
function C = padcat(varargin)
% PADCAT Pads and concatenates vectors and matrices column-wise
%   C = padcat(A, B, C, ...) pads inputs with NaNs to match the
%   maximum number of rows and concatenates them horizontally.

    narginchk(1, Inf);

    % Get number of rows for each input
    numRows = zeros(1, nargin);
    for k = 1:nargin
        v = varargin{k};
        if isvector(v)
            numRows(k) = numel(v);
        else
            numRows(k) = size(v, 1);
        end
    end

    maxRows = max(numRows);

    % Concatenate with padding
    C = [];
    for k = 1:nargin
        v = varargin{k};
        if isvector(v)
            v = v(:);  % make it column
        end

        [r, c] = size(v);
        if r < maxRows
            v = [v; num2cell(NaN(maxRows - r, c))];
        end
        if k == 1
            C = v;
        
        else
            C = [C, v];
        end
    end
end
%%
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
        mat = cell2mat(block);                % N×3 numeric matrix
        out = [out, mat.'];                   % Transpose → 3×N, then concatenate
    end
end

