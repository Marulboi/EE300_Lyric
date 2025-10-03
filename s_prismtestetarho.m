clc;
clear;

%% === INITIALIZATION ===
pig_num = num2str(5);
load(['C:\Users\emire\OneDrive\Desktop\resp\Pig\Pig',pig_num,'\ForwardM\ForwMat_Sock_HInTOut.mat'],'Trf_HT_coarse');
log_file = ['C:\Users\emire\OneDrive\Desktop\resp\Pig\Pig',pig_num,'\SignalsCSV\pngN\infolog.log'];
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
outdir = ['ExportedCSVsNEtarho',pig_num];

mkdir(outdir);


%% === MAIN LOOP ===
for i = 1:length(pace)
    path1 = ['C:\Users\emire\OneDrive\Desktop\resp\Pig\Pig',pig_num,'\Signals\', pace{i}, '.mat'];
    if length(top_vals{i}) < 2 || length(bot_vals{i}) < 2
        disp(['Skipped :', pace{i}]);
        continue
    end
    disp(['Computing :', pace{i}]);

    [folder, name, ~] = fileparts(path1);
    folder2 = strrep(folder, 'Signals', 'Signal_select_autoPan');
    filename2 = [name '.mat'];
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


    rho = zeros(1, N_time);
    eta = zeros(1, N_time);

    parfor t1 = 1:N_time
        b_t = Sig_sub(:, t1);
        
        norm_svd_func = @(lambda) norms_svd(U, S_diag, b_t, lambda);
        lambda_L_fix(t1) = my_l_curve(norm_svd_func, lambda_range);
        [r_tmp, e_tmp] = norms_svd(U, S_diag, b_t, lambda_range);
        lambdaid = find(lambda_range==lambda_L_fix(t1))  
        rho(t1) = r_tmp(lambdaid);
        eta(t1) = e_tmp(lambdaid);
        
    end


    rho_piecewise          = zeros(size(rho));
    eta_piecewise          = zeros(size(eta));
    
    if isempty(meanhelper)
        disp(['skipping:',paceCHU{i}])
        continue
    end
    for b = 1:meanhelper(end)

        rho_piecewise(meanhelper == b) = mean(rho(meanhelper == b));
        eta_piecewise(meanhelper == b) = mean(eta(meanhelper == b));
    end


    current_rho = [];
    current_eta = [];
    
    QRS_Onset_auto = Signal_select_Beat.QRS_Onset_auto';
    QRS_Onset_auto(Signal_select_Beat.QRS_Onset_auto'<0)=1;

    for k = QRS_Onset_auto
        idx = find(timeInterval == k + 2,1);

        current_rho(end+1)          = rho_piecewise(idx);
        current_eta(end+1)          = eta_piecewise(idx);
    end

    rho_data          = padcat(rho_data,          padcat(current_rho(respcycle==1)', current_rho(respcycle==2)', current_rho(respcycle==3)'));
    eta_data          = padcat(eta_data,          padcat(current_eta(respcycle==1)', current_eta(respcycle==2)', current_eta(respcycle==3)'));
    
    disp('==================')
    disp(eta_data)
    disp('==================')
    
    column_labels = [column_labels, pace{i}];
end

%% Save outputs to folder

rho_data_T = concat_blocks_by_3(rho_data);
eta_data_T = concat_blocks_by_3(eta_data);

toplabel = repelem(column_labels, size(eta_data, 1));


rho_data_T = [toplabel; num2cell(rho_data_T)];
rho_data_T = [leftlabel, rho_data_T];

eta_data_T = [toplabel; num2cell(eta_data_T)];
eta_data_T = [leftlabel, eta_data_T];
%%


writecell(rho_data_T,          fullfile(outdir, 'LcurveRho.csv'));

writecell(eta_data_T,          fullfile(outdir, 'LcurveEta.csv'));

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
