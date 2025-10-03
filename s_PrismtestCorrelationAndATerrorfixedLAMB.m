clc;
clear; 

%% === INITIALIZATION ===
pig_num = num2str(5);
load(['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig',pig_num,'\ForwardM\ForwMat_Sock_HInTOut.mat'],'Trf_HT_coarse');
log_file = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig',pig_num,'\SignalsCSV\pngN\infolog.log'];
[pace, top_vals, bot_vals,rest_vals] = parse_resp_log(log_file);

% Output matrices for each quantity
spaitelCC = [];
spaitelCC_F = [];
atError = [];
column_labels = {};


leftlabel = {' ';'Inhalation';'Exhalation';'rest'};
%%
% Output directory
outdir = ['ExportedCSVsCCmedianfixedlambda',pig_num];

mkdir(outdir);



%% === MAIN LOOP ===
for i = 1:length(pace)
    path1 = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig',pig_num,'\Signals\', pace{i}, '.mat'];
    if length(top_vals{i}) < 2 || length(bot_vals{i}) < 2
        disp(['Skipped :', pace{i}]);
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

    [U, S, V] = svd(Asub);
    S_diag = diag(S);
    solsize = size(S,2);

    lambda_L_fix = zeros(N_time, 1);

    real_values = rec.sock.Ve_filtered(:, timeInterval);
    crop = solsize-size(real_values,1);
    sizeof = size(real_values,1);
    interval = 1:sizeof;

    
    lambda_L_fix_piecewise = [];

    parfor t1 = 1:N_time

        b_t = Sig_sub(:, t1);
        norm_svd_func = @(lambda) norms_svd(U, S_diag, b_t, lambda);
        lambda_L_fix(t1) = my_l_curve(norm_svd_func, lambda_range);

    end

    
    if isempty(meanhelper)
        disp(['skipping:',pace{i}])
        continue
    end
  

     lambda_L_fix_piecewise = mean(lambda_L_fix);



    parfor t3 = 1:N_time
        b_t = Sig_sub(:,t3);
        
        X_my_Lcurve(:,t3) = tikhonov(U, S_diag, V, b_t, lambda_L_fix_piecewise);
        
        fprintf("Solutions mean %d // %d \n", [t3,N_time]);
    end

    X_my_Lcurve = medfilt2(X_my_Lcurve);


    current_lambda_L_fix = [];
    
    spaitelCC = [];
    for t = 1:meanhelper(end)
          corr_adaptive = zeros(1, meanhelper(end));
        for tt = 1:sizeof
            x_true = real_values(tt, (meanhelper==t));  % True signal at time t
            x_est_adapt = X_my_Lcurve(tt,(meanhelper==t));
    
            r_adapt = corrcoef(x_est_adapt, x_true);
    
            corr_adaptive(tt) = r_adapt(1,2);
    
    
           
        end
         fprintf("Error calc %d // %d \n", [t,N_time]);
        spaitelCC(t) = median(corr_adaptive(~isnan(corr_adaptive)));

    end

    % for b = 1:meanhelper(end)
    % 
    %     r_adapt_piecewise(meanhelper == b) = median(corr_adaptive(meanhelper == b));
    % 
    % end
    % 
    % spaitelCC = [];
    % for k = Signal_select_Beat.QRS_Onset_auto'
    %     idx = find(timeInterval == k + 2,1);
    % 
    %     spaitelCC(end+1) = r_adapt_piecewise(idx);
    % 
    % end
    if isempty(spaitelCC)
        disp('Skipping due to full NaN in spaitelCC for respcycle 1, 2, or 3');
        continue;
    end 

    spaitelCC_F = padcat(spaitelCC_F, padcat(spaitelCC(respcycle==1)', spaitelCC(respcycle==2)', spaitelCC(respcycle==3)'));

    disp('==================')
    disp(spaitelCC_F  )
    disp('==================')
    
    column_labels = [column_labels, pace{i}];


    % plot(spaitelCC);
    % pause;
end
%% Save outputs to folder
spaitelCC_T = concat_blocks_by_3(spaitelCC_F );


toplabel = repelem(column_labels, size(spaitelCC_F  , 1));

spaitelCC_T = [toplabel; num2cell(spaitelCC_T)];
spaitelCC_T = [leftlabel, spaitelCC_T];

%%

writecell(spaitelCC_T,     fullfile(outdir, 'temporal.csv'));


disp(['All files saved successfully to',outdir]);
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
