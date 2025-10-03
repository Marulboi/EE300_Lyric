clc;
clear;

%% Configuration

load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig1\ForwardM\ForwMat_Sock_HInTOut.mat','Trf_HT_coarse');
log_file = 'C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig1\SignalsCSV\png\infolog.log';
[pace, ~, ~] = parse_resp_log(log_file);

%% Main Loop
for i = 1:length(pace)
    path1 = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig1\Signals\',pace{i},'.mat'];

    disp(['Computing :',pace{i}]);

    [folder, name, ~] = fileparts(path1);
    folder2 = strrep(folder, 'Signals', 'Signal_select_auto');
    filename2 = [name '_Signal_select_Beat.mat'];
    path2 = fullfile(folder2, filename2);

    if exist(path1, 'file') && exist(path2, 'file')
        load(path1);
        load(path2);
    else
        fprintf('Skipping: %s or %s not found.\n', path1, path2);
        continue;
    end

    clear path1 folder name filename2 path2;

    [timeInterval, meanhelper] = timeintervalgenerator(Signal_select_Beat);

    piecewise = true;
    N_time = length(timeInterval);

    nLeadsBody = size(rec.vest.Ve_filtered, 1);
    brokenLeads = zeros(nLeadsBody, 1);
    brokenLeads(rec.vest.BadChannels) = 1;

    lambda_range = logspace(-6, 0, 500);
    apdcDegree = 5;
    robustness = 0.5;

    Asub = Trf_HT_coarse((brokenLeads == 0), :); 
    Sig_sub = rec.vest.Ve_filtered((brokenLeads == 0), timeInterval);

    
    filteredECG = rec.vest.Ve_filtered((brokenLeads == 0), :);
    
    Fs = 2048;                                     % Sampling frequency                    
    T = 1 / Fs;                                    % Sampling period       
    L = size(rec.vest.Ve_filtered, 2);             % Length of signal
    t = (0:L-1) * T;                               % Time vector
    
    % FFT of a selected ECG lead (e.g., lead 50)
    fftfiltECG = fft(filteredECG(50, :));
    
    % Full frequency mask (same size as fftfiltECG)
    full_mask = zeros(1, L);
    cutoff = 0.4;                                  % Cutoff frequency in Hz
    cutoff_idx = floor(cutoff * L / Fs);          % Index corresponding to 0.4 Hz
    
    % Preserve symmetric FFT by keeping both positive and negative frequencies
    full_mask(1:cutoff_idx) = 1;
    full_mask(end - cutoff_idx + 2:end) = 1;
    
    % Apply mask and perform inverse FFT to get low-pass filtered signal
    filteredECG_lp = real(ifft(fftfiltECG .* full_mask));


    [U, S, V] = svd(Asub);
    S_diag = diag(S); 
    solsize = size(S, 2);

    % Initialize storage
    X_Lcurve = zeros(solsize, N_time);
    X_GCV = zeros(solsize, N_time);
    X_RGCV = zeros(solsize, N_time);
    X_CRESO = zeros(solsize, N_time);
    X_Ucurve = zeros(solsize, N_time);
    X_ADPC = zeros(solsize, N_time);
    X_my_Lcurve = zeros(solsize, N_time);

    lambda_L = zeros(N_time, 1);
    lambda_GCV = zeros(N_time, 1);
    lambda_RGCV = zeros(N_time, 1);
    lambda_CRESO = zeros(N_time, 1);
    lambda_U = zeros(N_time, 1);
    lambda_L_fix = zeros(N_time, 1);

    cresomaxcurve = zeros(N_time, 1);

    disp("Initialization complete");

    lambda_ADPC_mean = apdc(U, S_diag, Sig_sub, apdcDegree);

    %% Adaptive Lambda Computation
    parfor t1 = 1:N_time
        b_t = Sig_sub(:, t1);

        [lambda_L(t1), ~, ~] = l_curve(U, S_diag, b_t);
        norm_svd_func = @(lambda) norms_svd(U, S_diag, b_t, lambda);
        [rho, eta] = norms_svd(U, S_diag, b_t, lambda_range);
        lambda_L_fix(t1) = my_l_curve(norm_svd_func, lambda_range);
        lambda_GCV(t1) = gcv(U, S_diag, b_t);
        lambda_RGCV(t1) = rgcv_svd(U, S, b_t, lambda_range, robustness);
        [lambda_CRESO(t1), cresomaxcurve(t1)] = creso(U, S_diag, b_t, lambda_range);
        lambda_U(t1) = newUcurve(U, S_diag, b_t, lambda_range);

        fprintf("Lambdas found %d/%d \n", [t1, N_time]);
    end

    if isempty(meanhelper)
        disp(['Skipping:', pace{i}])
        continue
    end

    %% Piecewise Lambda Averaging
    lambda_L_piecewise     = zeros(1, meanhelper(end));
    lambda_GCV_piecewise   = zeros(1, meanhelper(end));
    lambda_RGCV_piecewise  = zeros(1, meanhelper(end));
    lambda_CRESO_piecewise = zeros(1, meanhelper(end));
    lambda_U_piecewise     = zeros(1, meanhelper(end));
    lambda_L_fix_piecewise = zeros(1, meanhelper(end));

    for b = 1:meanhelper(end)
        lambda_L_piecewise(b)     = median(lambda_L(meanhelper == b));
        lambda_GCV_piecewise(b)   = median(lambda_GCV(meanhelper == b));
        lambda_RGCV_piecewise(b)  = median(lambda_RGCV(meanhelper == b));
        lambda_CRESO_piecewise(b) = median(lambda_CRESO(meanhelper == b));
        lambda_U_piecewise(b)     = median(lambda_U(meanhelper == b));
        lambda_L_fix_piecewise(b) = median(lambda_L_fix(meanhelper == b));
    end

    %% L-Curve Average Norms per Beat
    averageres = zeros(1, meanhelper(end));
    averagesol = zeros(1, meanhelper(end));
    k = 1;

    for onset = Signal_select_Beat.QRS_Onset_auto'
        start = onset;
        stop = onset + Signal_select_Beat.QRS_range;
        Sig_sub = rec.vest.Ve_filtered((brokenLeads == 0), start:stop);

        [sol_norm_avg, res_norm_avg, eta, rho] = avargeLcurve(U, S_diag, lambda_range, Sig_sub);
        averageres(k) = res_norm_avg(end);
        averagesol(k) = sol_norm_avg(end);
        k = k + 1;
    end


    % [~, locs] = findpeaks(rec.vest.Ve_filtered(50, timeInterval));
    % 
    % averageres = zeros(1, meanhelper(end));
    % averagesol = zeros(1, meanhelper(end));
    % k = 1;
    % 
    % for onset = locs
    %     b_t = rec.vest.Ve_filtered((brokenLeads == 0), onset);
    % 
    %     [ rho, eta] = norms_svd(U, S_diag, b_t, lambda_range);
    %     averageres(k) = rho(end);
    %     averagesol(k) = eta(end);
    %     k = k + 1;
    % end

    %% At this point:
    % - You have all computed data in memory:
    %   lambda_L_piecewise, lambda_GCV_piecewise, etc.
    %   averageres, averagesol
    % - You can optionally plot or analyze them.

    resppoints = filteredECG_lp(Signal_select_Beat.QRS_Onset_auto);
    
    
           % Collect all data rows to compare
    all_data = [
        lambda_L_piecewise;
        lambda_GCV_piecewise;
        lambda_RGCV_piecewise;
        lambda_CRESO_piecewise;
        lambda_U_piecewise;
        lambda_L_fix_piecewise;
        averageres;
        averagesol
    ];
    
    % Labels for legend (must match rows in all_data)
    labels = {
        'L-curve', ...
        'GCV', ...
        'RGCV', ...
        'CRESO', ...
        'U-curve', ...
        'L-curve (fixed)', ...
        'Avg Residual', ...
        'Avg Solution'
    };
    
    % Start figure
    figure; hold on;
    
    % Plot the reference ECG signal
    plot(resppoints, 'k', 'LineWidth', 2, 'DisplayName', 'Filtered ECG_{LP}');
    
    % Loop through each data row and plot with correlation
    n = size(all_data, 1);
    for l = 1:n
        data_row = all_data(l, :);
        R = corrcoef(data_row, resppoints);
        corrval = R(1, 2);  % Pearson correlation coefficient
    
        plot(data_row, 'DisplayName', sprintf('%s (r = %.2f)', labels{l}, corrval));
    end
    
    % Plot setup
    legend show;
    xlabel('Beat Index');
    ylabel('Value');
    title(['Correlation with Filtered ECG LP',pace{i}]);
    grid on;
        


end
