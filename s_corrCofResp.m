clc;

%close all;
clear;


%%

load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\ForwardM\ForwMat_Sock_HInTOut.mat','Trf_HT_coarse');
log_file = 'C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\SignalsCSV\png\infolog.log';

%%
[pace, ~, ~] = parse_resp_log(log_file);


%%
for i = 1:length(pace)
    path1 = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\Signals\',pace{i},'.mat'];
    
 
    disp(['Computing :',pace{i}]);
    
    % Derive the second path from the first 
    
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

    % Clear variables after use
    clear path1 folder name  filename2 path2;

%%
    [timeInterval,meanhelper] = timeintervalgenerator(Signal_select_Beat);
    
    
    %piecwise mean or median
    piecewise = true;
    
    N_time = length(timeInterval);
    
    nLeadsBody = size(rec.vest.Ve_filtered,1);
    brokenLeads = zeros(nLeadsBody,1);
    brokenLeads(rec.vest.BadChannels) = 1;
    
    
    
    % Define the lambda range
    lambda_range = logspace(-6, 0, 500);
    
    
    %apdc polynomial degree
    apdcDegree = 5 ;
    
    % Set robustness factor for RGCV
    robustness = 0.5;
    
    
    Asub = Trf_HT_coarse((brokenLeads == 0),:); 
    Sig_sub = rec.vest.Ve_filtered((brokenLeads == 0),timeInterval);
    %% Precompute SVD
    [U, S, V] = svd(Asub);
    S_diag = diag(S); 
    
    solsize = size(S,2);
    
    %% Initialize storage
    X_Lcurve = zeros(solsize, N_time);
    X_GCV = zeros(solsize, N_time);
    X_RGCV = zeros(solsize, N_time);
    X_CRESO = zeros(solsize, N_time);
    X_Ucurve = zeros(solsize, N_time);
    X_ADPC = zeros(solsize, N_time);
    X_my_Lcurve = zeros(solsize, N_time);
    
    % Lambdas at each time
    lambda_L = zeros(N_time,1);
    lambda_GCV = zeros(N_time,1);
    lambda_RGCV = zeros(N_time,1);
    lambda_CRESO = zeros(N_time,1);
    lambda_U = zeros(N_time,1);
    lambda_L_fix = zeros(N_time,1);
    
    
    cresomaxcurve = zeros(N_time,1);
    
    disp("initizalation complete")
    
    lambda_ADPC_mean = apdc(U, S_diag, Sig_sub,apdcDegree);  %Using my function as reg tools did not have it
    
    %% 1st Loop: Adaptive lambdas at each time step
    
    
    parfor t1 = 1:N_time
        b_t = Sig_sub(:,t1);
        %Compute lambdas for this b_t
        [lambda_L(t1), ~, ~] = l_curve(U, S_diag, b_t);  % Assuming l_curve returns (lambda, residual norm, solution norm)
    
        % Trying to fix L Curve with my function
    
        norm_svd_func = @(lambda) norms_svd(U, S_diag, b_t, lambda);
    
        [rho, eta] = norms_svd(U, S_diag, b_t, lambda_range);
        lambda_L_fix(t1) = my_l_curve(norm_svd_func, lambda_range);
        % [lambda_L_fix(t1),curve(t1)] = my_l_curve2(norm_svd_func, lambda_range);
        
        lambda_GCV(t1) = gcv(U, S_diag, b_t);
        lambda_RGCV(t1) = rgcv_svd(U, S, b_t, lambda_range, robustness); %Using my function as reg tools did not have it
        
        [lambda_CRESO(t1),cresomaxcurve(t1)] = creso(U, S_diag, b_t, lambda_range); %Using my function as reg tools did not have it
        
        %lambda_CRESO(t1) = creso(U, S_diag, b_t, lambda_range); %Using my function as reg tools did not have it
        
        lambda_U(t1) = newUcurve(U, S_diag, b_t, lambda_range); %Using my function as reg tools did not have it
        
        fprintf("lambdas found %d/%d \n", [t1,N_time]);
        
    end
    % Loop over beats
    %%
    if isempty(meanhelper)
        disp(['skipping:',pace{i}])
        continue
    end
    lambda_L_piecewise = zeros(1, meanhelper(end));
    lambda_GCV_piecewise = zeros(1, meanhelper(end));
    lambda_RGCV_piecewise = zeros(1, meanhelper(end));
    lambda_CRESO_piecewise = zeros(1, meanhelper(end));
    lambda_U_piecewise = zeros(1, meanhelper(end));
    lambda_L_fix_piecewise = zeros(1, meanhelper(end));
    for b = 1:meanhelper(end)
        % Compute median in this range
        lambda_L_m     = median(lambda_L(meanhelper == b));
        lambda_GCV_m   = median(lambda_GCV(meanhelper == b));
        lambda_RGCV_m  = median(lambda_RGCV(meanhelper == b));
        lambda_CRESO_m = median(lambda_CRESO(meanhelper == b));
        lambda_U_m     = median(lambda_U(meanhelper == b));
        lambda_L_fix_m = median(lambda_L_fix(meanhelper == b));
    
        % Assign the median value over that beat's time window
        lambda_L_piecewise(b)     = lambda_L_m;
        lambda_GCV_piecewise(b)   = lambda_GCV_m;
        lambda_RGCV_piecewise(b)  = lambda_RGCV_m;
        lambda_CRESO_piecewise(b) = lambda_CRESO_m;
        lambda_U_piecewise(b)     = lambda_U_m;
        lambda_L_fix_piecewise(b) = lambda_L_fix_m;
    end
    
    %%
    averageres = zeros(1, meanhelper(end));
    averagesol = zeros(1, meanhelper(end));
    k=1;
    for onset = Signal_select_Beat.QRS_Onset_auto' 
        
        % start = (i+Signal_select_Beat.BeatOffset);
        % stop = (i+Signal_select_Beat.BeatOffset+Signal_select_Beat.Beat_range);
        start = (onset);
        stop = (onset+Signal_select_Beat.QRS_range);
    
        Sig_sub = rec.vest.Ve_filtered((brokenLeads == 0),(start:stop));
        [sol_norm_avg, res_norm_avg, eta, rho] = avargeLcurve(U,S_diag,lambda_range,Sig_sub);
    
        averageres(k) = (res_norm_avg(end));
        averagesol(k) = (sol_norm_avg(end));
        k = k + 1 ;
    end
    %%
    data = [    
        num2cell(Signal_select_Beat.QRS_Onset_auto),...
        num2cell(lambda_L_piecewise'),...
        num2cell(lambda_GCV_piecewise'),...
        num2cell(lambda_RGCV_piecewise'),...
        num2cell(lambda_CRESO_piecewise'),...
        num2cell(lambda_U_piecewise'),...
        num2cell(lambda_L_fix_piecewise'),...
        num2cell(averageres'),...
        num2cell(averagesol'),...
    ];
    
    % Prepare the data for exporting
    headers = {'QRS Onset', 'Lambda L', 'Lambda GCV', 'Lambda RGCV', 'Lambda CRESO', 'Lambda U', 'Lambda L Fix', 'Average Residual', 'Average Solution'};
    data = [headers;data];  % Combine headers with the data
    
    % Check if the directory exists, if not, create it
    outputDir = fullfile(folder2, 'respcorcof');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    writetable(cell2table(data), fullfile(outputDir, [pace{i}, '_Respcheck.csv']));  % Save to CSV
end