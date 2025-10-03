clc;
%close all;
clear;

%%
PigNum = 5;
load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\ForwardM\ForwMat_Sock_HInTOut.mat','Trf_HT_coarse');
log_file = 'C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\SignalsCSV\png\infolog.log';

%%
[pace, top_vals, bot_vals] = parse_resp_log(log_file);


% Example: display first entry
disp(pace{1});
disp(top_vals{1});
disp(bot_vals{1});


csvdata = {'PigNum','Pacing Site',...
    'Lambda_L_inha','Lambda_L_exha','Lambda_L_rest',...
    'Lambda_GCV_inha','Lambda_GCV_exha','Lambda_GCV_rest',...
    'Lambda_RGCV_inha','Lambda_RGCV_exha','Lambda_RGCV_rest',...
    'Lambda_CRESO_inha','Lambda_CRESO_exha','Lambda_CRESO_rest',...
    'Lambda_U_inha','Lambda_U_exha','Lambda_U_rest',...
    'Lambda_L_fix_inha','Lambda_L_fix_exha','Lambda_L_fix_rest',...
    'Rho_inha','Rho_exha','Rho_rest',...
    'Eta_inha','Eta_exha','Eta_rest'};
%%
for i = 1:length(pace)
    path1 = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\Signals\',pace{i},'.mat'];
    
    if (length(top_vals{1})<2 || length(bot_vals{1})<2)
        disp(['Skipped :',pace{i}]);
        continue
    end   
    disp(['Computing :',pace{i}]);
    
    % Derive the second path from the first
    [folder, name, ~] = fileparts(path1);
    folder2 = strrep(folder, 'Signals', 'Signal_select_auto');
    filename2 = [name '_Signal_select_Beat.mat'];
    path2 = fullfile(folder2, filename2);
    load(path1);
    load(path2);
    
    % Clear variables after use
    clear path1 folder name folder2 filename2 path2;
    %%
    if(length(top_vals) > length(bot_vals))
        exhalationpeaks = top_vals{i};
        inhlationpeaks = bot_vals{i};
    
    else
        inhlationpeaks = top_vals{i};
        exhalationpeaks = bot_vals{i};
    end
    
    
    %%
    
    [timeInterval,~,respcycle] = timeintervalgenerator2(Signal_select_Beat,inhlationpeaks,exhalationpeaks);
    
    [~,meanhelper] = timeintervalgenerator(Signal_select_Beat);
    
    
    N_time = length(timeInterval);
    
    nLeadsBody = size(rec.vest.Ve_filtered,1);
    brokenLeads = zeros(nLeadsBody,1);
    brokenLeads(rec.vest.BadChannels) = 1;
    
    
    % Define the lambda range
    lambda_range = logspace(-6, 0, 500);
    
    
    % Set robustness factor for RGCV
    robustness = 0.5;
    
    
    Asub = Trf_HT_coarse((brokenLeads == 0),:); 
    Sig_sub = rec.vest.Ve_filtered((brokenLeads == 0),timeInterval);
    %% Precompute SVD
    [U, S, V] = svd(Asub);
    S_diag = diag(S); 
    
    solsize = size(S,2);
    
    % Lambdas at each time
    lambda_L = zeros(N_time,1);
    lambda_GCV = zeros(N_time,1);
    lambda_RGCV = zeros(N_time,1);
    lambda_CRESO = zeros(N_time,1);
    lambda_U = zeros(N_time,1);
    lambda_L_fix = zeros(N_time,1);
    
    
    disp("initizalation complete")
    
    
    %% 1st Loop: Adaptive lambdas at each time step
    
    rho = zeros(1,N_time); % Initialize solution norm array
    eta = zeros(1,N_time); % Initialize residual norm array
    
    
    parfor t1 = 1:N_time
        b_t = Sig_sub(:,t1);
        %Compute lambdas for this b_t
        [lambda_L(t1), ~, ~] = l_curve(U, S_diag, b_t);  % Assuming l_curve returns (lambda, residual norm, solution norm)
    
        % Trying to fix L Curve with my function
    
        norm_svd_func = @(lambda) norms_svd(U, S_diag, b_t, lambda);
    
        [r_tmp,e_tmp]  = norms_svd(U, S_diag, b_t, lambda_range);
        rho(t1) =  r_tmp(1);
        eta(t1) =  e_tmp(1);  
    
        lambda_L_fix(t1) = my_l_curve(norm_svd_func, lambda_range);
     
        lambda_GCV(t1) = gcv(U, S_diag, b_t);
        lambda_RGCV(t1) = rgcv_svd(U, S, b_t, lambda_range, robustness); %Using my function as reg tools did not have it
        
        [lambda_CRESO(t1),~] = creso(U, S_diag, b_t, lambda_range); %Using my function as reg tools did not have it
         
        lambda_U(t1) = newUcurve(U, S_diag, b_t, lambda_range); %Using my function as reg tools did not have it
        
        fprintf("lambdas found %d/%d \n", [t1,N_time]);
        
    end
    
    
    %%
    % Initialize the output vectors (same length as full signal)
    lambda_L_piecewise     = zeros(size(lambda_L));
    lambda_GCV_piecewise   = zeros(size(lambda_GCV));
    lambda_RGCV_piecewise  = zeros(size(lambda_RGCV));
    lambda_CRESO_piecewise = zeros(size(lambda_CRESO));
    lambda_U_piecewise     = zeros(size(lambda_U));
    lambda_L_fix_piecewise = zeros(size(lambda_L_fix));
    
    % Compute median in this range for rho and eta
    rho_piecewise = zeros(size(rho));
    eta_piecewise = zeros(size(eta));
        
    
    % Loop over beats
    for b = 1:meanhelper(end)
        % Compute median in this range
        % Assign the median value over that beat's time window directly
        lambda_L_piecewise(meanhelper == b)     = mean(lambda_L(meanhelper == b));
        lambda_GCV_piecewise(meanhelper == b)   = mean(lambda_GCV(meanhelper == b));
        lambda_RGCV_piecewise(meanhelper == b)  = mean(lambda_RGCV(meanhelper == b));
        lambda_CRESO_piecewise(meanhelper == b) = mean(lambda_CRESO(meanhelper == b));
        lambda_U_piecewise(meanhelper == b)     = mean(lambda_U(meanhelper == b));
        lambda_L_fix_piecewise(meanhelper == b) = mean(lambda_L_fix(meanhelper == b));
    
        rho_piecewise(meanhelper == b) = mean(rho(meanhelper == b));
        eta_piecewise(meanhelper == b) = mean(eta(meanhelper == b));
    end
    %%
    j = 1;
    % Initialize arrays for average lambdas and norms
    lambda_L_avgbeat = zeros(1, meanhelper(end));
    lambda_GCV_avgbeat = zeros(1, meanhelper(end));
    lambda_RGCV_avgbeat = zeros(1, meanhelper(end));
    lambda_CRESO_avgbeat = zeros(1, meanhelper(end));
    lambda_U_avgbeat = zeros(1, meanhelper(end));
    lambda_L_fix_avgbeat = zeros(1, meanhelper(end));
    rho_avgbeat = zeros(1, meanhelper(end));
    eta_avgbeat = zeros(1, meanhelper(end));
    
    for k = Signal_select_Beat.QRS_Onset_auto'
        lambda_L_avgbeat(j) = lambda_L_piecewise(find(timeInterval == k + 2,1));
        lambda_GCV_avgbeat(j) = lambda_GCV_piecewise(find(timeInterval == k + 2,1));
        lambda_RGCV_avgbeat(j) = lambda_RGCV_piecewise(find(timeInterval == k + 2,1));
        lambda_CRESO_avgbeat(j) = lambda_CRESO_piecewise(find(timeInterval == k + 2,1));
        lambda_U_avgbeat(j) = lambda_U_piecewise(find(timeInterval == k + 2,1));
        lambda_L_fix_avgbeat(j) = lambda_L_fix_piecewise(find(timeInterval == k + 2,1));
        rho_avgbeat(j) = rho_piecewise(find(timeInterval == k + 2,1));
        eta_avgbeat(j) = eta_piecewise(find(timeInterval == k + 2,1));
        j = j + 1; % Increment the index for the next beat
    end
    
    
    
    %%
    groupstr = @(vec, x) strjoin(string(vec(respcycle == x)), ';');
    
    appendData = {...
        PigNum, pace{i}, ...
        groupstr(lambda_L_avgbeat, 1), groupstr(lambda_L_avgbeat, 2), groupstr(lambda_L_avgbeat, 3), ...
        groupstr(lambda_GCV_avgbeat, 1), groupstr(lambda_GCV_avgbeat, 2), groupstr(lambda_GCV_avgbeat, 3), ...
        groupstr(lambda_RGCV_avgbeat, 1), groupstr(lambda_RGCV_avgbeat, 2), groupstr(lambda_RGCV_avgbeat, 3), ...
        groupstr(lambda_CRESO_avgbeat, 1), groupstr(lambda_CRESO_avgbeat, 2), groupstr(lambda_CRESO_avgbeat, 3), ...
        groupstr(lambda_U_avgbeat, 1), groupstr(lambda_U_avgbeat, 2), groupstr(lambda_U_avgbeat, 3), ...
        groupstr(lambda_L_fix_avgbeat, 1), groupstr(lambda_L_fix_avgbeat, 2), groupstr(lambda_L_fix_avgbeat, 3), ...
        groupstr(rho_avgbeat, 1), groupstr(rho_avgbeat, 2), groupstr(rho_avgbeat, 3), ...
        groupstr(eta_avgbeat, 1), groupstr(eta_avgbeat, 2), groupstr(eta_avgbeat, 3) ...
    };
    
    csvdata = [csvdata; appendData];
end


% Save the results to a CSV file
writecell(csvdata, 'BoxPlotResults.csv');

%% === Convert in-memory `csvdata` into long format for Prism ===
% Skip header
header = raw(1, :);
data = raw(2:end, :);

% The first two columns are PigNum and PacingSite
value_columns = header(3:end);
n_cols = length(value_columns);

long_data = {}; % {PigNum, PacingSite, Condition, Value}

for i = 1:size(data, 1)
    pig = data{i, 1};
    pacing_site = data{i, 2};
    
    for c = 1:n_cols
        label = value_columns{c};
        value_str = data{i, c + 2};  % offset by PigNum and Pacing Site
        
        % Skip if empty
        if isempty(value_str)
            continue;
        end

        % Ensure it's treated as a string before splitting
        value_str = string(value_str);  

        % Split and convert
        values = str2double(split(value_str, ';'));
        
        % Append each value as a new row
        for v = 1:numel(values)
            long_data(end + 1, :) = {pig, pacing_site, label, values(v)}; %#ok<AGROW>
        end
    end
end

% Convert to table and write
Tlong = cell2table(long_data, 'VariableNames', {'PigNum', 'PacingSite', 'Condition', 'Value'});
writetable(Tlong, 'BoxPlotResults_LongFormat_Prism.csv');

disp('Wrote long-format Prism CSV: BoxPlotResults_LongFormat_Prism.csv');