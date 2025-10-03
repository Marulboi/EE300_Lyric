figure;
% Create a tiled layout for the plots
tiledlayout(6, 2);

% Plot lambda_L_piecewise
nexttile;
plot(timeInterval, lambda_L_piecewise);
title('Lambda L Piecewise');
xlabel('Time Interval');
ylabel('Lambda L');

% Plot lambda_GCV_piecewise
nexttile;
plot(timeInterval, lambda_GCV_piecewise);
title('Lambda GCV Piecewise');
xlabel('Time Interval');
ylabel('Lambda GCV');

% Plot lambda_RGCV_piecewise
nexttile;
plot(timeInterval, lambda_RGCV_piecewise);
title('Lambda RGCV Piecewise');
xlabel('Time Interval');
ylabel('Lambda RGCV');

% Plot lambda_CRESO_piecewise
nexttile;
plot(timeInterval, lambda_CRESO_piecewise);
title('Lambda CRESO Piecewise');
xlabel('Time Interval');
ylabel('Lambda CRESO');

% Plot lambda_U_piecewise
nexttile;
plot(timeInterval, lambda_U_piecewise);
title('Lambda U Piecewise');
xlabel('Time Interval');
ylabel('Lambda U');

% Plot lambda_L_fix_piecewise
nexttile;
plot(timeInterval, lambda_L_fix_piecewise);
title('Lambda L Fix Piecewise');
xlabel('Time Interval');
ylabel('Lambda L Fix');

% Plot rho_piecewise
nexttile;
plot(timeInterval, rho_piecewise);
title('Rho Piecewise');
xlabel('Time Interval');
ylabel('Rho');

% Plot eta_piecewise
nexttile;
plot(timeInterval, eta_piecewise);
title('Eta Piecewise');
xlabel('Time Interval');
ylabel('Eta');

nexttile;
hold on;
plot(timeInterval, rec.vest.Ve_filtered(50, timeInterval), 'Color', [0.5, 0.5, 0.5, 0.2])
[env_up, ~] = envelope(rec.vest.Ve_filtered(5, timeInterval), 500, 'peak');
plot(timeInterval, env_up, 'r', 'LineWidth', 1.2);
xlabel('Time Index');
ylabel('Amplitude');
title('Respiration Signal + Envelope');

nexttile;
hold on;
plot(timeInterval, rec.vest.Ve_filtered(50, timeInterval), 'Color', [0.5, 0.5, 0.5, 0.2])
[env_up, ~] = envelope(rec.vest.Ve_filtered(5, timeInterval), 500, 'peak');
plot(timeInterval, env_up, 'r', 'LineWidth', 1.2);
xlabel('Time Index');
ylabel('Amplitude');
title('Respiration Signal + Envelope');


% Adjust layout
sgtitle('Piecewise Variables mean');