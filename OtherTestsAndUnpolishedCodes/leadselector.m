load('C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig2\Signals\EpiPacingSite-2.mat');
realvalues = rec.sock.Ve_filtered;  % Corrected variable name

% Initialize node lists
good_nodes = [];
bad_nodes = [];

n_nodes = size(realvalues, 1);

for i = 1:n_nodes
    figure(1); clf;
    
    % Plot original signal
    signal = realvalues(i, :);


    plot(signal, 'b'); hold on;
    
    % Calculate and plot envelopes
    try
        [env_upper, env_lower] = envelope(signal, 700, 'peak');
        plot(env_upper, 'r--', 'LineWidth', 1.2);  % Top envelope
        plot(env_lower, 'g--', 'LineWidth', 1.2);  % Bottom envelope
        legend('Signal', 'Upper Envelope', 'Lower Envelope');
    catch
        warning('Could not compute envelope for node %d (check signal length)', i);
    end

    title(sprintf('Node %d: Press ''a'' (good), ''r'' (bad), any other to skip', i));
    xlabel('Time');
    ylabel('Signal');

    % Get user input
    user_input = input('Input: ', 's');
    
    if strcmpi(user_input, 'a')
        good_nodes(end+1) = i;
    elseif strcmpi(user_input, 'r')
        bad_nodes(end+1) = i;
    end
end

% Display results
fprintf('\nGood nodes (a):\n');
disp(good_nodes);

fprintf('Bad nodes (r):\n');
disp(bad_nodes);