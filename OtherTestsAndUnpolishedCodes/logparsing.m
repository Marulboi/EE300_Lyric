filename = 'C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig5\SignalsCSV\png\infolog.log';  % Path to your log file

fid = fopen(filename, 'r');

pace = {};
top_values = {};
bottom_values = {};

current_pace = '';
reading_top = false;
reading_bottom = false;
top = [];
bottom = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));
    
    % Match line with file path
    if contains(line, 'File:') && contains(line, 'ECG lead selceted')
        if ~isempty(current_pace)
            pace{end+1} = current_pace;
            top_values{end+1} = top;
            bottom_values{end+1} = bottom;
        end
        top = [];
        bottom = [];
        reading_top = false;
        reading_bottom = false;

        % Extract pacing name from path (strip .csv and folders)
        tokens = regexp(line, 'File:\s*(.*?)([^\\/:]+)\.csv', 'tokens');
        if ~isempty(tokens)
            [~, baseName, ~] = fileparts(tokens{1}{2});
            current_pace = baseName;
        end

    elseif contains(line, 'rpeaks_top_20_amp')
        reading_top = true;
        reading_bottom = false;

    elseif contains(line, 'rpeaks_bottom_20_amp')
        reading_top = false;
        reading_bottom = true;

    elseif startsWith(line, 'INFO - ---------------------------------------------') || startsWith(line, 'INFO - Seleceted leads')
        reading_top = false;
        reading_bottom = false;

    else
        % Try to extract numeric values from rpeak lines
        tokens = regexp(line, '^\d+\s+(\d+)', 'tokens');
        if ~isempty(tokens)
            value = str2double(tokens{1}{1});
            if reading_top
                top(end+1) = value;
            elseif reading_bottom
                bottom(end+1) = value;
            end
        end
    end
end

% Append last block
if ~isempty(current_pace)
    pace{end+1} = current_pace;
    top_values{end+1} = top;
    bottom_values{end+1} = bottom;
end

fclose(fid);

%% Write to CSV
out_file = 'rpeak_summary.csv';
fid = fopen(out_file, 'w');
fprintf(fid, 'pace,top,bottom\n');
for i = 1:length(pace)
    top_str = strjoin(string(top_values{i}), ';');
    bottom_str = strjoin(string(bottom_values{i}), ';');
    fprintf(fid, '%s,"%s","%s"\n', pace{i}, top_str, bottom_str);
end
fclose(fid);

fprintf(' Data successfully written to %s\n', out_file);
