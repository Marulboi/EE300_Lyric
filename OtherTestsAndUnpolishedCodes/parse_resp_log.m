function [pace_names, top_peaks, bottom_peaks,rest_peaks] = parse_resp_log(log_path)
%PARSE_RPEAK_LOG Parse a log file to extract pacing site names and top/bottom rpeak values.
%
%   [pace_names, top_peaks, bottom_peaks] = parse_rpeak_log(log_path)
%
%   Inputs:
%       log_path      - Path to the .log file
%
%   Outputs:
%       pace_names    - Cell array of pacing site names
%       top_peaks     - Cell array of numeric arrays (top rpeak amplitudes)
%       bottom_peaks  - Cell array of numeric arrays (bottom rpeak amplitudes)

    fid = fopen(log_path, 'r');
    if fid == -1
        error('Could not open log file: %s', log_path);
    end

    pace_names = {};
    top_peaks = {};
    bottom_peaks = {};

    current_pace = '';
    reading_top = false;
    reading_bottom = false;
    top = [];
    bottom = [];
    
    rest_peaks = {};
    reading_outside = false;
    rest = [];

    while ~feof(fid)
        line = strtrim(fgetl(fid));

        % Match line with file path
        if contains(line, 'File:') && contains(line, 'ECG lead selceted')
            if ~isempty(current_pace)
                pace_names{end+1} = current_pace;
                top_peaks{end+1} = top;
                bottom_peaks{end+1} = bottom; rest_peaks{end+1} = rest;
            end
            top = [];
            bottom = [];
            rest = [];
            reading_top = false;
            reading_bottom = false;
            reading_outside = false;

            % Extract pacing name
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
        % after the existing "elseif contains(line, 'rpeaks_bottom_20_amp')" block
        elseif contains(line, 'rpeaks_outside')
            reading_top = false;
            reading_bottom = false;
            reading_outside = true;   % <-- ADD
        
        % extend your reset condition to also stop reading_outside
        elseif startsWith(line, 'INFO - ---------------------------------------------') || ...
               startsWith(line, 'INFO - Seleceted leads')
            reading_top = false;
            reading_bottom = false;
            reading_outside = false;   % <-- ADD


        else
            % Try to extract numeric value
            tokens = regexp(line, '^\d+\s+(\d+)', 'tokens');
            if ~isempty(tokens)
                value = str2double(tokens{1}{1});
                if reading_top
                    top(end+1) = value;
                elseif reading_bottom
                    bottom(end+1) = value;
                elseif reading_outside                 % <-- ADD
                    rest(end+1) = value;               % <-- ADD
                end
            end
        end
    end

    % Append final entry
    if ~isempty(current_pace)
        pace_names{end+1} = current_pace;
        top_peaks{end+1} = top;
        bottom_peaks{end+1} = bottom;
        rest_peaks{end+1} = rest;
    end

    fclose(fid);
end
