function export_variable_to_csv(inputFolder, variableName, outputFolder)
    
%%%
%Looks for multiple .mat files in a specified directory and extracts a given
%variable name from the .mat files and will save the variable under that
%given matfile
% Example usage 
%
% inputFolder = "C:\Users\user\Desktop\resp\Pig\Pig2\Signals";
% variableName = "rec.vest.Ve_filtered";
% outputFolder = "C:\Users\user\Desktop\resp\Pig\Pig2\SignalsCSV";
% export_variable_to_csv(inputFolder, variableName, outputFolder
%
%%%

    if nargin < 3
        outputFolder = fullfile(inputFolder, 'csv_output');
    end
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    matFiles = dir(fullfile(inputFolder, '**', '*.mat'));

    for i = 1:length(matFiles)
        matPath = fullfile(matFiles(i).folder, matFiles(i).name);
        try
            % Load only the top-level struct (first part of variableName)
            parts = strsplit(variableName, '.');
            topLevel = parts{1};
            data = load(matPath, topLevel);

            % Traverse nested fields
            variable = data;
            for p = 1:length(parts)
                if isstruct(variable) && isfield(variable, parts{p})
                    variable = variable.(parts{p});
                else
                    warning("Field '%s' not found in '%s'.", parts{p}, matPath);
                    variable = [];  % Mark invalid
                    break;
                end
            end

            if isempty(variable)
                continue;  % Skip if traversal failed
            end

            % Prepare output
            [~, baseName, ~] = fileparts(matFiles(i).name);
            relativePath = strrep(matFiles(i).folder, inputFolder, '');
            savePath = fullfile(outputFolder, relativePath);
            if ~exist(savePath, 'dir')
                mkdir(savePath);
            end
            csvFileName = fullfile(savePath, baseName + ".csv");

            % Write based on data type
            if isnumeric(variable) || islogical(variable)
                writematrix(variable, csvFileName);
            elseif ischar(variable) || isstring(variable)
                writecell(cellstr(variable), csvFileName);
            elseif isstruct(variable) && isfield(variable, 'name')
                writecell({variable.name}', csvFileName);
            elseif iscell(variable)
                writecell(variable, csvFileName);
            else
                warning("Unsupported data type in '%s'. Skipping.", matPath);
            end

        catch ME
            warning("Failed to process '%s': %s", matPath, ME.message);
        end
    end
end
