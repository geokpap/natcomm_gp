% List of the MAT files
files = {'ApAv_OLDfunction.mat','ApAv_em_1500ms.mat', 'ApAvApAp_20241114v2_1500ms.mat'};

% Specify the variables to combine
variablesToCombine = {'List2', 'stat_1_v2', 'Firing_rates_f', 'allPerTrialData', 'allSpikesPerBin', 'results','fromExcel1'};

% Initialize a structure to hold the combined data for these variables
combinedData = struct();
for varIndex = 1:length(variablesToCombine)
    combinedData.(variablesToCombine{varIndex}) = [];
end

% Initialize to track the max dimensions for variables that need padding
maxLayersStat1_v2 = 0;
maxColsFiringRates_f = 0;
maxColsAllPerTrialData = 0;
maxColsAllSpikesPerBin = 0;
maxColsResults = 0;  % Initialize if results is a dimensioned variable

% Loop through each file
for k = 1:length(files)
    % Load data from each MAT file using variable arguments
    data = load(files{k}, variablesToCombine{:});

    % Update the max dimensions for specific variables if necessary
    if isfield(data, 'stat_1_v2')
        maxLayersStat1_v2 = max(maxLayersStat1_v2, size(data.stat_1_v2, 3));
    end
    if isfield(data, 'Firing_rates_f')
        maxColsFiringRates_f = max(maxColsFiringRates_f, size(data.Firing_rates_f, 2));
    end
    if isfield(data, 'allPerTrialData')
        maxColsAllPerTrialData = max(maxColsAllPerTrialData, size(data.allPerTrialData, 2));
    end
    if isfield(data, 'allSpikesPerBin')
        maxColsAllSpikesPerBin = max(maxColsAllSpikesPerBin, size(data.allSpikesPerBin, 2));
    end
    if isfield(data, 'results')
        maxColsResults = max(maxColsResults, size(data.results, 2));  % Assuming results can have varying columns
    end
    if isfield(data, 'fromExcel1')
        maxColsResults = max(maxColsResults, size(data.results, 2));  % Assuming results can have varying columns
    end

    % Loop through each specified variable
    for varIndex = 1:length(variablesToCombine)
        currentVar = variablesToCombine{varIndex};

        % Check if the variable exists in the current file
        if isfield(data, currentVar)
            currentData = data.(currentVar);
            fprintf('Current dimensions of %s before appending from %s: %s\n', currentVar, files{k}, mat2str(size(currentData)));

            if strcmp(currentVar, 'stat_1_v2')
                if size(currentData, 3) < maxLayersStat1_v2
                    % Pad the third dimension with empty cells if it's a cell array, or NaN if numeric
                    if iscell(currentData)
                        currentData = cat(3, currentData, cell(size(currentData, 1), size(currentData, 2), maxLayersStat1_v2 - size(currentData, 3)));
                    else
                        padding = NaN(size(currentData, 1), size(currentData, 2), maxLayersStat1_v2 - size(currentData, 3));
                        currentData = cat(3, currentData, padding);
                    end
                end
            elseif strcmp(currentVar, 'Firing_rates_f') || strcmp(currentVar, 'allPerTrialData') || strcmp(currentVar, 'allSpikesPerBin') || strcmp(currentVar, 'results') || strcmp(currentVar, 'fromExcel1')
                % Determine the appropriate column count based on the variable
            % Determine the appropriate column count based on the variable
            if strcmp(currentVar, 'Firing_rates_f')
                maxCols = maxColsFiringRates_f;
            elseif strcmp(currentVar, 'allPerTrialData')
                maxCols = maxColsAllPerTrialData;
            elseif strcmp(currentVar, 'allSpikesPerBin')
                maxCols = maxColsAllSpikesPerBin;
            elseif strcmp(currentVar, 'results')
                maxCols = maxColsResults;
            elseif strcmp(currentVar, 'fromExcel1')
                maxCols = maxColsResults;    
            else
                maxCols = 0;  % Default case, should not occur unless new variables are added without updating this logic
            end

                if size(currentData, 2) < maxCols
                    if iscell(currentData)
                        currentData = [currentData, cell(size(currentData, 1), maxCols - size(currentData, 2))];
                    else
                        currentData = [currentData, NaN(size(currentData, 1), maxCols - size(currentData, 2))];
                    end
                end
            end

            if isempty(combinedData.(currentVar))
                combinedData.(currentVar) = currentData;
            else
                combinedData.(currentVar) = [combinedData.(currentVar); currentData];

                fprintf('New dimensions of %s after appending: %s\n', currentVar, mat2str(size(combinedData.(currentVar))));
            end
        else
            fprintf('%s not found in %s.\n', currentVar, files{k});
        end
    end
end

% Save the combined data to a new MAT file using version 7.3 for large files
save('combinedData_firstManuscript.mat', '-struct', 'combinedData', '-v7.3');
