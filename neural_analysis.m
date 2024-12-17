clear;
clc;
close all;

% Load from Windows
Base = 'X:\georgios\UROP\';  % adjust to your needs
List1 = dir(fullfile(Base, '*.*'));
List1 = List1([List1.isdir]);
SubFolder = {List1.name};
SubFolder(ismember(SubFolder, {'.', '..'})) = [];
List2  = dir(fullfile(Base, '**', '*.nex'));  % 
List3  = dir(fullfile(Base, '**', '*.nev'));  % 
List4  = dir(fullfile(Base, '**', '*.mat'));  % 

[numbers1, TEXT1, everything1]=xlsread('X:\georgios\UROP\Victoria\Mnks\LittleDebbie\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers2, TEXT2, everything2]=xlsread('X:\georgios\UROP\Victoria\Mnks\Prez\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers3, TEXT3, everything3]=xlsread('X:\georgios\UROP\Michelle\Mnks\LittleDebbie\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers4, TEXT4, everything4]=xlsread('X:\georgios\UROP\Michelle\Mnks\Prez\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers5, TEXT5, everything5]=xlsread('X:\georgios\UROP\Katie\Mnks\LittleDebbie\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers6, TEXT6, everything6]=xlsread('X:\georgios\UROP\Nick\Mnks\LittleDebbie\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers7, TEXT7, everything7]=xlsread('X:\georgios\UROP\Pip\Mnks\LittleDebbie\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers8, TEXT8, everything8]=xlsread('X:\georgios\UROP\Cooper\Mnks\LittleDebbie\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers9, TEXT9, everything9]=xlsread('X:\georgios\UROP\Cooper\Mnks\Prez\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers10, TEXT10, everything10]=xlsread('X:\georgios\UROP\Cristian\Mnks\LittleDebbie\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers11, TEXT11, everything11]=xlsread('X:\georgios\UROP\Cristian\Mnks\Prez\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');
[numbers12, TEXT12, everything12]=xlsread('X:\georgios\UROP\Katie\Mnks\Prez\Recordings\ApAv\Summary_ApAv.xlsx','Sheet2');

everything1(sum(cellfun(@(x) sum(isnan(x)), everything1), 2) >= 5, :) = [];
everything2(sum(cellfun(@(x) sum(isnan(x)), everything2), 2) >= 5, :) = [];
everything3(sum(cellfun(@(x) sum(isnan(x)), everything3), 2) >= 5, :) = [];
everything4(sum(cellfun(@(x) sum(isnan(x)), everything4), 2) >= 5, :) = [];
everything5(sum(cellfun(@(x) sum(isnan(x)), everything5), 2) >= 5, :) = [];
everything6(sum(cellfun(@(x) sum(isnan(x)), everything6), 2) >= 5, :) = [];
everything7(sum(cellfun(@(x) sum(isnan(x)), everything7), 2) >= 5, :) = [];
everything8(sum(cellfun(@(x) sum(isnan(x)), everything8), 2) >= 5, :) = [];
everything9(sum(cellfun(@(x) sum(isnan(x)), everything9), 2) >= 5, :) = [];
everything10(sum(cellfun(@(x) sum(isnan(x)), everything10), 2) >= 5, :) = [];
everything11(sum(cellfun(@(x) sum(isnan(x)), everything11), 2) >= 5, :) = [];
everything12(sum(cellfun(@(x) sum(isnan(x)), everything12), 2) >= 5, :) = [];

everything=[everything1;everything2;everything3;everything4;everything5;everything6;...
    everything7;everything8;everything9;everything10;everything11;everything12];

% Correct strings in column 17, column 18, and now column 19
for i = 1:size(everything, 1)
    % Trim spaces and correct strings in column 17
    trimmedStr17 = strtrim(everything{i, 17});
    if strcmpi(trimmedStr17, 'okunit')
        everything{i, 17} = 'OKunit';
    elseif strcmpi(trimmedStr17, 'goodunit')
        everything{i, 17} = 'goodunit';
    elseif strcmpi(trimmedStr17, 'poorunit')
        everything{i, 17} = 'poorunit';
    elseif strcmpi(trimmedStr17, 'mua')
        everything{i, 17} = 'MUA';
    end

    % Trim spaces and correct strings in column 18
    if ischar(everything{i, 18}) || isstring(everything{i, 18})
        trimmedStr18 = strtrim(everything{i, 18});
        correctStrings18 = {'cOFC', 'pACC', 'dpACC', 'claustrum', 'putamen', 'cOFC/putamen', 'inferior frontal gyrus', 'inferior frontal gyrus/cOFC', 'precentral gyrus'};
        for j = 1:length(correctStrings18)
            if strcmpi(trimmedStr18, correctStrings18{j})
                everything{i, 18} = correctStrings18{j};
            end
        end
    end

    % Trim spaces and correct strings in column 19
    trimmedStr19 = strtrim(everything{i, 19});
    if strcmpi(trimmedStr19, 'prez')
        everything{i, 19} = 'Prez';
    elseif strcmpi(trimmedStr19, 'littledebbie')
        everything{i, 19} = 'LittleDebbie';
    end
end

% Replace NaN or inf values in columns 4 to 15 with 0
for i = 1:size(everything, 1)
    for j = 4:15
        if isnumeric(everything{i, j}) && isscalar(everything{i, j})
            if isnan(everything{i, j}) || isinf(everything{i, j})
                everything{i, j} = 0;
            end
        end
    end
end

% Identify rows where the first cell is NaN (adjusting for mixed data types)
rowsWithNaNInFirstCell = cellfun(@(x) (isnumeric(x) && isnan(x(1))), everything(:, 1), 'UniformOutput', false);

% Convert rowsWithNaNInFirstCell to a logical array
rowsWithNaNInFirstCell_logical = cell2mat(rowsWithNaNInFirstCell);

% Remove those rows from the matrix
everything(rowsWithNaNInFirstCell_logical, :) = [];

% Iterate over each row in the second column to remove leading white spaces
for i = 1:size(everything, 1)
    if ischar(everything{i, 2}) || isstring(everything{i, 2})
        everything{i, 2} = strtrim(everything{i, 2});
    end
end

% Identify rows where the first cell contains the string 'Session'
rowsWithSessionInFirstCell = cellfun(@(x) ischar(x) && contains(x, 'Session'), everything(:, 1));

% Remove those rows from the matrix
everything(rowsWithSessionInFirstCell, :) = [];

% Identify rows where the second cell is numeric and contains NaN
rowsWithNaNInSecondColumn = cellfun(@(x) isnumeric(x) && isnan(x), everything(:, 2));

% Remove those rows from the matrix
everything(rowsWithNaNInSecondColumn, :) = [];


Files = fullfile(List2.folder, List2.name);

Practice=List2;
for i=length(Practice):-1:1
    if ~contains(List2(i).folder,'\Recordings\ApAv\')
        Practice(i) =[];
    end
end

List2=Practice;

% Convert the first column to datetime objects
dates = datetime(everything(:,1), 'InputFormat', 'MM/dd/yyyy');

% Find indices of rows where the year is 2018
indices2018 = year(dates) == 2018;

% Remove these rows from the cell array
everything(indices2018, :) = [];

% Logical indices for rows with 'cOFC' or 'pACC' in the 18th column of 'everything'
indices = strcmp(everything(:, 18), 'cOFC') | strcmp(everything(:, 18), 'pACC');

% Creating the new cell array 'everything_clean1' with only the filtered rows
everything_clean1 = everything(indices, :);

% Additional filtering for rows with 'OKunit' or 'goodunit' in the 17th column of 'everything_clean1'
additional_indices = strcmp(everything_clean1(:, 17), 'OKunit') | strcmp(everything_clean1(:, 17), 'goodunit');

% Applying the additional filter to 'everything_clean1'
everything_clean1 = everything_clean1(additional_indices, :);

everything=everything_clean1;

% Assuming 'everything' is your 1017x19 cell array
numRows = size(everything_clean1, 1);  % Get the number of rows

% Temporary variable to hold the data of column 3
tempColumn = everything_clean1(:, 3);

% Move column 19 to column 3
everything_clean1(:, 3) = everything_clean1(:, 19);

% Move the original column 3 data to column 19
everything_clean1(:, 19) = tempColumn;

% Now 'everything' has column 19 swapped with column 3



% 'everything' now has rows with dates from 2018 removed
% Initialize an array to hold indices of elements to keep
keepIndices = [];

% Loop through each element in List2
for i = 1:length(List2)
    % Extract the folder name
    folderName = List2(i).folder;

    % Find the position of the date in the folder name
    % Assuming the date format is 'yyyy-MM-dd' and comes right after '\ApAv\'
    dateStartPos = strfind(folderName, '\ApAv\') + 6;
    dateEndPos = dateStartPos + 9;

    % Extract the year from the date string
    yearStr = folderName(dateStartPos:dateStartPos+3);

    % Check if the year is not 2018, then add index to keepIndices
    if ~strcmp(yearStr, '2018')
        keepIndices(end + 1) = i;
    end
end

% Overwrite List2 with only the elements to keep
List2 = List2(keepIndices);

% Assuming List2 is your structure with fields 'name' and 'folder'
numEntries = numel(List2);  % Get the actual number of entries
malakas = cell(numEntries, 5);  % Initialize as a cell array with five columns

% Loop through each entry in List2
for i = 1:numEntries
    % Column 1 in malakas: Extract and format the date from the folder field
    dateMatch = regexp(List2(i).folder, '\\ApAv\\(\d{4}-\d{1,2}-\d{1,2})_', 'tokens');
    if ~isempty(dateMatch)
        % Using datetime to parse and format the date string directly
        dateValue = datetime(dateMatch{1}{1}, 'InputFormat', 'yyyy-MM-dd', 'Format', 'M/d/yyyy');
        malakas{i, 1} = char(dateValue);  % Convert datetime to string format directly using char
    end

    % Column 2 in malakas: name data from List2, remove '.nex' suffix if present
    malakas{i, 2} = regexprep(List2(i).name, '\.nex$', '');

    % Column 3 in malakas: Extract animal's name (moved from the original column 4 to 3)
    animalMatch = regexp(List2(i).folder, '\\Mnks\\(.*?)\\', 'tokens');
    if ~isempty(animalMatch)
        malakas{i, 3} = animalMatch{1}{1};
    end

    % Column 4 in malakas: Extract person's name (moved from the original column 3 to 4)
    personMatch = regexp(List2(i).folder, '\\UROP\\(.*?)\\', 'tokens');
    if ~isempty(personMatch)
        malakas{i, 4} = personMatch{1}{1};
    end

    % Column 5 in malakas: Directly copy the folder path from List2
    malakas{i, 5} = List2(i).folder;
end

% Concatenate the first three columns into single string keys for both matrices
getKey = @(x) strcat(x{1}, '_', x{2}, '_', x{3});
everythingKey = cellfun(getKey, num2cell(everything_clean1(:, 1:3), 2), 'UniformOutput', false);
malakasKey = cellfun(getKey, num2cell(malakas(:, 1:3), 2), 'UniformOutput', false);

% Find indices of common rows based on these keys
[commonKeys, idxMalakas] = ismember(malakasKey, everythingKey);

% Filter 'malakas' to keep only the common rows
asihtir = malakas(commonKeys, :);

% Initialize 'epitelous' with the same number of rows as 'asihtir' and 2 columns
numRows = size(asihtir, 1);
epitelous = cell(numRows, 2);

% Populate 'epitelous'
for i = 1:numRows
    % First column: Append '.nex' to each element from column 2 of 'asihtir'
    epitelous{i, 1} = [asihtir{i, 2}, '.nex'];

    % Second column: Copy directly from column 5 of 'asihtir'
    epitelous{i, 2} = asihtir{i, 5};
end

% Number of rows in epitelous
numRows = size(epitelous, 1);

% Initialize the struct array with the specified fields
List2 = struct('name', cell(numRows, 1), 'folder', cell(numRows, 1));

% Loop through each row of epitelous to populate the structure array
for i = 1:numRows
    List2(i).name = epitelous{i, 1};   % Assign data from the first column of epitelous to the 'name' field
    List2(i).folder = epitelous{i, 2}; % Assign data from the second column of epitelous to the 'folder' field
end

% Calculate the total number of channels once before the loop starts
totalChannels = length(List2);

for x = 1:totalChannels
    
    % Print the current channel and the total number of channels
    fprintf('Processing channel %d out of %d\n', x, totalChannels);
     
    generaltitle=([List2(x).folder,'\',List2(x).name]);
    nexFileData = readNexFile([List2(x).folder,'\',List2(x).name]);
   
        % Pattern for matching .nev files
    nevPattern = fullfile(List2(x).folder, '*.nev');
    pattern = '\d{4}'; % Pattern to match four consecutive digits
    yearMatch = regexp(nevPattern, pattern, 'match');
    yearStr = yearMatch{1}; % Extract the first match which should be the year
    yearNum = str2double(yearStr);

    % Find all .nev files that match the pattern in the current folder
    nevFiles = dir(nevPattern);

    % Check if there is exactly one .nev file
    if length(nevFiles) == 1
        % There's exactly one .nev file, so load it
        nevFilePath = fullfile(nevFiles(1).folder, nevFiles(1).name);
        ncs=read_neuralynx_nev(nevFilePath);
    elseif isempty(nevFiles)
        % No .nev files found
        disp(['No .nev files found in ', List2(x).folder]);
    else
        % More than one .nev file found, handle accordingly
        disp(['Multiple .nev files found in ', List2(x).folder, '. Please check the folder.']);
    end
  
            % Convert TimeStamps from microseconds to milliseconds
        for i = 1:length(ncs)
            ncs(i).TimeStamp = ncs(i).TimeStamp / 1000;
        end
        
        % Initialize variables for different trial types
        all_trials = {};
        choice_trials = {};
        approach_trials = {};
        avoidance_trials = {};
        pavlovian_trials = {};
        pavlovian_reward_trials = {};
        pavlovian_airpuff_trials = {};
        bad_trials = {};
        current_trial = []; % Initialize current_trial as empty

        % Initialize counters for the number of trials starting with 9 and 255
        num_trials_9 = 0;
        num_trials_255 = 0;
        num_trials_18=0;
        
        % Loop over all elements in the ncs structure
        for j = 1:length(ncs)
            if ncs(j).TTLValue == 9  % Check for the start of a trial with a 9
                num_trials_9 = num_trials_9 + 1;  % Increment the trial counter for 9
            elseif ncs(j).TTLValue == 255  % Check for the start of a trial with a 255
                num_trials_255 = num_trials_255 + 1;  % Increment the trial counter for 255
            end
        end
                
        % Process trials
        for j = 1:length(ncs)
            if ncs(j).TTLValue == 9  % Trial start
                current_trial = struct('TimeStamps', [], 'TTLValues', [], 'Type', '');
            end
        
                % Continue collecting data if we are within a trial
                if ~isempty(current_trial)
                    current_trial.TimeStamps(end+1) = ncs(j).TimeStamp;
                    current_trial.TTLValues(end+1) = ncs(j).TTLValue;
            
                    % Use switch-case for categorizing trials for better readability
                    switch ncs(j).TTLValue
                        case 233
                            current_trial.Type = 'Choice';
                        case 239
                            if strcmp(current_trial.Type, 'Choice')
                                current_trial.Type = 'Approach';
                            end
                        case 244
                            if strcmp(current_trial.Type, 'Choice')
                                current_trial.Type = 'Avoidance';
                            end
                        case {211, 215}
                            current_trial.Type = 'Pavlovian';
                        case 217
                            if strcmp(current_trial.Type, 'Pavlovian')
                                current_trial.Type = 'fReward';
                            end
                        case 238
                            if strcmp(current_trial.Type, 'Pavlovian')
                                current_trial.Type = 'fAirpuff';
                            end
                        case {202, 203, 204, 205, 206}
                            current_trial.Type = 'BadTrial';
                    end
                end
            
        
             if ncs(j).TTLValue == 18  % Trial end
                % Only save the trial if current_trial has been started properly
                if ~isempty(current_trial) && any(current_trial.TTLValues == 9)
                    all_trials{end+1} = current_trial; % Save all trials
            
                    % Categorize and save the trial based on its type
                    switch current_trial.Type
                        case 'Choice'
                            choice_trials{end+1} = current_trial;
                        case 'Approach'
                            approach_trials{end+1} = current_trial;
                        case 'Avoidance'
                            avoidance_trials{end+1} = current_trial;
                        case 'Pavlovian'
                            % Logic for categorizing Pavlovian trials
                            % (further refined into 'fReward' or 'fAirpuff' if necessary)
                            pavlovian_trials{end+1} = current_trial;
                        case 'fReward'
                            pavlovian_reward_trials{end+1} = current_trial;
                        case 'fAirpuff'
                            pavlovian_airpuff_trials{end+1} = current_trial;
                        case 'BadTrial'
                            bad_trials{end+1} = current_trial;
                        otherwise
                            % Handle any other unexpected trial types
                    end
                else
                    % Log if the start event is missing
                    fprintf('Warning: End of trial detected without start event for index %d\n', j);
                end
                current_trial = []; % Reset for the next trial
            end
        end
         % Initialize reward and punishment percentages for approach trials
        for i = 1:length(approach_trials)
            approach_trials{i}.RewardPercentage = NaN;
            approach_trials{i}.PunishmentPercentage = NaN;
        end
        
        % Initialize reward and punishment percentages for avoidance trials
        for i = 1:length(avoidance_trials)
            avoidance_trials{i}.RewardPercentage = NaN;
            avoidance_trials{i}.PunishmentPercentage = NaN;
        end
        
        % Initialize reward and punishment percentages for pavlovian reward trials
        for i = 1:length(pavlovian_reward_trials)
            pavlovian_reward_trials{i}.RewardPercentage = NaN;
            pavlovian_reward_trials{i}.PunishmentPercentage = NaN;
        end
        
        % Initialize reward and punishment percentages for pavlovian airpuff trials
        for i = 1:length(pavlovian_airpuff_trials)
            pavlovian_airpuff_trials{i}.RewardPercentage = NaN;
            pavlovian_airpuff_trials{i}.PunishmentPercentage = NaN;
        end
        
        % Initialize reward and punishment percentages for bad trials
        for i = 1:length(bad_trials)
            bad_trials{i}.RewardPercentage = NaN;
            bad_trials{i}.PunishmentPercentage = NaN;
        end
        
        for trial_type = {'approach_trials', 'avoidance_trials', 'pavlovian_reward_trials', 'pavlovian_airpuff_trials', 'bad_trials'}
            trials = eval(trial_type{1});
            for i = 1:length(trials)
                % Check and flag if it's a bad trial
                trials{i}.IsBadTrial = any(ismember(trials{i}.TTLValues, 202:206));
        
                % Find the index for the target event marker (251, 252, or 254)
                target_index = find(ismember(trials{i}.TTLValues, [251, 252, 254]), 1, 'last');
        
                if ~isempty(target_index) && target_index < length(trials{i}.TTLValues) - 1
                    % Extract reward and punishment sizes
                    trials{i}.RewardSize = trials{i}.TTLValues(target_index + 1);
                    trials{i}.PunishmentSize = trials{i}.TTLValues(target_index + 2);
        
                    % Calculate percentages
                    trials{i}.RewardPercentage = (trials{i}.RewardSize / 200) * 100;
                    trials{i}.PunishmentPercentage = (trials{i}.PunishmentSize / 200) * 100;
                else
                    % Set default values if target marker is missing or too close to end
                    trials{i}.RewardSize = 0;
                    trials{i}.PunishmentSize = 0;
                    trials{i}.RewardPercentage = 0;
                    trials{i}.PunishmentPercentage = 0;
                end
            end
            % Save back the processed trials
            eval([trial_type{1} ' = trials;']);
        end


        % ... (rest of the code, including the display of counts and percentages)
        % Initialize the results matrix with the total number of trials
        total_choice_trials = length(approach_trials) + length(avoidance_trials) + length(pavlovian_airpuff_trials)+...
            length(pavlovian_reward_trials) + length(bad_trials);
        results_matrix = zeros(total_choice_trials, 3);
        
        % Initialize counters for each trial type
        approach_index = 1;
        avoidance_index = 1;
        pavlovian_reward_index = 1;
        pavlovian_airpuff_index = 1;
        bad_trial_index = 1;
        result_index = 1; % Keep track of the position in the results_matrix
        
        % Total number of trials
        total_num_trials = length(approach_trials) + length(avoidance_trials) + ...
                           length(pavlovian_reward_trials) + length(pavlovian_airpuff_trials) + ...
                           length(bad_trials);
        results_matrix = zeros(total_num_trials, 3); % Initialize the results_matrix
        
        % Loop until we have placed all trials in the results_matrix
        while approach_index <= length(approach_trials) || avoidance_index <= length(avoidance_trials) || ...
              pavlovian_reward_index <= length(pavlovian_reward_trials) || pavlovian_airpuff_index <= length(pavlovian_airpuff_trials) || ...
              bad_trial_index <= length(bad_trials)
            
            % Determine the earliest trial among the remaining trials
            [~, min_index] = min([getFirstTimeStamp(approach_trials, approach_index), ...
                                  getFirstTimeStamp(avoidance_trials, avoidance_index), ...
                                  getFirstTimeStamp(pavlovian_reward_trials, pavlovian_reward_index), ...
                                  getFirstTimeStamp(pavlovian_airpuff_trials, pavlovian_airpuff_index), ...
                                  getFirstTimeStamp(bad_trials, bad_trial_index)]);
            
            % Add the earliest trial to the results matrix and increment the corresponding index
            switch min_index
                case 1 % Approach trial
                    results_matrix(result_index, :) = [approach_trials{approach_index}.RewardPercentage, ...
                                                       approach_trials{approach_index}.PunishmentPercentage, ...
                                                       1]; % 1 for approach
                    approach_index = approach_index + 1;
                case 2 % Avoidance trial
                    results_matrix(result_index, :) = [avoidance_trials{avoidance_index}.RewardPercentage, ...
                                                       avoidance_trials{avoidance_index}.PunishmentPercentage, ...
                                                       0]; % 0 for avoidance
                    avoidance_index = avoidance_index + 1;
                case 3 % Pavlovian reward trial
                    results_matrix(result_index, :) = [pavlovian_reward_trials{pavlovian_reward_index}.RewardPercentage, ...
                                                       pavlovian_reward_trials{pavlovian_reward_index}.PunishmentPercentage, ...
                                                       2]; % 2 for Pavlovian reward
                    pavlovian_reward_index = pavlovian_reward_index + 1;
                case 4 % Pavlovian airpuff trial
                    results_matrix(result_index, :) = [pavlovian_airpuff_trials{pavlovian_airpuff_index}.RewardPercentage, ...
                                                       pavlovian_airpuff_trials{pavlovian_airpuff_index}.PunishmentPercentage, ...
                                                       3]; % 3 for Pavlovian airpuff
                    pavlovian_airpuff_index = pavlovian_airpuff_index + 1;
                case 5 % Bad trial
                    % Check if RewardSize and PunishmentSize fields exist
                    if isfield(bad_trials{bad_trial_index}, 'RewardSize') && isfield(bad_trials{bad_trial_index}, 'PunishmentSize')
                        % Calculate the percentages for bad trials
                        rewardPercentage = (bad_trials{bad_trial_index}.RewardSize / 200) * 100;
                        punishmentPercentage = (bad_trials{bad_trial_index}.PunishmentSize / 200) * 100;
                    else
                        % Handle the case where fields are missing - set percentages to zero or another default value
                        rewardPercentage = 0;
                        punishmentPercentage = 0;
                        fprintf('Warning: RewardSize or PunishmentSize field missing in bad trial at index %d\n', bad_trial_index);
                    end
                    
                    % Add the calculated percentages to the results matrix
                    results_matrix(result_index, :) = [rewardPercentage, punishmentPercentage, 4]; % 4 for bad trial
                    bad_trial_index = bad_trial_index + 1;
                    end
        
            result_index = result_index + 1; % Move to the next row in the results matrix
        end
        
        % Filter the results_matrix to keep only rows where the 3rd column is 0 or 1
        choice_rows = results_matrix(:, 3) == 0 | results_matrix(:, 3) == 1;
        
        % Create results_choice_matrix with only the filtered rows
        results_choice_matrix = results_matrix(choice_rows, :);

      % 'avoidance_trials' with Avoidance trials within choice trials
        [b, ~, ~] = glmfit(results_choice_matrix(:,1:2), results_choice_matrix(:,3), 'binomial', 'link', 'logit');
        
        betaValues=b;
        plot_x = 0:100;
        plot_y = (-b(1) - b(2).*plot_x)./b(3);
        plot_y1=plot_y;
        
        % Extract all approach trials (where the third column is 1)
        approach_matrix = results_choice_matrix(results_choice_matrix(:,3) == 1, 1:2);
        approach_matrix1=approach_matrix;
        
        % Extract all avoidance trials (where the third column is 0)
        avoidance_matrix = results_choice_matrix(results_choice_matrix(:,3) == 0, 1:2);
        avoidance_matrix1=avoidance_matrix;
        
        % Define the event markers for start and end of reaction time for each type of trial
        start_markers =  [234, 247];
        end_marker_approach = 219;
        end_marker_avoidance = 243;
        
        % Initialize matrices to store reaction times
        approach_matrix_rt = [];
        avoidance_matrix_rt = [];
        
        % Loop through all trials to calculate reaction times
        for i = 1:length(all_trials)
            % Extract the event markers and timestamps for the current trial
            event_markers = all_trials{i}.TTLValues;
            timestamps = all_trials{i}.TimeStamps;
        
            % Find the indices of the start and end markers for approach and avoidance
            start_index = find(ismember(event_markers, start_markers), 1, 'first');
            end_index_approach = find(event_markers == end_marker_approach, 1, 'first');
            end_index_avoidance = find(event_markers == end_marker_avoidance, 1, 'first');
        
            % Calculate reaction time for approach if the markers are present
            if ~isempty(start_index) && ~isempty(end_index_approach)
                rt_approach = timestamps(end_index_approach) - timestamps(start_index);
                approach_matrix_rt = [approach_matrix_rt; rt_approach];
                approach_matrix_rt1=(approach_matrix_rt);
            end
            
            % Calculate reaction time for avoidance if the markers are present
            if ~isempty(start_index) && ~isempty(end_index_avoidance)
                rt_avoidance = timestamps(end_index_avoidance) - timestamps(start_index);
                avoidance_matrix_rt = [avoidance_matrix_rt; rt_avoidance];
                avoidance_matrix_rt1=(avoidance_matrix_rt);
            end
        end

    %     sprintf('%s%d',[List2(x).folder,'\',List2(x).name]);
            
        for ph=1:size(nexFileData.neurons,1)
        
            Events_timestamps=nexFileData.markers{1,1}.timestamps;
            Events_TTL_value=str2double(nexFileData.markers{1,1}.values{1,1}.strings);
            
            FindErrors=find(Events_TTL_value==214); %initial fixation onset for all trials
            FindErrors_plus1=FindErrors+1;
    
            if size(FindErrors_plus1,1)>0
                erotas=FindErrors_plus1;
                if FindErrors_plus1(end)>numel(Events_TTL_value)
                    erotas(end)=[];
                end
                FindErrors_plus1=erotas;
            end
    
            BadSpots=Events_TTL_value(FindErrors_plus1);
            BadSpots=[BadSpots,FindErrors_plus1];
            error=find(BadSpots(:,1)==18);
            error_idx=BadSpots(error,2);
            Events_TTL_value(error_idx)=204;
            Events_all=[Events_timestamps,Events_TTL_value];
     
            ev=[];
        
            rangeOfEvents=201:255;  
    
            for z=rangeOfEvents
                [ev_index,~] = find(Events_all(:,2) == z);
                ev{1,z}=Events_timestamps(ev_index);
            end
    
            startIndex = 1;
            endIndex = 200;
            ev(startIndex:endIndex) = [];
            
            yep = num2cell(201:255);
           
            villos=[ev;yep];
            event_names=villos(2,:);    
        
            Spikes=(nexFileData.neurons{ph,1}.timestamps)';
            Timestamps=[nexFileData.neurons{ph,1}.timestamps,ev];
            Spikes_per_bin=[];

            for s=1:size(ev,2)
            
                Timestamps_minus=(ev{1,s}-2)';
                Timestamps_plus=(ev{1,s}+2)';
                Timestamps_zero=(ev{1,s}+0)';
                
                z=[];
                % spike_index=[];
            
                for z=1:size(Timestamps_plus,2)
                    spike_index{z} = find(Spikes>Timestamps_minus(z) & Spikes<Timestamps_plus(z));
                    logicalIndexes{z} = (Spikes>Timestamps_minus(z) & Spikes<Timestamps_plus(z));
                    linearIndexes{z} = find(logicalIndexes{z});
                end
            
                k=[];
            
                for k=1:size(spike_index,2)
                    vector_size(k)= size(spike_index{1,k},2);
                end
            
                Spikes_num=sum(vector_size);
                
            full_range=2:-0.05:-2;
            
            j=[];
            full_bin=[];
        
            for j=1:size(Timestamps_minus,2)
                full_bin(j,:)=Timestamps_zero(1,j)-full_range;
            end
        
            u=[];
            p=[];
            
            spike_index_full=[];
            logicalIndexes_full=[];
            linearIndexes_full=[];
        
            for u=1:size(full_bin,1)
                for p=1:size(full_bin,2)-1
                    spike_index_full{u,p}=find(Spikes>full_bin(u,p) & Spikes<full_bin(u,p+1));
                    logicalIndexes_full = (Spikes>full_bin(u,p) & Spikes<full_bin(u,p+1));
                    linearIndexes_full{u,p} = find(logicalIndexes_full);
                end
            end

%%% find the spike timings in the specific time window
            spike_time_full=[];
            logicalTimeIndexes_full=[];
            linearTimeIndexes_full=[];
            SpikeTimes=Spikes;
        
            for u=1:size(full_bin,1)
                for p=1:size(full_bin,2)-1
                    isnan((Spikes>full_bin(u,p) & Spikes<full_bin(u,p+1)));
                end
            end
%%%%%%%%%%%%      
            k=[];
            n=[];
            vector_size_full=[];
        
            for k=1:size(spike_index_full,1)
                for n=1:size(spike_index_full,2)
                    vector_size_full(k,n)= size(spike_index_full{k,n},2);
                end
            end

            if size(vector_size_full,1) < 2
                Partial_sum=vector_size_full;
            else
                Partial_sum=sum(vector_size_full);
            end
        
            if size(Partial_sum,1)>0
                Spikes_per_bin(s,:)=Partial_sum;
                AllperTrialData{s}=vector_size_full;
                Spikes_full_sum{s}=sum(sum(vector_size_full));
                Partial_sum=[];
            else
                Partial_sum=nan(1, 80);
                Spikes_per_bin(s,:)=Partial_sum;
                AllperTrialData{s}=vector_size_full;
                Spikes_full_sum{s}=sum(sum(vector_size_full));
                Partial_sum=[];
            end
        
        end
              allPerTrialData{x,ph} = AllperTrialData;
              allSpikesPerBin{x,ph} = Spikes_per_bin;  % Stores the current 'Spikes_per_bin' in the cell array
            
             for l=1:size(Spikes_per_bin,1)
                Stats_1(x,l,ph)=sum(Spikes_per_bin(l,[41:70]))/sum(Spikes_per_bin(l,[11:40]));
                stat_1_v2(x,l,ph)=cellstr(num2str(dg_chi2test2([sum(Spikes_per_bin(l,[41:70])),sum(Spikes_per_bin(l,[11:40]))])));
                folder(x,ph)=cellstr(List2(x).folder);
                name(x,ph)=cellstr(List2(x).name);
             end

             % Check if the first three elements match 0, 214, and 233
            if length(Events_TTL_value) >= 3 && Events_TTL_value(1) == 0 && Events_TTL_value(2) == 214 && Events_TTL_value(3) == 233
                % Add two zeros at the beginning of the vector
                Events_TTL_value = [0; 0; Events_TTL_value];
            end


                % Define the event markers for start and end of reaction time for each type of trial
                start_marker = 234;
                end_marker_approach = 219;
                end_marker_avoidance = 243;
                
                % Initialize a matrix to store reaction times and trial type
                all_trials_rt = []; % Each row will have [RT, Trial Type], where Trial Type: 1 for approach, 0 for avoidance
                
                % Loop through all trials to calculate reaction times
                for i = 1:length(all_trials)
                    % Extract the event markers and timestamps for the current trial
                    event_markers = all_trials{i}.TTLValues;
                    timestamps = all_trials{i}.TimeStamps;
                
                    % Find the indices of the start and end markers for approach and avoidance
                    start_index = find(event_markers == start_marker, 1, 'first');
                    end_index_approach = find(event_markers == end_marker_approach, 1, 'first');
                    end_index_avoidance = find(event_markers == end_marker_avoidance, 1, 'first');
                
                    % Check and calculate reaction time for approach or avoidance
                    if ~isempty(start_index)
                        if ~isempty(end_index_approach)
                            rt_approach = timestamps(end_index_approach) - timestamps(start_index);
                            all_trials_rt = [all_trials_rt; rt_approach, 1]; % 1 for approach
                        elseif ~isempty(end_index_avoidance)
                            rt_avoidance = timestamps(end_index_avoidance) - timestamps(start_index);
                            all_trials_rt = [all_trials_rt; rt_avoidance, 0]; % 0 for avoidance
                        end
                    end
                end
                
                % all_trials_rt now contains the reaction time and trial type for each trial in order
                
                
                % Extract the data for the first set of events
                dataYellowBar = AllperTrialData{1,11}; % Pavlovian Yellow bar onset
                dataRedBar = AllperTrialData{1,15};    % Pavlovian Red bar onset
                dataChoiceCue = AllperTrialData{1,33}; % Choice cue onset
                
                % Define the bin range for the first calculation
                binStart1 = 41;
                binEnd1 = 70;
                numBins1 = binEnd1 - binStart1 + 1;
   
                if yearNum >2018
                    % Calculate average firing rates for each trial of each event (first set)
                    firingRatesYellowBar = sum(dataYellowBar(:, binStart1:binEnd1), 2) / numBins1 * (1000 / 50);
                    firingRatesRedBar = sum(dataRedBar(:, binStart1:binEnd1), 2) / numBins1 * (1000 / 50);
                    firingRatesChoiceCue = sum(dataChoiceCue(:, binStart1:binEnd1), 2) / numBins1 * (1000 / 50);
                    
                    % Assuming 'all_trials' is your array or cell array of trials
                    % and each trial has a field 'TTLValues' which is an array of event markers
                    
                    % Find indices of trials with both event markers 233 and 234
                    indices_233_and_234 = [];
                    for i = 1:length(all_trials)
                        if (any(all_trials{i}.TTLValues == 233) && any(all_trials{i}.TTLValues == 234) && any(all_trials{i}.TTLValues == 239)) ||...
                                (any(all_trials{i}.TTLValues == 233) && any(all_trials{i}.TTLValues == 234) && any(all_trials{i}.TTLValues == 244))
                            indices_233_and_234(end + 1) = i;
                        end
                    end
%                     fprintf('Number of trials with both event markers 233 and 234: %d\n', length(indices_233_and_234));
                    
                    % Find indices of trials with event markers 233, 234 and either 203 or 206
                    indices_233_234_and_203_or_206 = [];
                    for i = 1:length(all_trials)
                        if any(all_trials{i}.TTLValues == 233) && any(all_trials{i}.TTLValues == 234) && ...
                           (any(all_trials{i}.TTLValues == 203) || any(all_trials{i}.TTLValues == 206))
                            indices_233_234_and_203_or_206(end + 1) = i;
                        end
                    end
    % %                 fprintf('Number of trials with event markers 233 and 234 and either 203 or 206: %d\n', length(indices_233_234_and_203_or_206));
                    
                    % Find the relative indices of trials in indices_233_234_and_203_or_206 within indices_233_and_234
                    [~, relative_indices] = ismember(indices_233_234_and_203_or_206, indices_233_and_234);
                    relative_indices = relative_indices(relative_indices > 0); % Filter out zero values
                    
                    % Create a logical index array where all elements are initially set to true
                    keep_indices = true(length(firingRatesChoiceCue), 1);
                    
                    % Set the elements at relative_indices to false (these are the rows to be removed)
                    keep_indices(relative_indices) = false;
                    
                    % Create a new vector excluding the specified rows
                    firingRatesChoiceCue_filtered = firingRatesChoiceCue(keep_indices);
                    
                    % Extract the data for the second set of events
                    dataPavlovianReward = AllperTrialData{1,17}; % Pavlovian reward delivery
                    dataChoiceCue = AllperTrialData{1,33};   % Choice cue all

                    dataPavlovianAirpuff = AllperTrialData{1,38}; % Pavlovian airpuff delivery
                    dataChoiceAirpuff = AllperTrialData{1,39}; % Choice small reward delivery
                    dataChoiceSmallReward = AllperTrialData{1,44}; % Choice small reward delivery
                    
                    if length(dataChoiceSmallReward)+length(dataChoiceAirpuff) == length(results_choice_matrix(:,3)) && ...
                            sum(results_choice_matrix(:,3)==0)==length(AllperTrialData{1,44})&&...
                            sum(results_choice_matrix(:,3)==1)==length(AllperTrialData{1,39})

                        dataChoiceOutcome = zeros(length(results_choice_matrix(:,3)), 80);

                                            % Counters for approach and avoidance trials
                        approachCounter = 1;
                        avoidanceCounter = 1;
                        
                        for i = 1:length(results_choice_matrix(:,3))
                            if results_choice_matrix(i,3) == 1
                                % If it's an approach trial
                                dataChoiceOutcome(i,:) = dataChoiceAirpuff(approachCounter,:);
                                approachCounter = approachCounter + 1;
                            else
                                % If it's an avoidance trial
                                dataChoiceOutcome(i,:) = dataChoiceSmallReward(avoidanceCounter,:);
                                avoidanceCounter = avoidanceCounter + 1;
                            end
                        end
                    else

                                % Assuming Events_TTL_value and Events_timestamps are already defined
                                n = length(Events_TTL_value);
                                approach_sequence = [233, 234, 219, 201, 239, 240, 241, 242];
                                avoidance_sequence = [233, 234, 243, 207, 244, 245];
                                reaction_times = []; % To store reaction times
                                x1other = []; % To store x1other values
                                x2other = []; % To store x2other values
                                sequence_of_choices = []; % To store sequence choices (1 for approach, 0 for avoidance)
                                
                                i = 1;
                                while i <= n - max(length(approach_sequence), length(avoidance_sequence)) + 1
                                    % Check for approach sequence
                                    if isequal(Events_TTL_value(i:i+length(approach_sequence)-1), approach_sequence')
                                        % Calculate reaction time for approach
                                        idx_234_approach = i + 1;
                                        idx_219 = i + 2;
                                        reaction_time = Events_timestamps(idx_219) - Events_timestamps(idx_234_approach);
                                        reaction_times(end+1) = reaction_time;
                                
                                        % Calculate x1other and x2other
                                        if i > 2
                                            x1other(end+1) = Events_TTL_value(i-3);
                                            x2other(end+1) = Events_TTL_value(i-2);
                                        end
                                
                                        % Update sequence of choices
                                        sequence_of_choices(end+1) = 1; % 1 for approach
                                
                                        i = i + length(approach_sequence);
                                        continue;
                                    end
                                    
                                    % Check for avoidance sequence
                                    if isequal(Events_TTL_value(i:i+length(avoidance_sequence)-1), avoidance_sequence')
                                        % Calculate reaction time for avoidance
                                        idx_234_avoidance = i + 1;
                                        idx_243 = i + 2;
                                        reaction_time = Events_timestamps(idx_243) - Events_timestamps(idx_234_avoidance);
                                        reaction_times(end+1) = reaction_time;
                                
                                        % Calculate x1other and x2other
                                        if i > 2
                                            x1other(end+1) = Events_TTL_value(i-3);
                                            x2other(end+1) = Events_TTL_value(i-2);
                                        end
                                
                                        % Update sequence of choices
                                        sequence_of_choices(end+1) = 0; % 0 for avoidance
                                
                                        i = i + length(avoidance_sequence);
                                        continue;
                                    end
                                
                                    i = i + 1;
                                end
                         sequence_of_choices=sequence_of_choices';
                         dataChoiceOutcome = zeros(length(sequence_of_choices), 80);

                                            % Counters for approach and avoidance trials
                        approachCounter = 1;
                        avoidanceCounter = 1;
                        
                        for i = 1:length(sequence_of_choices)
                            if sequence_of_choices(i,1) == 1
                                % If it's an approach trial
                                dataChoiceOutcome(i,:) = dataChoiceAirpuff(approachCounter,:);
                                approachCounter = approachCounter + 1;
                            else
                                % If it's an avoidance trial
                                dataChoiceOutcome(i,:) = dataChoiceSmallReward(avoidanceCounter,:);
                                avoidanceCounter = avoidanceCounter + 1;
                            end
                        end

                    end

                    % Define the sequences
                            approach_sequence = [233, 234, 219, 201, 239, 240, 241, 242];
                            avoidance_sequence = [233, 234, 243, 207, 244, 245];
                            
                            % Find all occurrences of event marker 233
                            occurrences_233 = find(Events_TTL_value == 233);
                            
                            % Initialize a logical vector to mark bad trials
                            badTrials = false(length(occurrences_233), 1);
                            
                            % Loop through each occurrence of 233
                            for i = 1:length(occurrences_233)
                                idx = occurrences_233(i); % Current index of 233 event marker
                                
                                % Check if the approach_sequence or avoidance_sequence follows this occurrence
                                % Making sure not to exceed the bounds of Events_TTL_value
                                endIdxApproach = min(idx + length(approach_sequence) - 1, length(Events_TTL_value));
                                endIdxAvoidance = min(idx + length(avoidance_sequence) - 1, length(Events_TTL_value));
                                
                                sequenceAfter233Approach = Events_TTL_value(idx:endIdxApproach);
                                sequenceAfter233Avoidance = Events_TTL_value(idx:endIdxAvoidance);
                                
                                % Check for matching sequences
                                isApproach = isequal(sequenceAfter233Approach', approach_sequence(1:length(sequenceAfter233Approach)));
                                isAvoidance = isequal(sequenceAfter233Avoidance', avoidance_sequence(1:length(sequenceAfter233Avoidance)));
                                
                                % Mark as bad trial if neither sequence matches
                                if ~isApproach && ~isAvoidance
                                    badTrials(i) = true;
                                end
                            end
                            
                            if length(dataChoiceCue)==length(badTrials) 
                            % Filter out the bad trials from dataChoiceCue
                                dataChoiceCue = dataChoiceCue(~badTrials, :);
                            end

                    % Define the bin range for the second calculation
                    binStart2 = 41;
                    binEnd2 = 70;
                    numBins2 = binEnd2 - binStart2 + 1;
                    
                    % Calculate average firing rates for each trial of each event (second set)
                    firingRatesPavlovianReward = sum(dataPavlovianReward(:, binStart2:binEnd2), 2) / numBins2 * (1000 / 50);
                    firingRatesChoiceCue = sum(dataChoiceCue(:, binStart2:binEnd2), 2) / numBins2 * (1000 / 50);
                    firingRatesPavlovianAirpuff = sum(dataPavlovianAirpuff(:, binStart2:binEnd2), 2) / numBins2 * (1000 / 50);
                    firingRatesChoiceAirpuff = sum(dataChoiceAirpuff(:, binStart2:binEnd2), 2) / numBins2 * (1000 / 50);
                    firingRatesChoiceSmallReward = sum(dataChoiceSmallReward(:, binStart2:binEnd2), 2) / numBins2 * (1000 / 50);
                    firingRatesChoiceOutcome=sum(dataChoiceOutcome(:, binStart2:binEnd2), 2) / numBins2 * (1000 / 50);
                    
                    %%%%%%%% modeling                    
                    % For pPlus calculation - Option 1: Direct Observation
                    numberOfPositiveChoices = numel(approach_trials); % Enter the number of times the positive option was chosen
                    totalChoices = numel(approach_trials)+numel(avoidance_trials); % Enter the total number of choices made
                    
                    % Coefficients from your model
                    a1 = b(2); % Coefficient for the x-component in UPlus
                    a2 = b(3); % Coefficient for the y-component in UPlus
                    a3 = b(1); % Coefficient for UAV (constant negative value for USquare)
            
                if length(dataChoiceSmallReward)+length(dataChoiceAirpuff) == length(results_choice_matrix(:,3)) && ...
                        sum(results_choice_matrix(:,3)==0)==length(AllperTrialData{1,44})&&...
                        sum(results_choice_matrix(:,3)==1)==length(AllperTrialData{1,39})

                    % Variables for UPlus calculation
                    % Replace these with your actual data for each trial
                    x1 = results_choice_matrix(:,1); % Measure for x-component (e.g., reward size)
                    y = results_choice_matrix(:,2); % Measure for y-component (e.g., airpuff size) 
                else
                    x1=x1other;
                    y=x2other;
                end

                    % Calculating UPlus (Utility of Positive Option)
                    UPlus = a1 * x1 + a2 * y;
                    
                    % Utility of Square Option (UAV) is constant
                    USquare = -a3;
                                        
                    % Calculate probability of choosing the positive option (pPlus) based on utilities
                    pPlus = 1 ./ (1 + exp(-(UPlus - USquare))); % Utility-based probability
                    pPlus = max(min(pPlus, 1 - eps), eps); % Ensure pPlus is within (0, 1)

                    
                    % Calculate Eutil for Ap-Av Task based on pPlus and utilities
                    Eutil_ApAv = pPlus .* UPlus + (1 - pPlus) .* USquare;

                    
                    % Calculate Conflict (Conf) based on decision-making entropy for each trial
                    Conf = -pPlus .* log(pPlus) - (1 - pPlus) .* log(1 - pPlus);
                    
                    % --- Output Section ---
                     % Extract the choice data
                    choices = results_choice_matrix(:, 3); % Replace 'choice_column' with the actual column index

                    % Assuming Events_all is your 13100x2 matrix where the first column contains timings
                    % and the second column contains event markers
                    
                    % Initialize an empty vector to store reaction times
                    RT_hello = [];
                    
                    % Loop through the Events_all matrix to find matching event sequences
                    for i = 1:size(Events_all, 1)-1
                        % Check if the current event is 234 and the next event is either 243 or 219
                        if Events_all(i, 2) == 234 && (Events_all(i+1, 2) == 243 || Events_all(i+1, 2) == 219)
                            % Calculate the time difference (reaction time) and add it to the RT_hello vector
                            RT_hello = [RT_hello; Events_all(i+1, 1) - Events_all(i, 1)];
                        end
                    end
                    
                    % Now, Conf contains the conflict value for each trial
                    % Open or create a log file to append error messages
                    logFileID = fopen('errorLog_step_georgiosRoom.txt', 'a');
                    
                    % First check: size(all_trials_rt,1) vs. size(results_choice_matrix,1)

                    if length(dataChoiceSmallReward)+length(dataChoiceAirpuff) == length(results_choice_matrix(:,3)) && ...
                            sum(results_choice_matrix(:,3)==0)==length(AllperTrialData{1,44})&&...
                            sum(results_choice_matrix(:,3)==1)==length(AllperTrialData{1,39}) &&...
                            length(Conf)==length(all_trials_rt(:, 1))

                            % Second check: size(X1,1) vs. size(Y,1) (and similarly for Y1)
                            % Prepare the data for the second check
                            dataForSW = [results_choice_matrix(:,1:2), Eutil_ApAv, results_choice_matrix(:,3), ...
                                         results_choice_matrix(:,1).*results_choice_matrix(:,3), ...
                                         results_choice_matrix(:,2).*results_choice_matrix(:,3), Conf, all_trials_rt(:, 1)];
                    
                    elseif length(dataChoiceSmallReward)+length(dataChoiceAirpuff) == length(results_choice_matrix(:,3)) && ...
                            sum(results_choice_matrix(:,3)==0)==length(AllperTrialData{1,44})&&...
                            sum(results_choice_matrix(:,3)==1)==length(AllperTrialData{1,39}) &&...
                            length(Conf)==length(RT_hello)

                            dataForSW = [results_choice_matrix(:,1:2), Eutil_ApAv, results_choice_matrix(:,3), ...
                                         results_choice_matrix(:,1).*results_choice_matrix(:,3), ...
                                         results_choice_matrix(:,2).*results_choice_matrix(:,3), Conf, RT_hello];
                    else

                         dataForSW = [x1other',x2other', Eutil_ApAv', sequence_of_choices, ...
                                         x1other'.*sequence_of_choices, ...
                                         x2other'.*sequence_of_choices, Conf', reaction_times'];
                     end
                        X1 = dataForSW(:, 1:8); % Predictor variables
                        Y = firingRatesChoiceCue(:, 1); % Firing rates for Choice Cue
                        Y1 = firingRatesChoiceOutcome(:, 1); % Firing rates for Choice Outcome
                        
                    % Check if sizes match for the first check
                    if size(X1, 1) == size(Y, 1)
                        % Perform stepwise regression with stepwisefit for the first case
                        [b,SE,PVAL,inmodel,stats] = stepwisefit(X1, Y);
                        % Store the results in a structured way
                        results{x,ph}.Coefficients = b;
                        results{x,ph}.SE = SE;
                        results{x,ph}.PValues = PVAL;
                        results{x,ph}.InModel = inmodel;
                        results{x,ph}.Stats = stats;

                        % Perform the correlation analysis
                        [correlationCoefficients, pValues] = corr(X1, Y, 'Rows', 'complete');
                        
                        % Store the correlation results in a structured way similar to the stepwisefit results
                        results{x,ph}.CorrelationCoefficients = correlationCoefficients;
                        results{x,ph}.CorrelationPValues = pValues;
                        
                        mdl{x,ph} = stepwiselm(X1, Y, 'linear', 'Criterion', 'BIC');
                     else
                        % Log an error message if sizes do not match and skip this session
                        errMsg = sprintf('Error: Mismatch in sizes for X1 vs. Y for x=%d, ph=%d.\n', x, ph);
                        fprintf(logFileID, errMsg);
                        % Assuming there's a loop or further code, include necessary control to skip to the next iteration
                    end                    % Close the log file

                    if size(X1, 1) == size(Y1, 1)
                        % Perform stepwise regression with stepwisefit for the second case
                        [b,SE,PVAL,inmodel,stats] = stepwisefit(X1, Y1);
                        % Store the results in a structured way for the second case
                        results1{x,ph}.Coefficients = b;
                        results1{x,ph}.SE = SE;
                        results1{x,ph}.PValues = PVAL;
                        results1{x,ph}.InModel = inmodel;
                        results1{x,ph}.Stats = stats;

                        % Perform the correlation analysis
                        [correlationCoefficients, pValues] = corr(X1, Y1, 'Rows', 'complete');
                        
                        % Store the correlation results in a structured way similar to the stepwisefit results
                        results1{x,ph}.CorrelationCoefficients = correlationCoefficients;
                        results1{x,ph}.CorrelationPValues = pValues; 
                        mdl1{x,ph} = stepwiselm(X1, Y1, 'linear', 'Criterion', 'BIC');
                    else
                        % Log an error message if sizes do not match and skip this session
                        errMsg = sprintf('Error: Mismatch in sizes for X1 vs. Y1 for x=%d, ph=%d.\n', x, ph);
                        fprintf(logFileID, errMsg);
                        % Assuming there's a loop or further code, include necessary control to skip to the next iteration
                    end                    % Close the log file
                    fclose(logFileID);
                end

                % get the dates from the UROP sessions folder names on annex5
                % Pip's folder
                for iii=1:length(List2)
                    [B(iii)]=convertCharsToStrings(List2(iii).folder);
                    [C(iii)]=convertCharsToStrings(List2(iii).name);
                end
                
                Bnew=B';
                c=Bnew;
    
                Cnew=C';
                d=Cnew;
                
            formatOut = 'mm/dd/yyyy';
            NHP1='LittleDebbie';
            NHP2='Prez';
            matchstr='\LittleDebbie\Recordings\ApAv\';
            matchstr1='\Prez\Recordings\ApAv\';
            matchstr2='.mat';
            matchstr3='Katie';
            matchstr4='Pip';
            matchstr5='Cristian';
            matchstr6='Georgios';
            matchstr7='Elizabeth';
            matchstr8='Cooper';
            matchstr9='Michelle';
            matchstr10='Victoria';
            matchstr11='Nick';
    
            for kkk=1:length(Bnew)
                if  contains(Bnew(kkk),'Katie') && contains(Bnew(kkk),matchstr)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,58:67),'InputFormat','yyyy-MM-dd');
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr3;
                    Monkey{kkk}=NHP1;
                elseif contains(Bnew(kkk),'Pip') && contains(Bnew(kkk),matchstr)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,56:65),'InputFormat','yyyy-MM-dd');
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr4;
                    Monkey{kkk}=NHP1;
                 elseif contains(Bnew(kkk),'Cristian') && contains(Bnew(kkk),matchstr)
                     c=char(Bnew(kkk));
                     dt(kkk)=datetime(c(:,61:70),'InputFormat','yyyy-MM-dd');  
                     newStr(kkk) = erase(Cnew(kkk),matchstr2);
                     student{kkk}=matchstr5;
                     Monkey{kkk}=NHP1;
                 elseif contains(Bnew(kkk),'Georgios') && contains(Bnew(kkk),matchstr)
                     c=char(Bnew(kkk));
                     dt(kkk)=datetime(c(:,61:70),'InputFormat','yyyy-MM-dd');
                     newStr(kkk) = erase(Cnew(kkk),matchstr2);
                     student{kkk}=matchstr6;
                     Monkey{kkk}=NHP1;
                 elseif contains(Bnew(kkk),'Michelle') && contains(Bnew(kkk),matchstr)
                     c=char(Bnew(kkk));
                     dt(kkk)=datetime(c(:,61:70),'InputFormat','yyyy-MM-dd');
    %                      pt(kkk)=datetime(c(:,61:70),'ConvertFrom','yyyy-MM-dd');
                     newStr(kkk) = erase(Cnew(kkk),matchstr2);  
                     student{kkk}=matchstr9;
                     Monkey{kkk}=NHP1;
                 elseif contains(Bnew(kkk),'Elizabeth') && contains(Bnew(kkk),matchstr)
                     c=char(Bnew(kkk));
                     dt(kkk)=datetime(c(:,62:71),'InputFormat','yyyy-MM-dd');
                     newStr(kkk) = erase(Cnew(kkk),matchstr2);
                     student{kkk}=matchstr7;  
                     Monkey{kkk}=NHP1;  
                 elseif contains(Bnew(kkk),'Cooper') && contains(Bnew(kkk),matchstr)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,59:68),'InputFormat','yyyy-MM-dd');  
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr8;    
                    Monkey{kkk}=NHP1;        
                  elseif contains(Bnew(kkk),'Victoria') && contains(Bnew(kkk),matchstr)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,61:70),'InputFormat','yyyy-MM-dd');
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr10;  
                    Monkey{kkk}=NHP1;  
                elseif  contains(Bnew(kkk),'Nick') && contains(Bnew(kkk),matchstr)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,57:66),'InputFormat','yyyy-MM-dd');
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr11;
                    Monkey{kkk}=NHP1;  
                elseif  contains(Bnew(kkk),'Katie') && contains(Bnew(kkk),matchstr1)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,51:60),'InputFormat','yyyy-MM-dd');
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr3;
                    Monkey{kkk}=NHP2;        
                elseif contains(Bnew(kkk),'Pip') && contains(Bnew(kkk),matchstr1)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,48:57),'InputFormat','yyyy-MM-dd');
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr4;  
                    Monkey{kkk}=NHP2;      
                 elseif contains(Bnew(kkk),'Cristian') && contains(Bnew(kkk),matchstr1)
                     c=char(Bnew(kkk));
                     dt(kkk)=datetime(c(:,53:62),'InputFormat','yyyy-MM-dd');  
                     newStr(kkk) = erase(Cnew(kkk),matchstr2);
                     student{kkk}=matchstr5;
                     Monkey{kkk}=NHP2;
                 elseif contains(Bnew(kkk),'Georgios') && contains(Bnew(kkk),matchstr1)
                     c=char(Bnew(kkk));
                     dt(kkk)=datetime(c(:,53:62),'InputFormat','yyyy-MM-dd');
                     newStr(kkk) = erase(Cnew(kkk),matchstr2);
                     student{kkk}=matchstr6;
                     Monkey{kkk}=NHP2;
                 elseif contains(Bnew(kkk),'Michelle') && contains(Bnew(kkk),matchstr1)
                     c=char(Bnew(kkk));
                     dt(kkk)=datetime(c(:,53:62),'InputFormat','yyyy-MM-dd');
                     newStr(kkk) = erase(Cnew(kkk),matchstr2);
                     student{kkk}=matchstr9;
                     Monkey{kkk}=NHP2;
                 elseif contains(Bnew(kkk),'Elizabeth') && contains(Bnew(kkk),matchstr1)
                     c=char(Bnew(kkk));
                     dt(kkk)=datetime(c(:,54:63),'InputFormat','yyyy-MM-dd');
                     newStr(kkk) = erase(Cnew(kkk),matchstr2);
                     student{kkk}=matchstr7;
                     Monkey{kkk}=NHP2;      
                 elseif contains(Bnew(kkk),'Cooper') && contains(Bnew(kkk),matchstr1)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,52:61),'InputFormat','yyyy-MM-dd');  
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr8;
                    Monkey{kkk}=NHP2;        
                  elseif contains(Bnew(kkk),'Victoria') && contains(Bnew(kkk),matchstr1)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,53:62),'InputFormat','yyyy-MM-dd');
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr10;
                    Monkey{kkk}=NHP2;  
                 elseif contains(Bnew(kkk),'Nick') && contains(Bnew(kkk),matchstr1)
                    c=char(Bnew(kkk));
                    dt(kkk)=datetime(c(:,49:58),'InputFormat','yyyy-MM-dd');  
                    newStr(kkk) = erase(Cnew(kkk),matchstr2);
                    student{kkk}=matchstr11;
                    Monkey{kkk}=NHP2;          
                end
            end
    
            ApAv_dates=dt';
            ApAv_cluster_dates=rmmissing(ApAv_dates);
            ApAv_cluster_dates.Format= 'M/d/yyyy';
            ApAv_cluster_dates_cat = categorical(cellstr(ApAv_cluster_dates));
    
            ApAv_sesnames=newStr';
            ApAv_cluster_names=rmmissing(ApAv_sesnames);
    
            ApAv_students=(student(~cellfun('isempty',student)))';
    
            ApAv_monkeys=(Monkey(~cellfun('isempty',Monkey)))';
            
            newData8b=cellstr(ApAv_cluster_dates_cat); %datestr(ApAv_cluster_dates,formatOut);
    
    
            tripple_info=[ApAv_cluster_names,newData8b,ApAv_students,ApAv_monkeys];
            
            malaks=tripple_info;
            
            str2='_in2';
            str3='_in3';
            str4='_in4';
            str5='_in5';
            str6='_in6';
            str7='_in1';
    
               newData=tripple_info;
               pousti2=strfind(newData(:,1),[str2]);
               pousti2a=find(~cellfun(@isempty,pousti2));
           
          for mmm=1:size(newData,1)
                   malaks=newData;
                   pousti2=strfind(newData(:,1),[str2]);
                   pousti2a=find(~cellfun(@isempty,pousti2));
                 if size(pousti2a,1)>0
                   jjj1=repmat(malaks(pousti2a(1),:),2,1);
                   jjj2= strrep( jjj1 , '_in2' , '_inTWO');
                   newRows2 = jjj2;
                   newData = [malaks((1:pousti2a(1)-1), :); newRows2; malaks((pousti2a(1)+1):end, :)];
                  end
          end

               newData1=newData;
               pousti3=strfind(newData1(:,1),[str3]);
               pousti3a=find(~cellfun(@isempty,pousti3));

            for mmm=1:size(newData1,1)
                   malaks=newData1;
                   pousti3=strfind(newData1(:,1),[str3]);
                   pousti3a=find(~cellfun(@isempty,pousti3));
                 if size(pousti3a,1)>0
                   jjj1=repmat(malaks(pousti3a(1),:),3,1);
                   jjj2= strrep( jjj1 , '_in3' , '_inTHREE');
                   newRows3 = jjj2;  
                   newData1 = [malaks((1:pousti3a(1)-1), :); newRows3; malaks((pousti3a(1)+1):end, :)];
                  end
            end    

               newData2=newData1;
               pousti4=strfind(newData2(:,1),[str4]);
               pousti4a=find(~cellfun(@isempty,pousti4));
    
    
            for mmm=1:size(newData2,1)
                   malaks=newData2;
                   pousti4=strfind(newData2(:,1),[str4]);
                   pousti4a=find(~cellfun(@isempty,pousti4));
                 if size(pousti4a,1)>0
                   jjj1=repmat(malaks(pousti4a(1),:),4,1);
                   jjj2= strrep( jjj1 , '_in4' , '_inFOUR');
                   newRows4 = jjj2;  
                   newData2 = [malaks((1:pousti4a(1)-1), :); newRows4; malaks((pousti4a(1)+1):end, :)];
                  end
            end      
    
        newData3 = strrep(newData2,'_inTHREE','_in3');
        newData4 = strrep(newData3,'_inTWO','_in2');
        newData5 = strrep(newData4,'_inFOUR','_in4');
    
               newData_number=tripple_info;
               pousti2=strfind(newData_number(:,1),[str2]);
               pousti2a=find(~cellfun(@isempty,pousti2));
           
          for mmm=1:size(newData_number,1)
                   malaks=newData_number;
                   pousti2=strfind(newData_number(:,1),[str2]);
                   pousti2a=find(~cellfun(@isempty,pousti2));
                 if size(pousti2a,1)>0
                   jjj1=repmat(malaks(pousti2a(1),:),2,1);
                   jjj1_number=jjj1;
                   jjj1_number(:,1)=[1 2];
                   newRows2_number=jjj1_number;
                   newData_number = [malaks((1:pousti2a(1)-1), :); newRows2_number; malaks((pousti2a(1)+1):end, :)];
                  end
          end
    
               newData_number1=newData_number;
               pousti3=strfind(newData_number1(:,1),[str3]);
               pousti3a=find(~cellfun(@isempty,pousti3));
    
    
            for mmm=1:size(newData_number1,1)
                   malaks=newData_number1;
                   pousti3=strfind(newData_number1(:,1),[str3]);
                   pousti3a=find(~cellfun(@isempty,pousti3));
                 if size(pousti3a,1)>0
                   jjj1=repmat(malaks(pousti3a(1),:),3,1);
                   jjj1_number=jjj1;
                   jjj1_number(:,1)=[1 2 3];
                   newRows2_number=jjj1_number;
                   newData_number1 = [malaks((1:pousti3a(1)-1), :); newRows2_number; malaks((pousti3a(1)+1):end, :)];
                  end
            end    
    
               newData_number2=newData_number1;
               pousti4=strfind(newData_number2(:,1),[str4]);
               pousti4a=find(~cellfun(@isempty,pousti4));
    
    
            for mmm=1:size(newData_number2,1)
                   malaks=newData_number2;
                   pousti4=strfind(newData_number2(:,1),[str4]);
                   pousti4a=find(~cellfun(@isempty,pousti4));
                 if size(pousti4a,1)>0
                   jjj1=repmat(malaks(pousti4a(1),:),4,1);
                   jjj1_number=jjj1;
                   jjj1_number(:,1)=[1 2 3 4];
                   newRows2_number=jjj1_number;
                   newData_number2 = [malaks((1:pousti4a(1)-1), :); newRows2_number; malaks((pousti4a(1)+1):end, :)];
                  end
            end  

               newData_number3=newData_number2;
               pousti4=strfind(newData_number3(:,1),[str4]);
               pousti4a=find(~cellfun(@isempty,pousti4));

            for mmm=1:size(newData_number3,1)
                   malaks=newData_number3;
                   pousti5=strfind(newData_number3(:,1),[str7]);
                   pousti5a=find(~cellfun(@isempty,pousti5));
                 if size(pousti5a,1)>0
                   jjj1=repmat(malaks(pousti5a(1),:),1,1);
                   jjj1_number=jjj1;
                   jjj1_number(:,1)=[1];
                   newRows2_number=jjj1_number;
                   newData_number3 = [malaks((1:pousti5a(1)-1), :); newRows2_number; malaks((pousti5a(1)+1):end, :)];
                  end
            end    
    
        newData3 = strrep(newData2,'_inTHREE','_in3');
        newData4 = strrep(newData3,'_inTWO','_in2');
        newData5 = strrep(newData4,'_inFOUR','_in4');
        newData6=[newData5,newData_number3(:,1)];
        newData7=[newData6(:,2),newData6(:,1),strtrim(newData6(:,5)),newData6(:,4),newData6(:,3)];
        fromExcel=[strtrim(everything(:,1)),strtrim(everything(:,2)),everything(:,3),everything(:,19),everything(:,17),everything(:,18)];
        
    
        for uuu=1:size(newData7,1)
            newData8(uuu,:)={newData7{uuu,:}};
        end
    
         newData8a=strtrim(newData8(:,1));
         newData8b=datestr(newData8a(:,1),formatOut);
    
            asere=num2str(cell2mat(everything(:,3)));
            asere1={asere};
    
            for zzz=1:size(asere(:,1),1)
                    asere1(zzz,1)={asere(zzz,1)};
            end
    
            fromExcel1=[fromExcel(:,1),fromExcel(:,2),asere1,fromExcel(:,4),fromExcel(:,5),fromExcel(:,6)];
    
            AB=newData8(:,1:4);
            BC=fromExcel1(:,1:4);
    
            for i=1:size(AB,1)
                big(:,i)={strjoin(AB(i,:))};
            end
            big1=big';
    
            for j=1:size(BC,1)
                small(:,j)={strjoin(BC(j,:))};
            end
            small1=small';
            
            [Lia, Locb] = ismember(small1,big1); %, 'rows');
    
            Lia1=Lia;
            Lia1(Lia1==0)=[];
    
            Locb1=Locb;
            Locb1(Locb1==0)=[];
            
            newData9=newData8;
    
            for i=1:size(Locb1,1)
                for kk=1:size(fromExcel1,1)
    %                 if Locb>0
                        newData9(Locb1(i),6)=fromExcel1(i,5);
                        newData9(Locb1(i),7)=fromExcel1(i,6);
    %                 end
                end
            end
    
        Spikes_per_bin=(Spikes_per_bin)';
        TrialNumber_idx=find([Events_TTL_value] == 211 | [Events_TTL_value] == 215 | [Events_TTL_value] ==205 | [Events_TTL_value] ==204 | [Events_TTL_value] ==233);
        TrialNumber=numel(TrialNumber_idx);
        
        % outcome events onset in Pavlovian trials
        Outcome_f_idx=find([Events_TTL_value] ==217 |  [Events_TTL_value] ==238);
        Outcome_number=numel(Outcome_f_idx);

        % turns fixation point off, and displays airpuff and reward cue
        Before_outcome_f_idx=find([Events_TTL_value] ==215 |  [Events_TTL_value] ==211);
        % actually here the updated variable aligns with the outcome onset
        % period
        Before_outcome_f_idx=Before_outcome_f_idx+2;
        % time that whichever of the two outcomes is being delivered
        Timing_outcome_f_v2=Events_timestamps(Before_outcome_f_idx,:);
        % timestamps of the two outcomes going to be delivered
        Timestamps_outcome_f_v2=Events_TTL_value(Before_outcome_f_idx,:);

        % Find indices where the value is 214; for red and yellow bar sizes
        indicesOf214 = find(Events_TTL_value == 214);
        
        % For Red Bar Size (2 places before 214)
        indicesOfRedBarSize = indicesOf214 - 2;
        % Ensure indices are within bounds (greater than 0)
        indicesOfRedBarSize = indicesOfRedBarSize(indicesOfRedBarSize > 0);
        % Retrieve the Red_bar_size values
        Red_bar_size = Events_TTL_value(indicesOfRedBarSize);
        
        % For Yellow Bar Size (1 place before 214)
        indicesOfYellowBarSize = indicesOf214 - 1;
        % Ensure indices are within bounds (greater than 0)
        indicesOfYellowBarSize = indicesOfYellowBarSize(indicesOfYellowBarSize > 0);
        % Retrieve the Yellow_bar_size values
        Yellow_bar_size = Events_TTL_value(indicesOfYellowBarSize);

        Triatypes=[];
        Triatypes_v2=[];
        Triatypes_v3=[];
        Triatypes_v4=[];

        % Assuming Events_all is your matrix with size 1234x2
        % Find rows where any column contains the value 1
        rowsWithOne = any(Events_all == 1, 2);
        
        % Remove those rows from the matrix
        Events_all(rowsWithOne , :) = [];

        Triatypes=Events_all;
        Triatypes_v2=Events_all;
        Triatypes_v3=Events_all;
        Triatypes_v4=Events_all;

        Outcome_f=Events_all;

        % Filter Events_all to create Triatypes_v2
        % Values of interest for the second column
        valuesOfInterest = [212, 234, 216, 204, 205];
        
        % Logical index for rows with the second column matching the values of interest
        rowsToKeep = ismember(Events_all(:, 2), valuesOfInterest);
        
        % Keep only the rows that match the criteria
        Triatypes = Events_all(rowsToKeep, :);

        % Filter Events_all to create Triatypes_v2
        % Values of interest for the second column
        valuesOfInterest_v2 = [211, 233, 215, 204, 202, 205];
        
        % Logical index for rows with the second column matching the values of interest
        rowsToKeep_v2 = ismember(Events_all(:, 2), valuesOfInterest_v2);
        
        % Keep only the rows that match the criteria
        Triatypes_v2 = Events_all(rowsToKeep_v2, :);

        % Filter Events_all to create Triatypes_v3
        % Values of interest for the second column
        valuesOfInterest_v3 = [211, 233, 228, 204, 202, 205];
        
        % Logical index for rows with the second column matching the values of interest
        rowsToKeep_v3 = ismember(Events_all(:, 2), valuesOfInterest_v3);
        
        % Keep only the rows that match the criteria
        Triatypes_v3 = Events_all(rowsToKeep_v3, :);

        % Filter Events_all to create Triatypes_v4
        % Values of interest for the second column
        valuesOfInterest_v4 = [251, 252, 254];
        
        % Logical index for rows with the second column matching the values of interest
        rowsToKeep_v4 = ismember(Events_all(:, 2), valuesOfInterest_v4);
        
        % Keep only the rows that match the criteria
        Triatypes_v4 = Events_all(rowsToKeep_v4, :);
 
        for i=1:size(Events_all(:,2),1)
            if (Events_all(i,2)~=227 && Events_all(i,2)~=228)
                Outcome_f(i,:)=NaN;
            end
        end
        
        Outcome_f(any(isnan(Outcome_f), 2), :) = [];

if length(Triatypes) == length(Red_bar_size) && length(Triatypes) == length(Yellow_bar_size)
    %  for detailed subplot
    if length(Triatypes)==1200
        range1=([0:400:length(Triatypes)])';
    elseif length(Triatypes)==1350
        range1=([0:450:length(Triatypes)])';
    elseif length(Triatypes)>1200 && length(Triatypes)<1350
        range1=([0,450,900,length(Triatypes)])';
    elseif length(Triatypes)<1200 && length(Triatypes)>1100
        range1=([0,400,800,length(Triatypes)])';
    elseif length(Triatypes)<1100 && length(Triatypes)>1000
        range1=([0,350,700,length(Triatypes)])';    
    elseif length(Triatypes)<1000 && length(Triatypes)>900
        range1=([0,300,600,length(Triatypes)])';  
    elseif length(Triatypes)<900 && length(Triatypes)>800
        range1=([0,300,600,length(Triatypes)])';  
    elseif length(Triatypes)<800 && length(Triatypes)>700
        range1=([0,250,500,length(Triatypes)])';    
    elseif length(Triatypes)<700 && length(Triatypes)>600
        range1=([0,200,400,length(Triatypes)])';      
    elseif length(Triatypes)<600 && length(Triatypes)>500
        range1=([0,150,300,length(Triatypes)])';      
    elseif length(Triatypes)<500 && length(Triatypes)>400
        range1=([0,150,300,length(Triatypes)])';   
    elseif length(Triatypes)<400 && length(Triatypes)>300
        range1=([0,100,200,length(Triatypes)])'; 
    elseif length(Triatypes)<300 && length(Triatypes)>200
        range1=([0,50,100,length(Triatypes)])'; 
    elseif length(Triatypes)<200 && length(Triatypes)>100
        range1=([0,50,100,length(Triatypes)])'; 
    elseif length(Triatypes)==800 
        range1=([0,250,500,800])';
    elseif length(Triatypes)==1000 
        range1=([0,350,700,1000])';
    elseif length(Triatypes)==1100 
        range1=([0,350,700,length(Triatypes)])';       
    end
elseif length(Triatypes_v2) == length(Red_bar_size) && length(Triatypes) == length(Yellow_bar_size)
    if length(Triatypes_v2)==1200
        range1=([0:400:length(Triatypes_v2)])';
    elseif length(Triatypes_v2)==1350
        range1=([0:450:length(Triatypes_v2)])';
    elseif length(Triatypes_v2)>1200 && length(Triatypes_v2)<1350
        range1=([0,450,900,length(Triatypes_v2)])';
    elseif length(Triatypes_v2)<1200 && length(Triatypes_v2)>1100
        range1=([0,400,800,length(Triatypes_v2)])';
    elseif length(Triatypes_v2)<1100 && length(Triatypes_v2)>1000
        range1=([0,350,700,length(Triatypes_v2)])';    
    elseif length(Triatypes_v2)<1000 && length(Triatypes_v2)>900
        range1=([0,300,600,length(Triatypes_v2)])';  
    elseif length(Triatypes_v2)<900 && length(Triatypes_v2)>800
        range1=([0,300,600,length(Triatypes_v2)])';  
    elseif length(Triatypes_v2)<800 && length(Triatypes_v2)>700
        range1=([0,250,500,length(Triatypes_v2)])';  
    elseif length(Triatypes_v2)<700 && length(Triatypes_v2)>600
        range1=([0,200,400,length(Triatypes_v2)])';  
    elseif length(Triatypes_v2)<600 && length(Triatypes_v2)>500
        range1=([0,150,300,length(Triatypes_v2)])';      
    elseif length(Triatypes_v2)<500 && length(Triatypes_v2)>400
        range1=([0,150,300,length(Triatypes_v2)])';   
    elseif length(Triatypes_v2)<400 && length(Triatypes_v2)>300
        range1=([0,100,200,length(Triatypes_v2)])'; 
    elseif length(Triatypes_v2)<300 && length(Triatypes_v2)>200
        range1=([0,50,100,length(Triatypes_v2)])'; 
    elseif length(Triatypes_v2)<200 && length(Triatypes_v2)>100
        range1=([0,50,100,length(Triatypes_v2)])';   
    elseif length(Triatypes_v2)==800 
        range1=([0,250,500,800])';
    elseif length(Triatypes_v2)==1000 
        range1=([0,350,700,1000])';
    elseif length(Triatypes_v2)==1100 
        range1=([0,350,700,length(Triatypes_v2)])';       
    end
elseif length(Triatypes_v3) == length(Red_bar_size) && length(Triatypes_v3) == length(Yellow_bar_size)
    if length(Triatypes_v3)==1200
        range1=([0:400:length(Triatypes_v3)])';
    elseif length(Triatypes_v3)==1350
        range1=([0:450:length(Triatypes_v3)])';
    elseif length(Triatypes_v3)>1200 && length(Triatypes_v3)<1350
        range1=([0,450,900,length(Triatypes_v3)])';
    elseif length(Triatypes_v3)<1200 && length(Triatypes_v3)>1100
        range1=([0,400,800,length(Triatypes_v3)])';
    elseif length(Triatypes_v3)<1100 && length(Triatypes_v3)>1000
        range1=([0,350,700,length(Triatypes_v3)])';    
    elseif length(Triatypes_v3)<1000 && length(Triatypes_v3)>900
        range1=([0,300,600,length(Triatypes_v3)])';  
    elseif length(Triatypes_v3)<900 && length(Triatypes_v3)>800
        range1=([0,300,600,length(Triatypes_v3)])';  
    elseif length(Triatypes_v3)<800 && length(Triatypes_v3)>700
        range1=([0,250,500,length(Triatypes_v3)])';     
    elseif length(Triatypes_v3)<700 && length(Triatypes_v3)>600
        range1=([0,200,500,length(Triatypes_v3)])';    
    elseif length(Triatypes_v3)<600 && length(Triatypes_v3)>500
        range1=([0,150,300,length(Triatypes_v3)])';      
    elseif length(Triatypes_v3)<500 && length(Triatypes_v3)>400
        range1=([0,150,300,length(Triatypes_v3)])';   
    elseif length(Triatypes_v3)<400 && length(Triatypes_v3)>300
        range1=([0,100,200,length(Triatypes_v3)])'; 
    elseif length(Triatypes_v3)<300 && length(Triatypes_v3)>200
        range1=([0,50,100,length(Triatypes_v3)])'; 
    elseif length(Triatypes_v3)<200 && length(Triatypes_v3)>100
        range1=([0,50,100,length(Triatypes_v3)])';         
    elseif length(Triatypes_v3)==800 
        range1=([0,250,500,800])';
    elseif length(Triatypes_v3)==1000 
        range1=([0,350,700,1000])';
    elseif length(Triatypes_v3)==1100 
        range1=([0,350,700,length(Triatypes_v3)])';       
    end  
else
    Red_bar_size = [];
    Yellow_bar_size = [];
    
    for i = 1:length(Events_TTL_value)-2
        if ismember(Events_TTL_value(i), [251, 252, 254])
            Red_bar_size = [Red_bar_size; Events_TTL_value(i+1)];
            Yellow_bar_size = [Yellow_bar_size; Events_TTL_value(i+2)];
        end
    end

    if length(Triatypes_v4)==1200
        range1=([0:400:length(Triatypes_v4)])';
    elseif length(Triatypes_v4)==1350
        range1=([0:450:length(Triatypes_v4)])';
    elseif length(Triatypes_v4)>1350 && length(Triatypes_v4)<1400
        range1=([0,450,900,length(Triatypes_v4)])';
    elseif length(Triatypes_v4)>1200 && length(Triatypes_v4)<1350
        range1=([0,450,900,length(Triatypes_v4)])';
    elseif length(Triatypes_v4)<1200 && length(Triatypes_v4)>1100
        range1=([0,400,800,length(Triatypes_v4)])';
    elseif length(Triatypes_v4)<1100 && length(Triatypes_v4)>1000
        range1=([0,350,700,length(Triatypes_v4)])';    
    elseif length(Triatypes_v4)<1000 && length(Triatypes_v4)>900
        range1=([0,300,600,length(Triatypes_v4)])';  
    elseif length(Triatypes_v4)<900 && length(Triatypes_v4)>800
        range1=([0,300,600,length(Triatypes_v4)])';  
    elseif length(Triatypes_v4)<800 && length(Triatypes_v4)>700
        range1=([0,250,500,length(Triatypes_v4)])';    
    elseif length(Triatypes_v4)<700 && length(Triatypes_v4)>600
        range1=([0,200,500,length(Triatypes_v4)])';    
    elseif length(Triatypes_v4)<600 && length(Triatypes_v4)>500
        range1=([0,150,300,length(Triatypes_v4)])';      
    elseif length(Triatypes_v4)<500 && length(Triatypes_v4)>400
        range1=([0,150,300,length(Triatypes_v4)])';   
    elseif length(Triatypes_v4)<400 && length(Triatypes_v4)>300
        range1=([0,100,200,length(Triatypes_v4)])'; 
    elseif length(Triatypes_v4)<300 && length(Triatypes_v4)>200
        range1=([0,50,100,length(Triatypes_v4)])'; 
    elseif length(Triatypes_v4)<200 && length(Triatypes_v4)>100
        range1=([0,50,100,length(Triatypes_v4)])';              
    elseif length(Triatypes_v4)==800 
        range1=([0,250,500,800])';
    elseif length(Triatypes_v4)==1000 
        range1=([0,350,700,1000])';
    elseif length(Triatypes_v4)==1100 
        range1=([0,350,700,length(Triatypes_v4)])';       
    end         

end

    range2=range1+1;
    bothends=[range2,range1];

    % fixation breaks during the cue period
% Initial, first alternative, and second alternative codes for filtering
if length(Triatypes) == length(Red_bar_size) && length(Triatypes) == length(Yellow_bar_size)
    validCodes = [212, 201, 207, 206, 216, 205, 204, 202, 236];
elseif length(Triatypes_v2) == length(Red_bar_size) && length(Triatypes_v2) == length(Yellow_bar_size)
    validCodes = [211, 239, 244, 206, 215, 205, 204, 202, 236];
elseif length(Triatypes_v3) == length(Red_bar_size) && length(Triatypes_v3) == length(Yellow_bar_size)    
    validCodes = [211, 239, 244, 206, 228, 205, 204, 202, 236]; % Directly using these codes
elseif length(Triatypes_v4) == length(Red_bar_size) && length(Triatypes_v4) == length(Yellow_bar_size)    
    validCodes = [211, 239, 244, 206, 228, 205, 204, 202, 236]; % Directly using these codes    
end

% Filter Events_all based on the initial valid codes
CueBreakTrialCodes = ismember(Events_all(:,2), validCodes);
CueBreakTrials = Events_all(CueBreakTrialCodes, :);
CueBreakTrials = CueBreakTrials(CueBreakTrials(:,2) ~= 0, :); % Remove rows where the second column is 0

% Initial categorization
CueBreakTrials2 = zeros(size(CueBreakTrials, 1), size(CueBreakTrials, 2));
validFinalCodes = [203, 206];
for i = 1:size(CueBreakTrials, 1)
    CueBreakTrials2(i, :) = ismember(CueBreakTrials(i, 2), validFinalCodes);
end

% Check if length conditions are met, if not, re-filter using the first alternative codes
if length(CueBreakTrials2) ~= length(Red_bar_size) || length(CueBreakTrials2) ~= length(Yellow_bar_size)
    CueBreakTrialCodes = ismember(Events_all(:,2), validCodes);
    CueBreakTrials = Events_all(CueBreakTrialCodes, :);
    CueBreakTrials = CueBreakTrials(CueBreakTrials(:,2) ~= 0, :);
    % Recategorize
    CueBreakTrials2 = zeros(size(CueBreakTrials, 1), size(CueBreakTrials, 2));
    for i = 1:size(CueBreakTrials, 1)
        CueBreakTrials2(i, :) = ismember(CueBreakTrials(i, 2), validFinalCodes);
    end
end

% Check again if length conditions are not met, now use the second alternative codes
if length(CueBreakTrials2) ~= length(Red_bar_size) || length(CueBreakTrials2) ~= length(Yellow_bar_size)
    CueBreakTrialCodes = ismember(Events_all(:,2), validCodes);
    CueBreakTrials = Events_all(CueBreakTrialCodes, :);
    CueBreakTrials = CueBreakTrials(CueBreakTrials(:,2) ~= 0, :);
    % Recategorize
    CueBreakTrials2 = zeros(size(CueBreakTrials, 1), size(CueBreakTrials, 2));
    for i = 1:size(CueBreakTrials, 1)
        CueBreakTrials2(i, :) = ismember(CueBreakTrials(i, 2), validFinalCodes);
    end
end

% CueBreakTrials2 now contains 1 where event codes are 203 or 206, and 0 otherwise.


%%%%%%
%%%%%%%%   
%     Breaks_Decision=[];
    
    for k=1:size(range1,1)-1
        Breaks_Decision(x,k)=sum(CueBreakTrials2(range2(k):range1(k+1))==1);
    end
        
        
        % Part 2: Construct Trial_info based on the condition
        % Check if lengths are equal
        if length(Triatypes) == length(Red_bar_size) && length(Triatypes) == length(Yellow_bar_size)
            % Lengths are equal, use Triatypes
            Trial_info = [Triatypes, Red_bar_size, Yellow_bar_size];
        elseif length(Triatypes_v2) == length(Red_bar_size) && length(Triatypes_v2) == length(Yellow_bar_size)
            % Lengths are not equal, use Triatypes_v2
            Trial_info = [Triatypes_v2, Red_bar_size, Yellow_bar_size];
        elseif length(Triatypes_v3) == length(Red_bar_size) && length(Triatypes_v3) == length(Yellow_bar_size)
            Trial_info = [Triatypes_v3, Red_bar_size, Yellow_bar_size];
        else 
            Trial_info = [Triatypes_v4, Red_bar_size, Yellow_bar_size];
        end

        Trials_good=Trial_info;

    % fixation breaks during the fixation cue
% Initialize Precue_fix_breaks to an empty array
Precue_fix_breaks = [];

% Check if lengths of Triatypes, Red_bar_size, and Yellow_bar_size match
if length(Triatypes) == length(Red_bar_size) && length(Triatypes) == length(Yellow_bar_size)
    SelectedTriatypes = Triatypes;
elseif length(Triatypes_v2) == length(Red_bar_size) && length(Triatypes_v2) == length(Yellow_bar_size)
    SelectedTriatypes = Triatypes_v2;
elseif length(Triatypes_v3) == length(Red_bar_size) && length(Triatypes_v3) == length(Yellow_bar_size)
    SelectedTriatypes = Triatypes_v3;
else
    SelectedTriatypes = Triatypes_v4;
end

% Vectorized operation to set values based on the condition
% This replaces the entire for-loop with a single, efficient line.
Precue_fix_breaks = zeros(size(SelectedTriatypes, 1), size(SelectedTriatypes, 2));
Precue_fix_breaks(SelectedTriatypes(:,2) == 204 | SelectedTriatypes(:,2) == 205, :) = 1;

%     Breaks_Fix=[];
    for k=1:size(range1,1)-1
        Breaks_Fix(x,k)=sum(Precue_fix_breaks(range2(k):range1(k+1))==1);
    end        
% Assuming Events_all is already defined as a 15060x2 matrix
% Initialize an empty matrix for the new data
new_matrix = [];

% Iterate through each row in Events_all
for i = 1:size(Events_all, 1)
    % Check if the current event is 227 or 228
    if Events_all(i, 2) == 227 || Events_all(i, 2) == 228
        % Ensure there are at least 6 events before
        if i > 6
            % Extract the required information
            timing_3_before = Events_all(i-3, 1); % The timing 3 events before the 227 or 228 event
            event_code = Events_all(i, 2); % The event code
            event_code_6_before = Events_all(i-6, 2); % Event code 6 positions before
            event_code_5_before = Events_all(i-5, 2); % Event code 5 positions before
            
            % Append to the new matrix
            new_matrix = [new_matrix; timing_3_before, event_code, event_code_6_before, event_code_5_before];
        end
    end
end

        % finding the bar sizes when error trials
        Precue_fix_breaks1=[Precue_fix_breaks,Red_bar_size,Yellow_bar_size];
%         CueBreakTrials3=[CueBreakTrials2,Red_bar_size,Yellow_bar_size];
        
        Trials_good2=[new_matrix,Outcome_f];
        
        Rew_30_CO=Trials_good2(Trials_good2(:,3)==30);
        Rew_60_CO=Trials_good2(Trials_good2(:,3)==60);
        Rew_90_CO=Trials_good2(Trials_good2(:,3)==90);
        Rew_120_CO=Trials_good2(Trials_good2(:,3)==120);
        Rew_150_CO=Trials_good2(Trials_good2(:,3)==150);
        Rew_180_CO=Trials_good2(Trials_good2(:,3)==180);
        
        Rew_30_60_CO=sort([Rew_30_CO;Rew_60_CO]);
        Rew_90_120_CO=sort([Rew_90_CO;Rew_120_CO]);
        Rew_150_180_CO=sort([Rew_150_CO;Rew_180_CO]);
        
        Rew_30_60_90_CO=sort([Rew_30_CO;Rew_60_CO;Rew_90_CO]);
        Rew_120_150_180_CO=sort([Rew_120_CO;Rew_150_CO;Rew_180_CO]);
        
        Rew_30_RD=Trials_good2(Trials_good2(:,3)==30,5);
        Rew_60_RD=Trials_good2(Trials_good2(:,3)==60,5);
        Rew_90_RD=Trials_good2(Trials_good2(:,3)==90,5);
        Rew_120_RD=Trials_good2(Trials_good2(:,3)==120,5);
        Rew_150_RD=Trials_good2(Trials_good2(:,3)==150,5);
        Rew_180_RD=Trials_good2(Trials_good2(:,3)==180,5);
        
        Rew_30_60_RD=sort([Rew_30_RD;Rew_60_RD]);
        Rew_90_120_RD=sort([Rew_90_RD;Rew_120_RD]);
        Rew_150_180_RD=sort([Rew_150_RD;Rew_180_RD]);
        
        Rew_30_60_90_RD=sort([Rew_30_RD;Rew_60_RD;Rew_90_RD]);
        Rew_120_150_180_RD=sort([Rew_120_RD;Rew_150_RD;Rew_180_RD]);
        
        Air_30_CO=Trials_good2(Trials_good2(:,4)==30);
        Air_60_CO=Trials_good2(Trials_good2(:,4)==60);
        Air_90_CO=Trials_good2(Trials_good2(:,4)==90);
        Air_120_CO=Trials_good2(Trials_good2(:,4)==120);
        Air_150_CO=Trials_good2(Trials_good2(:,4)==150);
        Air_180_CO=Trials_good2(Trials_good2(:,4)==180);
        
        Air_30_60_CO=sort([Air_30_CO;Air_60_CO]);
        Air_90_120_CO=sort([Air_90_CO;Air_120_CO]);
        Air_150_180_CO=sort([Air_150_CO;Air_180_CO]);
        
        Air_30_60_90_CO=sort([Air_30_CO;Air_60_CO;Air_90_CO]);
        Air_120_150_180_CO=sort([Air_120_CO;Air_150_CO;Air_180_CO]);

        Air_30_RD=Trials_good2(Trials_good2(:,4)==30,5);
        Air_60_RD=Trials_good2(Trials_good2(:,4)==60,5);
        Air_90_RD=Trials_good2(Trials_good2(:,4)==90,5);
        Air_120_RD=Trials_good2(Trials_good2(:,4)==120,5);
        Air_150_RD=Trials_good2(Trials_good2(:,4)==150,5);
        Air_180_RD=Trials_good2(Trials_good2(:,4)==180,5);
        
        Air_30_60_RD=sort([Air_30_RD;Air_60_RD]);
        Air_90_120_RD=sort([Air_90_RD;Air_120_RD]);
        Air_150_180_RD=sort([Air_150_RD;Air_180_RD]);
        
        Air_30_60_90_RD=sort([Air_30_RD;Air_60_RD;Air_90_RD]);
        Air_120_150_180_RD=sort([Air_120_RD;Air_150_RD;Air_180_RD]);
        
        % plotting the forced reward/punishment trials
        for s=1:44
       
        ev_f={Rew_30_CO,Rew_60_CO,Rew_90_CO,Rew_120_CO,Rew_150_CO,Rew_180_CO,...
            Air_30_CO,Air_60_CO,Air_90_CO,Air_120_CO,Air_150_CO,Air_180_CO,...
            Rew_30_RD,Rew_60_RD,Rew_90_RD,Rew_120_RD,Rew_150_RD,Rew_180_RD,...
            Air_30_RD,Air_60_RD,Air_90_RD,Air_120_RD,Air_150_RD,Air_180_RD,...
            Rew_30_60_CO,Rew_90_120_CO,Rew_150_180_CO,...
            Rew_30_60_RD,Rew_90_120_RD,Rew_150_180_RD,...
            Air_30_60_CO,Air_90_120_CO,Air_150_180_CO,...
            Air_30_60_RD,Air_90_120_RD,Air_150_180_RD,...
            Rew_30_60_90_CO,Rew_120_150_180_CO,...
            Air_30_60_90_CO,Air_120_150_180_CO,...
            Rew_30_60_90_RD,Rew_120_150_180_RD,...
            Air_30_60_90_RD,Air_120_150_180_RD};

                        % Skip the loop iteration if the current cell is empty
            if isempty(ev_f{1, s})
                continue;  % Skip to the next iteration of the loop
            end

        Timestamps_minus_f=(ev_f{1,s}-2)';
        Timestamps_plus_f=(ev_f{1,s}+2)';
        Timestamps_zero_f=(ev_f{1,s}+0)';
        
        z=[];
        % spike_index=[];
        
        for z=1:size(Timestamps_plus_f,2)
            spike_index_f{z} = find(Spikes>Timestamps_minus_f(z) & Spikes<Timestamps_plus_f(z));
            logicalIndexes_f{z} = (Spikes>Timestamps_minus_f(z) & Spikes<Timestamps_plus_f(z));
            linearIndexes_f{z} = find(logicalIndexes_f{z});
        end
        
        k=[];
        
        for k=1:size(spike_index_f,2)
            vector_size_f(k)= size(spike_index_f{1,k},2);
        end
        
        Spikes_num_f=sum(vector_size_f);
        
        full_range_f=2:-0.05:-2;
        
        j=[];
        full_bin_f=[];
        
        for j=1:size(Timestamps_minus_f,2)
            full_bin_f(j,:)=Timestamps_zero_f(1,j)-full_range_f;
        end
        
        u=[];
        p=[];
        
        spike_index_full_f=[];
        logicalIndexes_full_f=[];
        linearIndexes_full_f=[];
        
        for u=1:size(full_bin_f,1)
            for p=1:size(full_bin_f,2)-1
                spike_index_full_f{u,p}=find(Spikes>full_bin_f(u,p) & Spikes<full_bin_f(u,p+1));
                logicalIndexes_full_f = (Spikes>full_bin_f(u,p) & Spikes<full_bin_f(u,p+1));
                linearIndexes_full_f{u,p} = find(logicalIndexes_full_f);
            end
        end
        
        k=[];
        n=[];
        vector_size_full_f=[];
        
        for k=1:size(spike_index_full_f,1)
            for n=1:size(spike_index_full_f,2)
                vector_size_full_f(k,n)= size(spike_index_full_f{k,n},2);
            end
        end
        
        Partial_sum_f=sum(vector_size_full_f);
        Spikes_per_bin_f(s,:)=Partial_sum_f;
        Spikes_full_sum_f{s}=sum(sum(vector_size_full_f));
        Partial_sum_f=[];
        
        end
        
        event_names_f={'RewCO=30','RewCO=60','RewCO=90','RewCO=120','RewCO=150','RewCO=180',...
            'AirCO=30','AirCO=60','AirCO=90','AirCO=120','AirCO=150','AirCO=180',...
            'RewRD=30','RewRD=60','RewRD=90','RewRD=120','RewRD=150','RewRD=180',...
            'AirRD=30','AirRD=60','AirRD=90','AirRD=120','AirRD=150','AirRD=180',...
            'Rew 30 60_CO','Rew 90 120 CO','Rew 150 180 CO',...
            'Rew 30 60_RD','Rew 90 120 RD','Rew 150 180 RD',...
            'Air 30 60_CO','Air 90 120 CO','Air 150 180 CO',...
            'Air 30 60_RD','Air 90 120 RD','Air 150 180 RD',...
            'Rew 30 60 90 CO','Rew 120 150 180 CO',...
            'Air 30 60 90 CO','Air 120 150 180 CO',...
            'Rew 30 60 90 RD','Rew 120 150 180 RD',...
            'Air 30 60 90 RD','Air 120 150 180 RD'};

        % Bin duration in seconds (50ms)
        bin_duration_seconds = 0.05;
        
        % Convert spike counts to firing rates (in spikes per second)
        Firing_rates_f{x,ph} = Spikes_per_bin_f / bin_duration_seconds;
            RewC0_30=sum(Spikes_per_bin_f(1,41:70))/numel(Rew_30_CO);
            RewC0_60=sum(Spikes_per_bin_f(2,41:70))/numel(Rew_60_CO);
            RewC0_90=sum(Spikes_per_bin_f(3,41:70))/numel(Rew_90_CO);
            RewC0_120=sum(Spikes_per_bin_f(4,41:70))/numel(Rew_120_CO);
            RewC0_150=sum(Spikes_per_bin_f(5,41:70))/numel(Rew_150_CO);
            RewC0_180=sum(Spikes_per_bin_f(6,41:70))/numel(Rew_180_CO);
        
            RewardCues_store{x,ph}=[RewC0_30,RewC0_60,RewC0_90,RewC0_120,RewC0_150,RewC0_180];
            RewardCues=[RewC0_30,RewC0_60,RewC0_90,RewC0_120,RewC0_150,RewC0_180];

        
            AirC0_30=sum(Spikes_per_bin_f(7,41:70))/numel(Air_30_CO);
            AirC0_60=sum(Spikes_per_bin_f(8,41:70))/numel(Air_60_CO);
            AirC0_90=sum(Spikes_per_bin_f(9,41:70))/numel(Air_90_CO);
            AirC0_120=sum(Spikes_per_bin_f(10,41:70))/numel(Air_120_CO);
            AirC0_150=sum(Spikes_per_bin_f(11,41:70))/numel(Air_150_CO);
            AirC0_180=sum(Spikes_per_bin_f(12,41:70))/numel(Air_180_CO);
        
            AirpuffCues_store{x,ph}=[AirC0_30,AirC0_60,AirC0_90,AirC0_120,AirC0_150,AirC0_180];
            AirpuffCues=[AirC0_30,AirC0_60,AirC0_90,AirC0_120,AirC0_150,AirC0_180];
        
            x12=[30,60,90,120,150,180];
            y1=[RewardCues];
            y2=[AirpuffCues];
        
            RewC0_3060=sum(sum(Spikes_per_bin_f(1:2,41:70)))/(numel(Rew_30_CO)+numel(Rew_60_CO));
            RewC0_90120=sum(sum(Spikes_per_bin_f(3:4,41:70)))/(numel(Rew_90_CO)+numel(Rew_120_CO));
            RewC0_150180=sum(sum(Spikes_per_bin_f(5:6,41:70)))/(numel(Rew_150_CO)+numel(Rew_180_CO));
        
            RewardCues_pairs_store{x,ph}=[RewC0_3060,RewC0_90120,RewC0_150180];
            RewardCues_pairs=[RewC0_3060,RewC0_90120,RewC0_150180];

        
            AirC0_3060=sum(sum(Spikes_per_bin_f(7:8,41:70)))/(numel(Air_30_CO)+numel(Air_60_CO));
            AirC0_90120=sum(sum(Spikes_per_bin_f(9:10,41:70)))/(numel(Air_90_CO)+numel(Air_120_CO));
            AirC0_150180=sum(sum(Spikes_per_bin_f(11:12,41:70)))/(numel(Air_150_CO)+numel(Air_180_CO));
        
            AirpuffCues_pairs_store{x,ph}=[AirC0_3060,AirC0_90120,AirC0_150180];
            AirpuffCues_pairs=[AirC0_3060,AirC0_90120,AirC0_150180];

        
            x2=categorical({'30+60','90+120','150+180'});
            x2 = reordercats(x2,{'30+60','90+120','150+180'});
        
            y3=[RewardCues_pairs];
            y4=[AirpuffCues_pairs];

            x1=categorical({'30','60','90','120','150','180'});
            x1 = reordercats(x1,{'30','60','90','120','150','180'});

            h1=[1,2,3,4,5,6];
            ju=[h1;y1];
            is_exclude1 = ju(end,:) ==0;
            ju( :, is_exclude1 ) = [];
            rew_order=(sortrows(ju',2'))';
            rew_order_single{x,ph}=rew_order(1,:);

            ku=[h1;y2];
            is_exclude2 = ku(end,:) ==0;
            ku( :, is_exclude2 ) = [];
            air_order=(sortrows(ku',2'))';
            air_order_single{x,ph}=air_order(1,:);

            h2=[1,2,3];
            tu=[h2;y3];
            is_exclude3 = tu(end,:) ==0;
            tu( :, is_exclude3 ) = [];
            rew_order_two=(sortrows(tu',2'))';
            rew_order_pair{x,ph}=rew_order_two(1,:);

            zu=[h2;y4];
            is_exclude4 = zu(end,:) ==0;
            zu( :, is_exclude4 ) = [];
            air_order_two=(sortrows(zu',2'))';
            air_order_pair{x,ph}=air_order_two(1,:);

   
      end
end

% load('name_normal.mat')
nexnamelist=name;

numRows = size(nexnamelist, 1);
nexfileslist = cell(numRows * 2, 1);
index = 1;

for row = 1:numRows
    rowElements = nexnamelist(row, :);
    for col = 1:numel(rowElements)
        nexfileslist{index} = rowElements{col};
        index = index + 1;
    end
end

nexfileslist = nexfileslist(~cellfun(@isempty, nexfileslist));

% for the name of the folders
% load('folder_normal.mat')
folderlist=folder;

numRows = size(folderlist, 1);
foldersofnexfileslist = cell(numRows * 2, 1);
index = 1;

for row = 1:numRows
    rowElements = folderlist(row, :);
    for col = 1:numel(rowElements)
        foldersofnexfileslist{index} = rowElements{col};
        index = index + 1;
    end
end

foldersofnexfileslist = foldersofnexfileslist(~cellfun(@isempty, foldersofnexfileslist));

% for the stats list
% load('stats_normal.mat')
statslist=stat_1_v2;

P1 = permute(statslist,[3 2 1]);
result = reshape(permute(P1,[1,3,2]),[],size(P1,2));
brik = result(~all(cellfun(@isempty, result), 2), :);

allInfoTogether=[foldersofnexfileslist,nexfileslist,brik];

% Extract the name of the student from each folder name
pattern = '(?<=\\UROP\\).*?(?=\\Mnks\\)';
cellStrings=allInfoTogether(:,1);
% Initialize an empty cell array to store the substrings
student_folder = cell(size(cellStrings));

% Loop through each cell string and extract the substring
        for i = 1:numel(cellStrings)
            % Extract the substring using regexp
            substring = regexp(cellStrings{i}, pattern, 'match');
            
            % Store the extracted substring in the cell array
            if ~isempty(substring)
                student_folder{i} = substring{1};
            end
        end

% Extract the name of the monkey from each folder name 
pattern = '(?<=\\Mnks\\).*?(?=\\Recordings\\)';
cellStrings=allInfoTogether(:,1);
% Initialize an empty cell array to store the substrings
monkey_name = cell(size(cellStrings));

% Loop through each cell string and extract the substring
        for i = 1:numel(cellStrings)
            % Extract the substring using regexp
            substring = regexp(cellStrings{i}, pattern, 'match');
            
            % Store the extracted substring in the cell array
            if ~isempty(substring)
                monkey_name{i} = substring{1};
            end
        end

% Extract the dates from each folder name
% Regular expression pattern
pattern = '\\(\d{4}-\d{2}-\d{2})';

% Initialize an empty cell array to store the dates
dates = cell(size(cellStrings));

        % Loop through each cell string and extract the date
        for i = 1:numel(cellStrings)
            % Extract the date using regexp
            dateMatch = regexp(cellStrings{i}, pattern, 'tokens');
            
            % Store the extracted date in the cell array
             if ~isempty(dateMatch)
                dateStr = dateMatch{1}{1};  % Extracted date string in 'YYYY-MM-DD' format
                dateVec = datevec(dateStr); % Convert date string to a date vector
                dateStrFormatted = sprintf('%d/%d/%d', dateVec(2), dateVec(3), dateVec(1)); % Format the date as 'M/D/YYYY'
                dates{i} = dateStrFormatted;
             end            
        end

% Count occurrences of each file name
fileNames=allInfoTogether(:,2);

        for i = 1:numel(fileNames)
            fileNames{i} = strrep(fileNames{i}, '.nex', '');
        end

% Replace repeated strings with numbers
% Replace consecutive repetitions with numeric values
strings=allInfoTogether(:,1);

% Create a map to store the counts for each string
countMap = containers.Map('KeyType', 'char', 'ValueType', 'double');

% Replace the strings with numeric values
result = cell(size(strings));
        for i = 1:numel(strings)
            str = strings{i};
            if isKey(countMap, str)
                count = countMap(str) + 1;
            else
                count = 1;
            end
            countMap(str) = count;
            result{i} = num2str(count);
        end

        fromStW=[dates,fileNames,result,monkey_name,allInfoTogether];       
        
        matrix1=fromStW;
        matrix2=fromExcel1;
        % Create a logical matrix to mark the matching rows
        matching_rows = false(size(matrix2, 1), 1);
        
        for i = 1:size(matrix2, 1)
            % Convert the first four cells in each row of matrix2 to a string
            row2_string = strjoin(cellfun(@num2str, matrix2(i, 1:4), 'UniformOutput', false), ' ');
            
            % Search for matching rows in matrix1
            matching_idx = [];
            for j = 1:size(matrix1, 1)
                % Convert the first four cells in each row of matrix1 to a string
                row1_string = strjoin(cellfun(@num2str, matrix1(j, 1:4), 'UniformOutput', false), ' ');
                
                % Check if the row2_string matches row1_string
                if isequal(row2_string, row1_string)
                    matching_idx = j;
                    break;
                end
            end
            
            % If a match is found, append the cells in columns 7 to 51 of matrix1 to matrix2
            if ~isempty(matching_idx)
                matrix2(i, 7:61) = matrix1(matching_idx, 7:61);
                matching_rows(i) = true;
            end
        end
        
        % Replace non-matching rows in matrix2 with NaN
        matrix2(~matching_rows, :) = num2cell(NaN);

        % Find rows that contain at least one non-NaN value
        nonNaNRows = any(cellfun(@(x) any(~isnan(x)), matrix2), 2);
        
        % Remove rows that contain only NaN values
        matrix2 = matrix2(nonNaNRows, :);

        % Find rows with 'MUA' or 'poorunit' in the 5th column
        matching_rows = ismember(string(matrix2(:, 5)), ["MUA", "poorunit"]);
%         matching_rows = ismember(string(matrix2(:, 5)), ["poorunit", "goodunit","goodunit"]);
     
        % Remove matching rows from matrix2
        matrix2 = matrix2(~matching_rows, :); 

        % Convert the elements in column 6 of matrix2 to strings
        strings = cellfun(@num2str, matrix2(:, 6), 'UniformOutput', false);
        
        % Remove repeated strings
        uniqueStrings = unique(strings);
        
        % Get the unique strings and their counts
        stringCounts = zeros(size(uniqueStrings));
        for i = 1:numel(uniqueStrings)
            stringCounts(i) = sum(strcmp(strings, uniqueStrings{i}));
        end
        
        % Create a cell array to store the unique strings and their counts
        num_per_ROI = cell(numel(uniqueStrings), 2);
        for i = 1:numel(uniqueStrings)
            num_per_ROI{i, 1} = uniqueStrings{i};
            num_per_ROI{i, 2} = stringCounts(i);
        end

        data=num_per_ROI;    
        second_column = cell2mat(data(:, 2));  % Convert second column to numeric array
        total_sum = sum(second_column);  % Calculate the sum of values in the second column
        disp(total_sum);  % Display the total sum

        % Define column titles
        columnTitles = {'Date', 'CSC#', 'Cluster#', 'NHP', 'Category', 'Location','201','202','203','204','205','206','207',...
            '208','209','210','211','212','213','214','215','216','217','218','219','220','221','222','223','224','225','226',...
            '227','228','229','230','231','232','233','234','235','236','237','238','239','240','241','242','243','244','245',...
            '246','247','248','249','250','251','252','253','254','255'};

        % Insert the row at the top of matrix2
        matrix2 = [columnTitles; matrix2];

        % Find columns with only one non-NaN value and all other values as NaN
        singleNonNaNValueCols = false(1, size(matrix2, 2));
        for col = 1:size(matrix2, 2)
            colValues = matrix2(:, col);
            if all(cellfun(@(x) isequal(size(x), size(colValues{1})), colValues))  % Check if all cells have the same size
                colValues = cell2mat(colValues);
                if sum(~isnan(colValues)) == 1 & all(isnan(colValues))  % Use element-wise logical operators
                    singleNonNaNValueCols(col) = true;
                end
            end
        end
        
        % Remove the columns from matrix2
        matrix2 = matrix2(:, ~singleNonNaNValueCols);

        % Convert matrix2 to a numeric matrix
        numericMatrix2 = cellfun(@str2double, matrix2);
        
        % Find columns with only one non-NaN value and all other values as NaN
        singleNonNaNValueCols = sum(~isnan(numericMatrix2), 1) == 1 & sum(isnan(numericMatrix2), 1) == size(numericMatrix2, 1) - 1;
        
%         % Remove the columns from matrix2
%         matrix2 = matrix2(:, ~singleNonNaNValueCols);

        %%% to find the significant values and the not significant
        % Loop through the rows (excluding the first row)
        for row = 2:size(matrix2, 1)
            % Loop through the columns (from column 7 until the last column)
            for col = 7:size(matrix2, 2)
                % Check if the cell value is numeric and less than 0.05
                if isnumeric(matrix2{row, col}) || ischar(matrix2{row, col})
                    value = str2double(matrix2{row, col}); % Convert to numeric value
                    if ~isnan(value) && value < 0.01
                        matrix2{row, col} = 'YES'; % Replace with 'significant' string
                    elseif ~isnan(value) && value > 0.01
                        matrix2{row, col} = 'NO'; % Replace with 'significant' string
                    end
                end
            end
        end

                        cOFC_numbers = zeros(2, size(matrix2, 2)-6); % Initialize column counts array
                         % Find the indices of rows where column 6 aligns with 'cOFC'
                        cofc_rows = find(strcmp(matrix2(:, 6), 'cOFC'));
                        
                        % Loop through the columns (from column 7 to the last column)
                        for col = 7:size(matrix2, 2)
                            % Count the number of cells with 'YES' in the current column that align with 'cOFC'
                            cOFC_numbers(1, col-6) = sum(strcmp(matrix2(cofc_rows, col), 'YES'));
                        
                            % Count the number of cells with 'NO' in the current column that align with 'cOFC'
                            cOFC_numbers(2, col-6) = sum(strcmp(matrix2(cofc_rows, col), 'NO'));
                        
                            % Calculate the percentage of 'YES' in the current column that align with 'cOFC'
                            cOFC_numbers(3, col-6) = (cOFC_numbers(1, col-6) / (cOFC_numbers(1, col-6) + cOFC_numbers(2, col-6))) * 100;
                        end
                indices = [30, 20, 22, 44, 26, 28, 8, 15, 31, 10];

% 14:	fixation onset for all type of trials
% 33:	simultaneously turns off fixation point and displays choice cue
% 39:	airpuff delivery onset after choosing approach in choice trials
% 44:	small reward delivery onset for choosing avoidance in choice trials
% 19:	turns off square and turns on circle around cross for choosing approach in choice trials
% 43: 	turning off cross and highlight on around square for choosing avoidance in choice trials
% 15:	simultaneously turns off fixation point and displays red bar
% 11:	simultaneously turns off fixation point and displays yellow bar
% 17: 	reward delivery onset in Pavlovian trials
% 38:	airpuff delivery onset in Pavlovian trials


% Events of interest
Events_of_interest = {'Fix point ON';...
                     'Choice cue ON';...
                     'Ap Air ON';...
                     'Av rew ON';...
                     'Ap highlight';...
                     'Av highlight';...
                     'Red bar ON';...
                     'Yellow bar ON';...
                     'Red bar Air ON';...
                     'Yellow bar Rew ON'};

                values_cOFC = (round(cOFC_numbers(3, indices)))';
                
                % for pACC
                          
                        pACC_numbers = zeros(2, size(matrix2, 2)-6); % Initialize column counts array
                         % Find the indices of rows where column 6 aligns with 'pACC'
                        pACC_rows = find(strcmp(matrix2(:, 6), 'pACC'));
                        
                        % Loop through the columns (from column 7 to the last column)
                        for col = 7:size(matrix2, 2)
                            % Count the number of cells with 'YES' in the current column that align with 'pACC'
                            pACC_numbers(1, col-6) = sum(strcmp(matrix2(pACC_rows, col), 'YES'));
                        
                            % Count the number of cells with 'NO' in the current column that align with 'pACC'
                            pACC_numbers(2, col-6) = sum(strcmp(matrix2(pACC_rows, col), 'NO'));
                        
                            % Calculate the percentage of 'YES' in the current column that align with 'pACC'
                            pACC_numbers(3, col-6) = (pACC_numbers(1, col-6) / (pACC_numbers(1, col-6) + pACC_numbers(2, col-6))) * 100;
                        end
                
                values_pACC = (round(pACC_numbers(3, indices)))';
                % for dpACC
                          
                        dpACC_numbers = zeros(2, size(matrix2, 2)-6); % Initialize column counts array
                         % Find the indices of rows where column 6 aligns with 'dpACC'
                        dpACC_rows = find(strcmp(matrix2(:, 6), 'dpACC'));
                        
                        % Loop through the columns (from column 7 to the last column)
                        for col = 7:size(matrix2, 2)
                            % Count the number of cells with 'YES' in the current column that align with 'dpACC'
                            dpACC_numbers(1, col-6) = sum(strcmp(matrix2(dpACC_rows, col), 'YES'));
                        
                            % Count the number of cells with 'NO' in the current column that align with 'dpACC'
                            dpACC_numbers(2, col-6) = sum(strcmp(matrix2(dpACC_rows, col), 'NO'));
                        
                            % Calculate the percentage of 'YES' in the current column that align with 'dpACC'
                            dpACC_numbers(3, col-6) = (dpACC_numbers(1, col-6) / (dpACC_numbers(1, col-6) + dpACC_numbers(2, col-6))) * 100;
                        end
                values_dpACC = (round(dpACC_numbers(3, indices)))';
                
                % for putamen
                          
                        putamen_numbers = zeros(2, size(matrix2, 2)-6); % Initialize column counts array
                         % Find the indices of rows where column 6 aligns with 'putamen'
                        putamen_rows = find(strcmp(matrix2(:, 6), 'putamen'));
                        
                        % Loop through the columns (from column 7 to the last column)
                        for col = 7:size(matrix2, 2)
                            % Count the number of cells with 'YES' in the current column that align with 'putamen'
                            putamen_numbers(1, col-6) = sum(strcmp(matrix2(putamen_rows, col), 'YES'));
                        
                            % Count the number of cells with 'NO' in the current column that align with 'putamen'
                            putamen_numbers(2, col-6) = sum(strcmp(matrix2(putamen_rows, col), 'NO'));
                        
                            % Calculate the percentage of 'YES' in the current column that align with 'putamen'
                            putamen_numbers(3, col-6) = (putamen_numbers(1, col-6) / (putamen_numbers(1, col-6) + putamen_numbers(2, col-6))) * 100;
                        end   
                values_putamen = (round(putamen_numbers(3, indices)))';
                values_all=[values_cOFC,values_pACC,values_dpACC,values_putamen];  

                        % Concatenate Events_of_interest with values_all
        data_summary = horzcat(Events_of_interest, num2cell(values_all));
        
        % Create a table with appropriate column names
        table_headers = {'Event', 'cOFC', 'pACC', 'dpACC', 'putamen'};
        table_data = [table_headers; data_summary];
% Helper function to get the first timestamp of a trial, or Inf if the index is out of range
function ts = getFirstTimeStamp(trialArray, index)
    if index <= length(trialArray)
        ts = trialArray{index}.TimeStamps(1);
    else
        ts = Inf;
    end
end