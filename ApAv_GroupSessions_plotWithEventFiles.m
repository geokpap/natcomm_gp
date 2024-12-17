clear;
clc;
close all;
currentFigure = 1;
subplotIndex = 1;

% Define the root path of the directory tree to search
rootPath = '/Users/georgiospapageorgiou/Dropbox (MIT)/Analyses/Mnks/LittleDebbie/Recordings/ApAv/RT_bothMonkeys/';

% Generate paths to all the subfolders
allSubFolders = genpath(rootPath);

% Split into a list of folder names
folderList = strsplit(allSubFolders, ':');
folderList = folderList(~cellfun('isempty',folderList)); % Remove any empty cells that may exist

% Loop over all folders
for x = 1:length(folderList)
    thisFolder = folderList{x};
    
    % Skip the rootPath
    if strcmp(thisFolder, rootPath)
        continue; % Skip this iteration
    end
    
    fprintf('Processing folder: %s\n', thisFolder);
    
    % Extract the date and animal name from the folder name
    % Assuming the folder name is in the format: ...YYYY-MM-DD...
    datePattern = '\d{4}-\d{2}-\d{2}';
    tokens = regexp(thisFolder, datePattern, 'match');
    
    if ~isempty(tokens)
        dateString = tokens{1};
        formattedDate = strrep(dateString, '-', ''); % Replace hyphens to get 'YYYYMMDD'
    else
        formattedDate = 'DateNotFound'; % Or handle as appropriate
    end
    
    % Assuming the animal name is static as 'LittleDebbie'
    % animalName = 'LittleDebbie';
    
    fprintf('    Date extracted: %s\n', formattedDate);
    % fprintf('    Animal name: %s\n', animalName);
    
    % Find all NEV files in this folder
    filePattern = fullfile(thisFolder, '*.nev');
    nevFiles = dir(filePattern);

    
    % Loop over all NEV files
    for k = 1:length(nevFiles)
        baseFileName = nevFiles(k).name;
        fullFileName = fullfile(thisFolder, baseFileName);
        
        % Display the current file name and date
        fprintf('    Processing file: %s\n', fullFileName);
        fprintf('    On date: %s\n', formattedDate);
        % fprintf('    For animal: %s\n', animalName);
   

        % Load NEV data
        ncs = read_neuralynx_nev(thisFolder);
      
        % Convert TimeStamps from microseconds to milliseconds
        for i = 1:length(ncs)
            ncs(i).TimeStamp = ncs(i).TimeStamp / 1000;
        end
        
        % Initialize variables for different trial types
        all_trials = {};
        choice_trials = {};
        approach_trials = {};
        avoidance_trials = {};
        current_trial = []; % Initialize current_trial as empty
        
        % Process trials
        for i = 1:length(ncs)
            if ncs(i).TTLValue == 9  % Trial start
                current_trial = struct('TimeStamps', [], 'TTLValues', [], 'Type', '');
            end
        
            % If we are within a trial, keep collecting TimeStamps and TTLValues
            if ~isempty(current_trial)
                current_trial.TimeStamps(end+1) = ncs(i).TimeStamp;
                current_trial.TTLValues(end+1) = ncs(i).TTLValue;
                
                % Check for Choice trial marker
                if ncs(i).TTLValue == 233
                    current_trial.Type = 'Choice';
                elseif ncs(i).TTLValue == 239 && strcmp(current_trial.Type, 'Choice')
                    current_trial.Type = 'Approach';
                elseif ncs(i).TTLValue == 244 && strcmp(current_trial.Type, 'Choice')
                    current_trial.Type = 'Avoidance';
                end
            end
        
            if ncs(i).TTLValue == 18  % Trial end
                % Only save the trial if current_trial has been started properly
                if ~isempty(current_trial) && any(current_trial.TTLValues == 9)
                    all_trials{end+1} = current_trial; % Save all trials
                    
                    % Categorize and save the trial based on its type
                    if strcmp(current_trial.Type, 'Choice')
                        choice_trials{end+1} = current_trial;
                    elseif strcmp(current_trial.Type, 'Approach')
                        approach_trials{end+1} = current_trial;
                    elseif strcmp(current_trial.Type, 'Avoidance')
                        avoidance_trials{end+1} = current_trial;
                    end
                else
                    % Here, you can handle the case when the start event is missing
                    % For example, you can log this occurrence or simply ignore it
                    fprintf('Warning: End of trial detected without start event for index %d\n', i);
                end
                current_trial = []; % Reset for the next trial regardless
            end
        end        
        % Initialize variables for reward and punishment percentages
        for i = 1:length(approach_trials)
            approach_trials{i}.RewardPercentage = NaN;
            approach_trials{i}.PunishmentPercentage = NaN;
        end
        
        for i = 1:length(avoidance_trials)
            avoidance_trials{i}.RewardPercentage = NaN;
            avoidance_trials{i}.PunishmentPercentage = NaN;
        end
        
        % Process Approach and Avoidance trials to extract reward and punishment sizes and calculate percentages
        for trial_type = {'approach_trials', 'avoidance_trials'}
            trials = eval(trial_type{1});
            for i = 1:length(trials)
                % Find the index of the peripheral target eventmarker
                target_index = find(trials{i}.TTLValues == 251 | trials{i}.TTLValues == 252, 1, 'last');
                if ~isempty(target_index)
                    % Look for the next two values which are less than 200 and not 19 or 8
                    subsequent_values = trials{i}.TTLValues(target_index+1:end);
                    size_values = subsequent_values(subsequent_values < 200 & ~ismember(subsequent_values, [19, 8]));
                    if length(size_values) >= 2
                        trials{i}.RewardSize = size_values(1);
                        trials{i}.PunishmentSize = size_values(2);
                        % Calculate and store the percentages
                        trials{i}.RewardPercentage = (size_values(1) / 200) * 100;
                        trials{i}.PunishmentPercentage = (size_values(2) / 200) * 100;
                    else
                        warning(['Trial ' num2str(i) ' in ' trial_type{1} ' does not have two subsequent size values.']);
                    end
                else
                    warning(['No peripheral target eventmarker found in trial ' num2str(i) ' in ' trial_type{1} '.']);
                end
            end
            % Save back the processed trials
            eval([trial_type{1} ' = trials;']);
        end

        % Initialize the results matrix with the total number of trials
        total_trials = length(approach_trials) + length(avoidance_trials);
        results_matrix = zeros(total_trials, 3);
        
        % Initialize counters for approach and avoidance trials
        approach_index = 1;
        avoidance_index = 1;
        result_index = 1; % Keep track of the position in the results_matrix
        
        % Loop until we have placed all trials in the results_matrix
        while approach_index <= length(approach_trials) || avoidance_index <= length(avoidance_trials)
            
            if approach_index <= length(approach_trials) && ...
               (avoidance_index > length(avoidance_trials) || ...
                (approach_index <= length(approach_trials) && avoidance_index <= length(avoidance_trials) && ...
                 approach_trials{approach_index}.TimeStamps(1) <= avoidance_trials{avoidance_index}.TimeStamps(1)))
                
                % Add approach trial to the matrix and increment index
                results_matrix(result_index, :) = [approach_trials{approach_index}.RewardPercentage, ...
                                                   approach_trials{approach_index}.PunishmentPercentage, ...
                                                   1]; % 1 for approach
                approach_index = approach_index + 1;
            elseif avoidance_index <= length(avoidance_trials)
                % Add avoidance trial to the matrix and increment index
                results_matrix(result_index, :) = [avoidance_trials{avoidance_index}.RewardPercentage, ...
                                                   avoidance_trials{avoidance_index}.PunishmentPercentage, ...
                                                   0]; % 0 for avoidance
                avoidance_index = avoidance_index + 1;
            end
            
            result_index = result_index + 1; % Move to the next row in the results matrix
        end
      % 'avoidance_trials' with Avoidance trials within choice trials
        [b, ~, ~] = glmfit(results_matrix(:,1:2), results_matrix(:,3), 'binomial', 'link', 'logit');
        
        betaValues{1,x}=b;
        plot_x = 0:100;
        plot_y = (-b(1) - b(2).*plot_x)./b(3);
        plot_y1{1,x}=plot_y;
        
        % Extract all approach trials (where the third column is 1)
        approach_matrix = results_matrix(results_matrix(:,3) == 1, 1:2);
        approach_matrix1{1,x}=approach_matrix;
        
        % Extract all avoidance trials (where the third column is 0)
        avoidance_matrix = results_matrix(results_matrix(:,3) == 0, 1:2);
        avoidance_matrix1{1,x}=avoidance_matrix;
        
        % Define the event markers for start and end of reaction time for each type of trial
        start_marker = 234;
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
            start_index = find(event_markers == start_marker, 1, 'first');
            end_index_approach = find(event_markers == end_marker_approach, 1, 'first');
            end_index_avoidance = find(event_markers == end_marker_avoidance, 1, 'first');
        
            % Calculate reaction time for approach if the markers are present
            if ~isempty(start_index) && ~isempty(end_index_approach)
                rt_approach = timestamps(end_index_approach) - timestamps(start_index);
                approach_matrix_rt = [approach_matrix_rt; rt_approach];
                approach_matrix_rt1{1,x}=(approach_matrix_rt);
            end
            
            % Calculate reaction time for avoidance if the markers are present
            if ~isempty(start_index) && ~isempty(end_index_avoidance)
                rt_avoidance = timestamps(end_index_avoidance) - timestamps(start_index);
                avoidance_matrix_rt = [avoidance_matrix_rt; rt_avoidance];
                avoidance_matrix_rt1{1,x}=(avoidance_matrix_rt);
            end
        end

        if subplotIndex > 64
            currentFigure = currentFigure + 1;
            subplotIndex = 1;
            figure(currentFigure);
        end
        subplot(8, 8, subplotIndex);
        subplotIndex = subplotIndex + 1;
        
        sz=20;
        h(1)=scatter(approach_matrix(:,1),approach_matrix(:,2),sz,'+','g','LineWidth',1);
        hold on
        h(2)=scatter(avoidance_matrix(:,1),avoidance_matrix(:,2),sz,'s','r');
        hold on
        h(3)=plot(plot_x,plot_y,'-','LineWidth',2,'color','b');
        xlabel('Reward amount (%)','FontSize',10);
        ylabel('Airpuff amount (%)','FontSize',10);
        xlim([0 100]);
        ylim([0 100]);
        % lgd=legend(h,'Approach','Avoidance','Location','SouthEast');
        % lgd.FontSize = 8;
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 8)
        % legend('boxoff') 
        set(gcf,'color','w');
        title([num2str(formattedDate)]);
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
avoidance_matrix_rt1(1) = [];
approach_matrix_rt1(1) = [];
avoidance_matrix1(1)=[];
approach_matrix1(1)=[];

num_sessions=numel(avoidance_matrix1);
% Extend your grid boundaries to ensure coverage
x_edges = [-2.5:5:102.5];
y_edges = [-2.5:5:102.5];

% Initialize matrices to store counts for each session
num_sessions = length(avoidance_matrix1);
all_Approach_green_counts = zeros(length(y_edges)-1, length(x_edges)-1, num_sessions);
all_Avoidance_red_counts = zeros(length(y_edges)-1, length(x_edges)-1, num_sessions);

% Loop for counting
for i = 1:num_sessions
    for x_idx = 1:length(x_edges)-1
        for y_idx = 1:length(y_edges)-1
            x_center = (x_edges(x_idx) + x_edges(x_idx+1))/2;
            y_center = (y_edges(y_idx) + y_edges(y_idx+1))/2;

            all_Approach_green_counts(y_idx, x_idx, i) = sum(sqrt((approach_matrix1{1,i}(:,1)-x_center).^2 + (approach_matrix1{1,i}(:,2)-y_center).^2) <= 2.5);
            all_Avoidance_red_counts(y_idx, x_idx, i) = sum(sqrt((avoidance_matrix1{1,i}(:,1)-x_center).^2 + (avoidance_matrix1{1,i}(:,2)-y_center).^2) <= 2.5);
        end
    end
end

% Since we extended our grid, we'll remove the first and last rows/columns to get our original size back
all_Approach_green_counts = all_Approach_green_counts(2:end-1, 2:end-1, :);
all_Avoidance_red_counts = all_Avoidance_red_counts(2:end-1, 2:end-1, :);

% Average the count matrices across sessions
avg_Approach_green = mean(all_Approach_green_counts, 3);
avg_Avoidance_red = mean(all_Avoidance_red_counts, 3);

% Visualization
figure;
imagesc(x_edges(1:end-1), y_edges(1:end-1), avg_Approach_green - avg_Avoidance_red);  % Example: showing the difference between the two
colormap('jet');  % Or any other colormap of your choice
c = colorbar;
xlabel('Reward amount (%)');
ylabel('Airpuff amount (%)');
title('Average across all sessions');
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100);
set(c, 'YTickLabel', {'Av', '', '', '', '', '', '', '', '', '', 'Ap'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11));
box off;  % Removes the box around the figure

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  % Assuming b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 2, 'color', 'b');  % Average decision boundary
hold off;





%%%% for the reaction times
% Initialize matrices to store reaction time sums and counts for each trial type and each session
all_Approach_green_rt_sums = zeros(length(y_edges)-1, length(x_edges)-1, num_sessions);
all_Avoidance_red_rt_sums = zeros(length(y_edges)-1, length(x_edges)-1, num_sessions);
all_Approach_green_rt_counts = zeros(length(y_edges)-1, length(x_edges)-1, num_sessions);
all_Avoidance_red_rt_counts = zeros(length(y_edges)-1, length(x_edges)-1, num_sessions);

% Loop through each session and each bin to sum reaction times and count occurrences
for i = 1:num_sessions
    for x_idx = 1:length(x_edges)-1
        for y_idx = 1:length(y_edges)-1
            % Define the bin center
            x_center = (x_edges(x_idx) + x_edges(x_idx+1))/2;
            y_center = (y_edges(y_idx) + y_edges(y_idx+1))/2;

            % Find indices of approach and avoidance within this bin
            Approach_bin_indices = sqrt((approach_matrix1{1,i}(:,1)-x_center).^2 + (approach_matrix1{1,i}(:,2)-y_center).^2) <= 2.5;
            Avoidance_bin_indices = sqrt((avoidance_matrix1{1,i}(:,1)-x_center).^2 + (avoidance_matrix1{1,i}(:,2)-y_center).^2) <= 2.5;

            % Sum reaction times in the approach bins and count them
            all_Approach_green_rt_sums(y_idx, x_idx, i) = sum(approach_matrix_rt1{1,i}(Approach_bin_indices));
            all_Approach_green_rt_counts(y_idx, x_idx, i) = sum(Approach_bin_indices);

            % Sum reaction times in the avoidance bins and count them
            all_Avoidance_red_rt_sums(y_idx, x_idx, i) = sum(avoidance_matrix_rt1{1,i}(Avoidance_bin_indices));
            all_Avoidance_red_rt_counts(y_idx, x_idx, i) = sum(Avoidance_bin_indices);
        end
    end
end

% Calculate the average reaction times per bin, using nanmean to handle NaNs
avg_Approach_green_rt = nanmean(all_Approach_green_rt_sums ./ all_Approach_green_rt_counts, 3);
avg_Avoidance_red_rt = nanmean(all_Avoidance_red_rt_sums ./ all_Avoidance_red_rt_counts, 3);

% Trim the extended grid for the average reaction times
avg_Approach_green_rt = avg_Approach_green_rt(2:end-1, 2:end-1);
avg_Avoidance_red_rt = avg_Avoidance_red_rt(2:end-1, 2:end-1);

% Visualization with subplots
figure;

% Define the colormap to be used for both heatmaps
colormap('jet'); % A colormap that transitions from blue to red

% Subplot for Approach Reaction Times
subplot(1, 2, 1);
imagesc(x_edges(2:end-2), y_edges(2:end-2), avg_Approach_green_rt); % Use trimmed grid
caxis([300, 600]); % Set color scale from 0 to 500 milliseconds
colorbar; % Include a colorbar
xlabel('Reward amount (%)');
ylabel('Airpuff amount (%)');
title('Average Approach Reaction Time');
set(gca, 'YDir', 'normal'); % Correct the Y-axis direction

% Subplot for Avoidance Reaction Times
subplot(1, 2, 2);
imagesc(x_edges(2:end-2), y_edges(2:end-2), avg_Avoidance_red_rt); % Use trimmed grid
caxis([300, 600]); % Set color scale from 0 to 500 milliseconds
colorbar; % Include a colorbar
xlabel('Reward amount (%)');
ylabel('Airpuff amount (%)');
title('Average Avoidance Reaction Time');
set(gca, 'YDir', 'normal'); % Correct the Y-axis direction

% Adjust the subplots for better spacing
set(gcf, 'Position', [100, 100, 1200, 400]); % Set figure size and position

% Calculate the difference in average reaction times per bin between approach and avoidance
rt_difference = avg_Approach_green_rt - avg_Avoidance_red_rt;

% Visualization of the combined reaction time differences
figure;

% Heatmap for the difference in reaction times
imagesc(x_edges(1:end-1), y_edges(1:end-1), rt_difference); % Use trimmed grid
colormap('jet'); % A colormap that transitions from blue to red
c = colorbar; % Include a colorbar
caxis([-50 50]); % Set color scale range
% c.Label.String = 'RT Difference (Approach - Avoidance)';
xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Difference in Average Reaction Time (Approach - Avoidance)','FontSize', 24, 'FontName', 'Arial'); 
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100, 'FontSize', 20);
set(c, 'YTickLabel', {'Av', '', '', '', '', '', '', '', '', '', 'Ap'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11), 'FontSize', 20);
box off;  % Removes the box around the figure

% Adjust figure size and position
set(gcf, 'Position', [100, 100, 800, 600]); % Set figure size and position

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  % Assuming b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 5, 'color', 'k');  % Average decision boundary
hold off;


%%% for sum
rt_sum = (avg_Approach_green_rt + avg_Avoidance_red_rt)./2;

% Visualization of the combined reaction time differences
figure;

% Heatmap for the difference in reaction times
imagesc(x_edges(2:end-2), y_edges(2:end-2), rt_sum); % Use trimmed grid
colormap('jet'); % A colormap that transitions from blue to red
c = colorbar; % Include a colorbar
c.Label.String = 'RT Average (Approach + Avoidance)';
xlabel('Reward amount (%)');
ylabel('Airpuff amount (%)');
title('Sum of Average Reaction Time (Approach + Avoidance)');
set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100);
box off; % Removes the box around the figure

% Adjust figure size and position
set(gcf, 'Position', [100, 100, 800, 600]); % Set figure size and position

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  % Assuming b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 2, 'color', 'b');  % Average decision boundary
hold off;

% Visualization
figure;

% Define the Gaussian smoothing parameters
kernelSize = 5; % Size of the Gaussian kernel
sigma = 2; % Standard deviation of the Gaussian kernel
kernel = fspecial('gaussian', kernelSize, sigma); % Create the Gaussian kernel

% Apply Gaussian smoothing
smoothedData = imfilter(avg_Approach_green - avg_Avoidance_red, kernel, 'replicate');

% Display the smoothed data
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedData);  % Display the smoothed heatmap
colormap('jet');  % Or any other colormap of your choice
c = colorbar;

% Labeling and aesthetic settings
xlabel('Reward amount (%)', 'FontName', 'Arial', 'FontSize', 18); 
ylabel('Airpuff amount (%)', 'FontName', 'Arial', 'FontSize', 18); 
title('Average across all sessions', 'FontName', 'Arial', 'FontSize', 24); 
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100, 'FontName', 'Arial', 'FontSize', 18);
set(c, 'YTickLabel', {'Av', '', '', '', '', '', '', '', '', '', 'Ap'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11), 'FontName', 'Arial', 'FontSize', 18);
box off;  % Removes the box around the figure

% Decision boundary
hold on;
plot_x = linspace(0, 100, 100);  % More precise control over plotted range
avg_b = mean(cell2mat(betaValues), 2);  % Assuming b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 5, 'Color', 'k');  % Average decision boundary plotted in blue for visibility
hold off;

