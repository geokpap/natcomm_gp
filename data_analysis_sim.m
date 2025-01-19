close all;

% Clear potential variable conflict
clear strings;

load('sample_data.mat');
filtered_manipulation=matrix2;


% Filter the matrix to keep only the desired rows
matrix2=filtered_manipulation;

% Now filtered_manipulation contains only the rows that do not have 'MUA', 'poorunit' in the 5th column, 'Error' in the 62nd column, 
% and have either 'pACC' or 'cOFC' in the 6th column

% Initialize counters for each predictor in each brain structure
count_pACC = zeros(1, 8); % Count of predictors x1 to x8 in pACC
count_cOFC = zeros(1, 8); % Count of predictors x1 to x8 in cOFC
total_pACC = 0; % Total occurrences of pACC
total_cOFC = 0; % Total occurrences of cOFC

for row = 1:size(filtered_manipulation, 1)
    brain_structure = filtered_manipulation{row, 6};
    content = filtered_manipulation{row, 62};
    
    % Check if the content is a cell array before proceeding
    if ~iscell(content)
        disp(['Skipping row ' num2str(row) ' due to unexpected content: ' class(content)]);
        continue; % Skip this iteration and move to the next row
    end
    
    predictors = content; % Since we've confirmed it's a cell array, we can proceed
    
    % Update total counts for pACC or cOFC
    if strcmp(brain_structure, 'pACC')
        total_pACC = total_pACC + 1;
        for i = 1:length(predictors)
            predictor = predictors{i, 1};
            if ischar(predictor) || isstring(predictor)
                idx = str2double(extractAfter(predictor, 'x')); % Extract the number from 'x1', 'x2', etc.
                if ~isnan(idx) && idx > 0
                    count_pACC(idx) = count_pACC(idx) + 1;
                end
            end
        end
    elseif strcmp(brain_structure, 'cOFC')
        total_cOFC = total_cOFC + 1;
        for i = 1:length(predictors)
            predictor = predictors{i, 1};
            if ischar(predictor) || isstring(predictor)
                idx = str2double(extractAfter(predictor, 'x'));
                if ~isnan(idx) && idx > 0
                    count_cOFC(idx) = count_cOFC(idx) + 1;
                end
            end
        end
    end
end


% Calculate percentages
percentage_pACC = (count_pACC / total_pACC) * 100;
percentage_cOFC = (count_cOFC / total_cOFC) * 100;

% Display results
disp('Count and Percentage of Predictors in pACC:');
for i = 1:8
    fprintf('Predictor x%d: Count = %d, Percentage = %.2f%%\n', i, count_pACC(i), percentage_pACC(i));
end

disp('Count and Percentage of Predictors in cOFC:');
for i = 1:8
    fprintf('Predictor x%d: Count = %d, Percentage = %.2f%%\n', i, count_cOFC(i), percentage_cOFC(i));
end

% Create a matrix for the bar chart
barData = [percentage_pACC; percentage_cOFC]';

% Define the new labels for the predictors
labels = {'Rew', 'Ave', 'Eutil', 'Cho', 'Rew*Cho', 'Ave*Cho', 'Conf', 'RT'};
positions = 1:length(labels);

% Create a grouped bar chart
figure;
bar(barData);

% Customize the plotYou 
set(gca, 'XTick', positions, 'XTickLabel', labels, 'XTickLabelRotation', 45);
ylabel('Neurons (%)', 'FontSize', 20, 'FontWeight', 'bold');
title('Percentage of Predictors for Firing Rates in pACC and cOFC - Cue Period', 'FontSize', 22, 'FontWeight', 'bold');
legend({'pACC', 'cOFC'}, 'Location', 'bestoutside', 'Box', 'off', 'FontSize', 18);
set(gca, 'FontSize', 18, 'FontName', 'Arial');
set(gcf, 'Color', 'w'); % Set background color to white
box off;

% Ensure the plot is displayed appealingly
axis tight;
grid on;

% Change the plot to landscape orientation
set(gcf, 'PaperOrientation', 'landscape');

% Initialize variables to store the relevant data for Cue Period
pACC_data_cuePeriod = [];
cOFC_data_cuePeriod = [];

% Iterate through each row of filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Convert the value in the 67th column to a number
    stats_for_cuePeriod = filtered_manipulation{i, 67};

    % Check if the value is numeric, non-NaN, and less than 0.05
    if ~isnan(stats_for_cuePeriod) && stats_for_cuePeriod < 0.05
        % Check if the cell in the 63rd column is non-empty and numeric
        if ~isempty(filtered_manipulation{i, 63}) && isnumeric(filtered_manipulation{i, 63})
            dataValue = filtered_manipulation{i, 63};
            if strcmp(filtered_manipulation{i, 6}, 'pACC')
                pACC_data_cuePeriod = [pACC_data_cuePeriod; dataValue];
            elseif strcmp(filtered_manipulation{i, 6}, 'cOFC')
                cOFC_data_cuePeriod = [cOFC_data_cuePeriod; dataValue];
            end
        end
    end
end

% Calculate means if there is data
mean_pACC_cuePeriod = NaN;
mean_cOFC_cuePeriod = NaN;
if ~isempty(pACC_data_cuePeriod)
    mean_pACC_cuePeriod = mean(pACC_data_cuePeriod);
end
if ~isempty(cOFC_data_cuePeriod)
    mean_cOFC_cuePeriod = mean(cOFC_data_cuePeriod);
end


% Initialize variables to store the relevant data for Cue Period
pACC_data_AirpuffPeriod = [];
cOFC_data_AirpuffPeriod = [];

% Iterate through each row of filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Convert the value in the 68th column to a number
    stats_for_AirpuffPeriod = str2double(filtered_manipulation{i, 68});

    % Check if the value is numeric, non-NaN, and less than 0.05
    if ~isnan(stats_for_AirpuffPeriod) && stats_for_AirpuffPeriod < 0.05
        % Check if the cell in the 88rd column is non-empty and numeric
        if ~isempty(filtered_manipulation{i, 88}) && isnumeric(filtered_manipulation{i, 88})
            dataValue = filtered_manipulation{i, 88};
            if strcmp(filtered_manipulation{i, 6}, 'pACC')
                pACC_data_AirpuffPeriod = [pACC_data_AirpuffPeriod; dataValue];
            elseif strcmp(filtered_manipulation{i, 6}, 'cOFC')
                cOFC_data_AirpuffPeriod = [cOFC_data_AirpuffPeriod; dataValue];
            end
        end
    end
end

% Calculate means if there is data
mean_pACC_AirpuffPeriod = NaN;
mean_cOFC_AirpuffPeriod = NaN;
if ~isempty(pACC_data_AirpuffPeriod)
    mean_pACC_AirpuffPeriod = mean(pACC_data_AirpuffPeriod);
end
if ~isempty(cOFC_data_AirpuffPeriod)
    mean_cOFC_AirpuffPeriod = mean(cOFC_data_AirpuffPeriod);
end


%%%%%%% for small reward period
%%%%%%%%%%%
% Initialize variables to store the relevant data for Cue Period
pACC_data_SmallRewardPeriod = [];
cOFC_data_SmallRewardPeriod = [];

% Iterate through each row of filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Convert the value in the 69th column to a number
    stats_for_SmallRewardPeriod = str2double(filtered_manipulation{i, 69});

    % Check if the value is numeric, non-NaN, and less than 0.05
    if ~isnan(stats_for_SmallRewardPeriod) && stats_for_SmallRewardPeriod < 0.05
        % Check if the cell in the 89th column is non-empty and numeric
        if ~isempty(filtered_manipulation{i, 89}) && isnumeric(filtered_manipulation{i, 89})
            dataValue = filtered_manipulation{i, 89};
            if strcmp(filtered_manipulation{i, 6}, 'pACC')
                pACC_data_SmallRewardPeriod = [pACC_data_SmallRewardPeriod; dataValue];
            elseif strcmp(filtered_manipulation{i, 6}, 'cOFC')
                cOFC_data_SmallRewardPeriod = [cOFC_data_SmallRewardPeriod; dataValue];
            end
        end
    end
end

% Calculate means if there is data
mean_pACC_SmallRewardPeriod = NaN;
mean_cOFC_SmallRewardPeriod = NaN;
if ~isempty(pACC_data_SmallRewardPeriod)
    mean_pACC_SmallRewardPeriod = mean(pACC_data_SmallRewardPeriod);
end
if ~isempty(cOFC_data_SmallRewardPeriod)
    mean_cOFC_SmallRewardPeriod = mean(cOFC_data_SmallRewardPeriod);
end


%%%%%%%%%%%
%% for big reward
%%%%%%% for big reward period
%%%%%%%%%%%
% Initialize variables to store the relevant data for Cue Period
pACC_data_BigRewardPeriod = [];
cOFC_data_BigRewardPeriod = [];

% Iterate through each row of filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Convert the value in the 76th column to a number
    stats_for_BigRewardPeriod = str2double(filtered_manipulation{i, 76});

    % Check if the value is numeric, non-NaN, and less than 0.05
    if ~isnan(stats_for_BigRewardPeriod) && stats_for_BigRewardPeriod < 0.05
        % Check if the cell in the 90th column is non-empty and numeric
        if ~isempty(filtered_manipulation{i, 90}) && isnumeric(filtered_manipulation{i, 90})
            dataValue = filtered_manipulation{i, 90};
            if strcmp(filtered_manipulation{i, 6}, 'pACC')
                pACC_data_BigRewardPeriod = [pACC_data_BigRewardPeriod; dataValue];
            elseif strcmp(filtered_manipulation{i, 6}, 'cOFC')
                cOFC_data_BigRewardPeriod = [cOFC_data_BigRewardPeriod; dataValue];
            end
        end
    end
end

% Calculate means if there is data
mean_pACC_BigRewardPeriod = NaN;
mean_cOFC_BigRewardPeriod = NaN;
if ~isempty(pACC_data_BigRewardPeriod)
    mean_pACC_BigRewardPeriod = mean(pACC_data_BigRewardPeriod);
end
if ~isempty(cOFC_data_BigRewardPeriod)
    mean_cOFC_BigRewardPeriod = mean(cOFC_data_BigRewardPeriod);
end
%%%%%%%%%%

%%%%%%%%%%%
%% for red cue
%%%%%%% for Reward Cue Period
%%%%%%%%%%%
% Initialize variables to store the relevant data for Cue Period
pACC_data_RedCuePeriod = [];
cOFC_data_RedCuePeriod = [];

% Iterate through each row of filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Convert the value in the 72th column to a number
    stats_for_RedCuePeriod = str2double(filtered_manipulation{i, 72});

    % Check if the value is numeric, non-NaN, and less than 0.05
    if ~isnan(stats_for_RedCuePeriod) && stats_for_RedCuePeriod < 0.05
        % Check if the cell in the 91th column is non-empty and numeric
        if ~isempty(filtered_manipulation{i, 91}) && isnumeric(filtered_manipulation{i, 91})
            dataValue = filtered_manipulation{i, 91};
            if strcmp(filtered_manipulation{i, 6}, 'pACC')
                pACC_data_RedCuePeriod = [pACC_data_RedCuePeriod; dataValue];
            elseif strcmp(filtered_manipulation{i, 6}, 'cOFC')
                cOFC_data_RedCuePeriod = [cOFC_data_RedCuePeriod; dataValue];
            end
        end
    end
end

% Calculate means if there is data
mean_pACC_RedCuePeriod = NaN;
mean_cOFC_RedCuePeriod = NaN;
if ~isempty(pACC_data_RedCuePeriod)
    mean_pACC_RedCuePeriod = mean(pACC_data_RedCuePeriod);
end
if ~isempty(cOFC_data_RedCuePeriod)
    mean_cOFC_RedCuePeriod = mean(cOFC_data_RedCuePeriod);
end

%%%%%%%%%%%
%% for red cue
%%%%%%% for Reward Cue Period
% Initialize variables to store the relevant data for Cue Period
pACC_data_YellowCuePeriod = [];
cOFC_data_YellowCuePeriod = [];

% Iterate through each row of filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Convert the value in the 73th column to a number
    stats_for_YellowCuePeriod = str2double(filtered_manipulation{i, 73});

    % Check if the value is numeric, non-NaN, and less than 0.05
    if ~isnan(stats_for_YellowCuePeriod) && stats_for_YellowCuePeriod < 0.05
        % Check if the cell in the 92th column is non-empty and numeric
        if ~isempty(filtered_manipulation{i, 92}) && isnumeric(filtered_manipulation{i, 92})
            dataValue = filtered_manipulation{i, 92};
            if strcmp(filtered_manipulation{i, 6}, 'pACC')
                pACC_data_YellowCuePeriod = [pACC_data_YellowCuePeriod; dataValue];
            elseif strcmp(filtered_manipulation{i, 6}, 'cOFC')
                cOFC_data_YellowCuePeriod = [cOFC_data_YellowCuePeriod; dataValue];
            end
        end
    end
end

% Calculate means if there is data
mean_pACC_YellowCuePeriod = NaN;
mean_cOFC_YellowCuePeriod = NaN;
if ~isempty(pACC_data_YellowCuePeriod)
    mean_pACC_YellowCuePeriod = mean(pACC_data_YellowCuePeriod);
end
if ~isempty(cOFC_data_YellowCuePeriod)
    mean_cOFC_YellowCuePeriod = mean(cOFC_data_YellowCuePeriod);
end
%%%%%%%%%%
%%%%%%%%%%%
%% for red reward
%%%%%%% for Reward Cue Period
%%%%%%%%%%%
% Initialize variables to store the relevant data for Cue Period
pACC_data_RedRewardPeriod = [];
cOFC_data_RedRewardPeriod = [];

% Iterate through each row of filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Convert the value in the 74th column to a number
    stats_for_RedRewardPeriod = str2double(filtered_manipulation{i, 74});

    % Check if the value is numeric, non-NaN, and less than 0.05
    if ~isnan(stats_for_RedRewardPeriod) && stats_for_RedRewardPeriod < 0.05
        % Check if the cell in the 93th column is non-empty and numeric
        if ~isempty(filtered_manipulation{i, 93}) && isnumeric(filtered_manipulation{i, 93})
            dataValue = filtered_manipulation{i, 93};
            if strcmp(filtered_manipulation{i, 6}, 'pACC')
                pACC_data_RedRewardPeriod = [pACC_data_RedRewardPeriod; dataValue];
            elseif strcmp(filtered_manipulation{i, 6}, 'cOFC')
                cOFC_data_RedRewardPeriod = [cOFC_data_RedRewardPeriod; dataValue];
            end
        end
    end
end

% Calculate means if there is data
mean_pACC_RedRewardPeriod = NaN;
mean_cOFC_RedRewardPeriod = NaN;
if ~isempty(pACC_data_RedRewardPeriod)
    mean_pACC_RedRewardPeriod = mean(pACC_data_RedRewardPeriod);
end
if ~isempty(cOFC_data_RedRewardPeriod)
    mean_cOFC_RedRewardPeriod = mean(cOFC_data_RedRewardPeriod);
end
%%%%%%%%%%
%%%%%%%%%%%
%% for red reward
%%%%%%% for Reward Cue Period
%%%%%%%%%%%
% Initialize variables to store the relevant data for Cue Period
pACC_data_YellowAirpuffPeriod = [];
cOFC_data_YellowAirpuffPeriod = [];

% Iterate through each row of filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Convert the value in the 75th column to a number
    stats_for_YellowAirpuffPeriod = str2double(filtered_manipulation{i, 75});

    % Check if the value is numeric, non-NaN, and less than 0.05
    if ~isnan(stats_for_YellowAirpuffPeriod) && stats_for_YellowAirpuffPeriod < 0.05
        % Check if the cell in the 94th column is non-empty and numeric
        if ~isempty(filtered_manipulation{i, 94}) && isnumeric(filtered_manipulation{i, 94})
            dataValue = filtered_manipulation{i, 94};
            if strcmp(filtered_manipulation{i, 6}, 'pACC')
                pACC_data_YellowAirpuffPeriod = [pACC_data_YellowAirpuffPeriod; dataValue];
            elseif strcmp(filtered_manipulation{i, 6}, 'cOFC')
                cOFC_data_YellowAirpuffPeriod = [cOFC_data_YellowAirpuffPeriod; dataValue];
            end
        end
    end
end

% Calculate means if there is data
mean_pACC_YellowAirpuffPeriod = NaN;
mean_cOFC_YellowAirpuffPeriod = NaN;
if ~isempty(pACC_data_YellowAirpuffPeriod)
    mean_pACC_YellowAirpuffPeriod = mean(pACC_data_YellowAirpuffPeriod);
end
if ~isempty(cOFC_data_YellowAirpuffPeriod)
    mean_cOFC_YellowAirpuffPeriod = mean(cOFC_data_YellowAirpuffPeriod);
end
%%%%%%%%%%

% Plotting
figure;
subplot(2, 1, 1); % First subplot in a 2x1 grid
bar(mean_pACC_cuePeriod, 'FaceColor', [0.8500 0.3250 0.0980]); % A distinct color for the bar
title('Mean Firing Rates for pACC (<0.05 in 39th Column)');
xlabel('Bin');
ylabel('Mean Firing Rate');

subplot(2, 1, 2); % Second subplot in a 2x1 grid
bar(mean_cOFC_cuePeriod, 'FaceColor', [0 0.4470 0.7410]); % A distinct color for the bar
title('Mean Firing Rates for cOFC (<0.05 in 39th Column)');
xlabel('Bin');
ylabel('Mean Firing Rate');

% Plotting in a single figure with overlaid bars
figure;

% Bar plot for mean_pACC_cuePeriod
bar(mean_pACC_cuePeriod, 'FaceColor', [0.8500 0.3250 0.0980]); % Distinct color for pACC
hold on; % Hold the current plot

% Bar plot for mean_cOFC_cuePeriod overlaid
bar(mean_cOFC_cuePeriod, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.5); % Distinct color for cOFC and make it slightly transparent

% Add a vertical line at index 40
xline(40, 'LineWidth', 2, 'Color', 'k'); % Black line with a width of 2

hold off; % Release the hold

% Enhancements for better visibility
title('Overlayed Mean Firing Rates for pACC and cOFC (p < 0.05 at Cue/Decision Period)','FontSize',30);
xlabel('Time (bin size = 50ms)','FontSize',24);
ylabel('Mean Firing Rate (Hz)','FontSize',24);
legend('pACC', 'cOFC');
set(gca, 'FontSize', 22); % Set font size for axis tick labels
box off
% Make sure all parts of the figure are updated
set(findall(gcf,'-property','FontSize'),'FontSize',18);

% Plotting in a single figure with overlaid bars
figure;

% Bar plot for mean_pACC_cuePeriod
bar(mean_pACC_cuePeriod, 'FaceColor', [0.8500 0.3250 0.0980]); % Distinct color for pACC
hold on; % Hold the current plot

% Bar plot for mean_cOFC_cuePeriod overlaid
bar(mean_cOFC_cuePeriod, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.5); % Distinct color for cOFC and make it slightly transparent

% Add a vertical line at index 40
xline(40, 'LineWidth', 2, 'Color', 'k'); % Black line with a width of 2

hold off; % Release the hold

% Enhancements for better visibility
title('Overlayed Mean Firing Rates for pACC and cOFC (p < 0.05 at Cue/Decision Period)', 'FontSize', 40);
xlabel('Time (bin size = 50ms)', 'FontSize', 32);
ylabel('Mean Firing Rate (Hz)', 'FontSize', 32);

% Set legend with increased font size
legend('pACC', 'cOFC', 'FontSize', 24);

% Set font size for axis tick labels and other figure elements
set(gca, 'FontSize', 28); % Set font size for axis tick labels

% Turn off the box around the plot to make it less cluttered
box off;

% Ensure consistent font size throughout all figure elements
set(findall(gcf,'-property','FontSize'),'FontSize',26); % This might be redundant after explicitly setting other font sizes



% Initialize counters for "cOFC" and "pACC" for both cue and airpuff periods and small reward periods 
count_cOFC_1_cue = 0;
count_cOFC_minus1_cue = 0;
count_cOFC_1_airpuff = 0;
count_cOFC_minus1_airpuff = 0;
count_cOFC_1_smallreward = 0;
count_cOFC_minus1_smallreward = 0;
count_cOFC_1_bigreward = 0;
count_cOFC_minus1_bigreward = 0;
count_cOFC_1_redcue = 0;
count_cOFC_minus1_redcue = 0;
count_cOFC_1_yellowcue = 0;
count_cOFC_minus1_yellowcue = 0;
count_cOFC_1_redreward = 0;
count_cOFC_minus1_redreward = 0;
count_cOFC_1_yellowairpuff = 0;
count_cOFC_minus1_yellowairpuff = 0;

count_pACC_1_cue = 0;
count_pACC_minus1_cue = 0;
count_pACC_1_airpuff = 0;
count_pACC_minus1_airpuff = 0;
count_pACC_1_smallreward = 0;
count_pACC_minus1_smallreward = 0;
count_pACC_1_bigreward = 0;
count_pACC_minus1_bigreward = 0;
count_pACC_1_redcue = 0;
count_pACC_minus1_redcue = 0;
count_pACC_1_yellowcue = 0;
count_pACC_minus1_yellowcue = 0;
count_pACC_1_redreward = 0;
count_pACC_minus1_redreward = 0;
count_pACC_1_yellowairpuff = 0;
count_pACC_minus1_yellowairpuff = 0;

% Number of rows in filtered_manipulation
nRows = size(filtered_manipulation, 1);

% Iterate through each row
for i = 1:nRows
    % Check for "cOFC"
    if strcmp(filtered_manipulation{i, 6}, 'cOFC')
        % Cue period (column 67)
        if filtered_manipulation{i, 67} < 0.05
            if filtered_manipulation{i, 64} == 1
                count_cOFC_1_cue = count_cOFC_1_cue + 1;
            elseif filtered_manipulation{i, 64} == -1
                count_cOFC_minus1_cue = count_cOFC_minus1_cue + 1;
            end
        end
        % Airpuff period (column 68)
        if filtered_manipulation{i, 68} < 0.05
            if filtered_manipulation{i, 65} == 1
                count_cOFC_1_airpuff = count_cOFC_1_airpuff + 1;
            elseif filtered_manipulation{i, 65} == -1
                count_cOFC_minus1_airpuff = count_cOFC_minus1_airpuff + 1;
            end
        end
    % Small reward period (column 69)
        if filtered_manipulation{i, 69} < 0.05
            if filtered_manipulation{i, 80} == 1
                count_cOFC_1_smallreward = count_cOFC_1_smallreward + 1;
            elseif filtered_manipulation{i, 80} == -1
                count_cOFC_minus1_smallreward = count_cOFC_minus1_smallreward + 1;
            end
        end

    % big reward period (column 76)
        if filtered_manipulation{i, 76} < 0.05
            if filtered_manipulation{i, 87} == 1
                count_cOFC_1_bigreward = count_cOFC_1_bigreward + 1;
            elseif filtered_manipulation{i, 87} == -1
                count_cOFC_minus1_bigreward = count_cOFC_minus1_bigreward + 1;
            end
        end

    % Reward Cue Period (column 72)
        if filtered_manipulation{i, 72} < 0.05
            if filtered_manipulation{i, 83} == 1
                count_cOFC_1_redcue = count_cOFC_1_redcue + 1;
            elseif filtered_manipulation{i, 83} == -1
                count_cOFC_minus1_redcue = count_cOFC_minus1_redcue + 1;
            end
        end

            % yellow cue period (column 73)
        if filtered_manipulation{i, 73} < 0.05
            if filtered_manipulation{i, 84} == 1
                count_cOFC_1_yellowcue = count_cOFC_1_yellowcue + 1;
            elseif filtered_manipulation{i, 84} == -1
                count_cOFC_minus1_yellowcue = count_cOFC_minus1_yellowcue + 1;
            end
        end

            % red reward period (column 74)
        if filtered_manipulation{i, 74} < 0.05
            if filtered_manipulation{i, 85} == 1
                count_cOFC_1_redreward = count_cOFC_1_redreward + 1;
            elseif filtered_manipulation{i, 85} == -1
                count_cOFC_minus1_redreward = count_cOFC_minus1_redreward + 1;
            end
        end

            % big reward period (column 75)
        if filtered_manipulation{i, 75} < 0.05
            if filtered_manipulation{i, 86} == 1
                count_cOFC_1_yellowairpuff = count_cOFC_1_yellowairpuff + 1;
            elseif filtered_manipulation{i, 86} == -1
                count_cOFC_minus1_yellowairpuff = count_cOFC_minus1_yellowairpuff + 1;
            end
        end

    % Check for "pACC"
    elseif strcmp(filtered_manipulation{i, 6}, 'pACC')
        % Cue period (column 67)
        if filtered_manipulation{i, 67} < 0.05
            if filtered_manipulation{i, 64} == 1
                count_pACC_1_cue = count_pACC_1_cue + 1;
            elseif filtered_manipulation{i, 64} == -1
                count_pACC_minus1_cue = count_pACC_minus1_cue + 1;
            end
        end
        % Outcome period (column 68)
        if filtered_manipulation{i, 68} < 0.05
            if filtered_manipulation{i, 65} == 1
                count_pACC_1_airpuff = count_pACC_1_airpuff + 1;
            elseif filtered_manipulation{i, 65} == -1
                count_pACC_minus1_airpuff = count_pACC_minus1_airpuff + 1;
            end
        end

           % Small reward period (column 69)
        if filtered_manipulation{i, 69} < 0.05
            if filtered_manipulation{i, 80} == 1
                count_pACC_1_smallreward = count_pACC_1_smallreward + 1;
            elseif filtered_manipulation{i, 80} == -1
                count_pACC_minus1_smallreward = count_pACC_minus1_smallreward + 1;
            end
        end

           % big reward period (column 76)
        if filtered_manipulation{i, 72} < 0.05
            if filtered_manipulation{i, 83} == 1
                count_pACC_1_redcue = count_pACC_1_redcue + 1;
            elseif filtered_manipulation{i, 83} == -1
                count_pACC_minus1_redcue = count_pACC_minus1_redcue + 1;
            end
        end

           % big reward period (column 76)
        if filtered_manipulation{i, 73} < 0.05
            if filtered_manipulation{i, 84} == 1
                count_pACC_1_yellowcue = count_pACC_1_yellowcue + 1;
            elseif filtered_manipulation{i, 84} == -1
                count_pACC_minus1_yellowcue = count_pACC_minus1_yellowcue + 1;
            end
        end

                   % big reward period (column 76)
        if filtered_manipulation{i, 74} < 0.05
            if filtered_manipulation{i, 85} == 1
                count_pACC_1_redreward = count_pACC_1_redreward + 1;
            elseif filtered_manipulation{i, 85} == -1
                count_pACC_minus1_redreward = count_pACC_minus1_redreward + 1;
            end
        end

                   % big reward period (column 76)
        if filtered_manipulation{i, 75} < 0.05
            if filtered_manipulation{i, 86} == 1
                count_pACC_1_yellowairpuff = count_pACC_1_yellowairpuff + 1;
            elseif filtered_manipulation{i, 86} == -1
                count_pACC_minus1_yellowairpuff = count_pACC_minus1_yellowairpuff + 1;
            end
        end

                   % big reward period (column 76)
        if filtered_manipulation{i, 76} < 0.05
            if filtered_manipulation{i, 87} == 1
                count_pACC_1_bigreward = count_pACC_1_bigreward + 1;
            elseif filtered_manipulation{i, 87} == -1
                count_pACC_minus1_bigreward = count_pACC_minus1_bigreward + 1;
            end
        end


    end
end

% Calculate the ratios for cue and airpuff periods
ratio_cOFC_cue = calculateRatio(count_cOFC_1_cue, count_cOFC_minus1_cue);
ratio_cOFC_airpuff = calculateRatio(count_cOFC_1_airpuff, count_cOFC_minus1_airpuff);
ratio_cOFC_smallreward = calculateRatio(count_cOFC_1_smallreward, count_cOFC_minus1_smallreward);
ratio_cOFC_bigreward = calculateRatio(count_cOFC_1_bigreward, count_cOFC_minus1_bigreward);
ratio_cOFC_redcue = calculateRatio(count_cOFC_1_redcue, count_cOFC_minus1_redcue);
ratio_cOFC_yellowcue = calculateRatio(count_cOFC_1_yellowcue, count_cOFC_minus1_yellowcue);
ratio_cOFC_redreward = calculateRatio(count_cOFC_1_redreward, count_cOFC_minus1_redreward);
ratio_cOFC_yellowairpuff = calculateRatio(count_cOFC_1_yellowairpuff, count_cOFC_minus1_yellowairpuff);


ratio_pACC_cue = calculateRatio(count_pACC_1_cue, count_pACC_minus1_cue);
ratio_pACC_airpuff = calculateRatio(count_pACC_1_airpuff, count_pACC_minus1_airpuff);
ratio_pACC_smallreward = calculateRatio(count_pACC_1_smallreward, count_pACC_minus1_smallreward);
ratio_pACC_bigreward = calculateRatio(count_pACC_1_bigreward, count_pACC_minus1_bigreward);
ratio_pACC_redcue = calculateRatio(count_pACC_1_redcue, count_pACC_minus1_redcue);
ratio_pACC_yellowcue = calculateRatio(count_pACC_1_yellowcue, count_pACC_minus1_yellowcue);
ratio_pACC_redreward = calculateRatio(count_pACC_1_redreward, count_pACC_minus1_redreward);
ratio_pACC_yellowairpuff = calculateRatio(count_pACC_1_yellowairpuff, count_pACC_minus1_yellowairpuff);

% Display the results
disp('Cue Period:');
disp(['cOFC (1/-1) with p < 0.05: ', num2str(ratio_cOFC_cue)]);
disp(['pACC (1/-1) with p < 0.05: ', num2str(ratio_pACC_cue)]);
disp('Airpuff Period:');
disp(['cOFC (1/-1) with p < 0.05: ', num2str(ratio_cOFC_airpuff)]);
disp(['pACC (1/-1) with p < 0.05: ', num2str(ratio_pACC_airpuff)]);
disp('Small reward Period:');
disp(['cOFC (1/-1) with p < 0.05: ', num2str(ratio_cOFC_smallreward)]);
disp(['pACC (1/-1) with p < 0.05: ', num2str(ratio_pACC_smallreward)]);
disp('Big reward Period:');
disp(['cOFC (1/-1) with p < 0.05: ', num2str(ratio_cOFC_bigreward)]);
disp(['pACC (1/-1) with p < 0.05: ', num2str(ratio_pACC_bigreward)]);
disp('Reward Cue Period:');
disp(['cOFC (1/-1) with p < 0.05: ', num2str(ratio_cOFC_redcue)]);
disp(['pACC (1/-1) with p < 0.05: ', num2str(ratio_pACC_redcue)]);
disp('Yellow cue Period:');
disp(['cOFC (1/-1) with p < 0.05: ', num2str(ratio_cOFC_yellowcue)]);
disp(['pACC (1/-1) with p < 0.05: ', num2str(ratio_pACC_yellowcue)]);
disp('Red Reward Period:');
disp(['cOFC (1/-1) with p < 0.05: ', num2str(ratio_cOFC_redreward)]);
disp(['pACC (1/-1) with p < 0.05: ', num2str(ratio_pACC_redreward)]);
disp('Yellow airpuff Period:');
disp(['cOFC (1/-1) with p < 0.05: ', num2str(ratio_cOFC_yellowairpuff)]);
disp(['pACC (1/-1) with p < 0.05: ', num2str(ratio_pACC_yellowairpuff)]);

% Event labels with specified names
events = {'Fix point ON', 'Choice cue ON', 'Ap air ON', 'Av rew ON', 'Ap highlight', ...
          'Av highlight', 'Red bar ON', 'Yellow bar ON', 'Red bar REW', 'Yellow air ON', 'Ap rew ON'};

%%%%%% cleaning up
% Initialize a logical vector to mark rows to keep
keepRows = false(size(filtered_manipulation, 1), 1);

% Iterate through each row to check the condition
for i = 1:size(filtered_manipulation, 1)
    % Extract the relevant columns for the current row
    rowData = cell2mat(filtered_manipulation(i, 66:76));
    
    % Check if the row is not all NaNs or all zeros
    if ~all(isnan(rowData) | rowData == 0)
        keepRows(i) = true; % Mark row to be kept
    end
end

% Filter the matrix to keep only the desired rows
filtered_manipulation= filtered_manipulation(keepRows, :);

%%%%%%%% for the event related plots
% Initialize counters
total_cOFC = 0;
total_pACC = 0;
significant_cOFC = zeros(1, 11); % For 10 events (columns 66 to 75)
significant_pACC = zeros(1, 11);

% Loop through each row
for i = 1:size(filtered_manipulation, 1)
    unitType = filtered_manipulation{i, 6};
    
    % Initialize temporary counters for significant events per row
    tempSignificant_cOFC = zeros(1, 11);
    tempSignificant_pACC = zeros(1, 11);
    
    % Check unit type
    if strcmp(unitType, 'cOFC')
        total_cOFC = total_cOFC + 1; % Increment total counter for cOFC
    elseif strcmp(unitType, 'pACC')
        total_pACC = total_pACC + 1; % Increment total counter for pACC
    end
    
    % Iterate through columns 66 to 75 to check significance
    for col = 66:76
        pValue = filtered_manipulation{i, col};
        if pValue < 0.05
            if strcmp(unitType, 'cOFC')
                tempSignificant_cOFC(col-65) = 1; % Mark as significant for cOFC
            elseif strcmp(unitType, 'pACC')
                tempSignificant_pACC(col-65) = 1; % Mark as significant for pACC
            end
        end
    end
    
    % Update the significant counters
    significant_cOFC = significant_cOFC + tempSignificant_cOFC;
    significant_pACC = significant_pACC + tempSignificant_pACC;
end

% Calculate percentages
percent_significant_cOFC = (significant_cOFC / total_cOFC) * 100;
percent_significant_pACC = (significant_pACC / total_pACC) * 100;

% Display the results
disp('Percentage of significant cOFC units across events (66 to 75):');
disp(percent_significant_cOFC);
disp('Percentage of significant pACC units across events (66 to 75):');
disp(percent_significant_pACC);

% Event labels with specified names
events = {'Fix point ON', 'Choice cue ON', 'Ap air ON', 'Ap rew ON', 'Av rew ON', 'Ap highlight', ...
          'Av highlight', 'Red bar ON', 'Yellow bar ON', 'Red bar REW', 'Yellow air ON'};

% Assuming cOFC and pACC contain the relevant data for plotting
% Placeholder data for cOFC and pACC; replace these with the actual data
cOFC_percent = percent_significant_cOFC; % the actual cOFC data here
pACC_percent = percent_significant_pACC; % the actual pACC data here

% Set up the figure size for a landscape style
figure('Units', 'inches', 'Position', [0.5, 0.5, 14, 8]); % 14x8 inches figure

% Calculate positions and width for the bars to accommodate two groups
numEvents = length(events);
positions = 1:numEvents;
width = 0.3; % Width for two groups

% Define subtle colors for cOFC and pACC
colors = [0.7, 0.7, 0.9; 0.8, 0.8, 0.7]; % Subtle colors for two datasets

hold on;

% Create the bar plots for cOFC and pACC
bar(positions - width/2, [cOFC_percent(1:3),cOFC_percent(11),cOFC_percent(4:10)], width, 'FaceColor', colors(1,:));
bar(positions + width/2, [pACC_percent(1:3),pACC_percent(11),pACC_percent(4:10)], width, 'FaceColor', colors(2,:));

hold off;

% Customize the plot
set(gca, 'XTick', positions, 'XTickLabel', events, 'XTickLabelRotation', 45);
ylabel('Neurons (%)', 'FontSize', 20, 'FontWeight', 'bold');
title('Neural Activity by Event and Brain Region during the Ap-Av task', 'FontSize', 22, 'FontWeight', 'bold');
legend({'cOFC', 'pACC'}, 'Location', 'bestoutside', 'Box', 'off', 'FontSize', 18);
set(gca, 'FontSize', 18, 'FontName', 'Arial');
set(gcf, 'Color', 'w'); % Set background color to white
box off;

% Ensure the plot is displayed appealingly
axis tight;
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assuming cOFC and pACC contain the relevant data for plotting
% Placeholder data for cOFC and pACC; replace these with the actual data
cOFC = significant_cOFC; % the actual cOFC data here
pACC = significant_pACC; % the actual pACC data here

% Set up the figure size for a landscape style
figure('Units', 'inches', 'Position', [0.5, 0.5, 14, 8]); % 14x8 inches figure

% Calculate positions and width for the bars to accommodate two groups
numEvents = length(events);
positions = 1:numEvents;
width = 0.3; % Width for two groups

% Define subtle colors for cOFC and pACC
colors = [0.7, 0.7, 0.9; 0.8, 0.8, 0.7]; % Subtle colors for two datasets

hold on;

% Create the bar plots for cOFC and pACC
bar(positions - width/2, cOFC, width, 'FaceColor', colors(1,:));
bar(positions + width/2, pACC, width, 'FaceColor', colors(2,:));

hold off;

% Customize the plot
set(gca, 'XTick', positions, 'XTickLabel', events, 'XTickLabelRotation', 45);
ylabel('Neurons (%)', 'FontSize', 20, 'FontWeight', 'bold');
title('Neural Activity by Event and Brain Region during the Ap-Av task', 'FontSize', 22, 'FontWeight', 'bold');
legend({'cOFC', 'pACC'}, 'Location', 'bestoutside', 'Box', 'off', 'FontSize', 18);
set(gca, 'FontSize', 18, 'FontName', 'Arial');
set(gcf, 'Color', 'w'); % Set background color to white
box off;

% Ensure the plot is displayed appealingly
axis tight;
grid on;

%%%%%PLOT THE RATIO PLOTS
% Number of events
numEvents = 10;

% Arrays to store the counts
excitationCounts_cOFC = zeros(1, numEvents);
suppressionCounts_cOFC = zeros(1, numEvents);
totalSignificant_cOFC = zeros(1, numEvents);

excitationCounts_pACC = zeros(1, numEvents);
suppressionCounts_pACC = zeros(1, numEvents);
totalSignificant_pACC = zeros(1, numEvents);

% Loop through each row in filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    unitType = filtered_manipulation{i, 6}; % cOFC or pACC
    
    % Loop through events (columns 66 to 76 for significance, 76 to 85 for ratios)
    for j = 1:numEvents
        significance = filtered_manipulation{i, 65 + j}; % Significance
        ratio = filtered_manipulation{i, 76 + j}; % Ratio
        
        % Check if the event is significant
        if significance < 0.05
            if ratio == 1 || ratio == -1 % Count only significant excitation or suppression
                if strcmp(unitType, 'cOFC')
                    totalSignificant_cOFC(j) = totalSignificant_cOFC(j) + 1;
                    if ratio == 1
                        excitationCounts_cOFC(j) = excitationCounts_cOFC(j) + 1;
                    else
                        suppressionCounts_cOFC(j) = suppressionCounts_cOFC(j) + 1;
                    end
                elseif strcmp(unitType, 'pACC')
                    totalSignificant_pACC(j) = totalSignificant_pACC(j) + 1;
                    if ratio == 1
                        excitationCounts_pACC(j) = excitationCounts_pACC(j) + 1;
                    else
                        suppressionCounts_pACC(j) = suppressionCounts_pACC(j) + 1;
                    end
                end
            end
        end
    end
end

% Calculate the percentages
percentExcitation_cOFC = (excitationCounts_cOFC ./ totalSignificant_cOFC) * 100;
percentSuppression_cOFC = (suppressionCounts_cOFC ./ totalSignificant_cOFC) * 100;
percentExcitation_pACC = (excitationCounts_pACC ./ totalSignificant_pACC) * 100;
percentSuppression_pACC = (suppressionCounts_pACC ./ totalSignificant_pACC) * 100;

% Handling divisions by zero
percentExcitation_cOFC(isnan(percentExcitation_cOFC)) = 0;
percentSuppression_cOFC(isnan(percentSuppression_cOFC)) = 0;
percentExcitation_pACC(isnan(percentExcitation_pACC)) = 0;
percentSuppression_pACC(isnan(percentSuppression_pACC)) = 0;

% Now, you can plot these percentages for cOFC and pACC for each event

% Plotting
events = {'Fix point ON', 'Choice cue ON', 'Ap air ON', 'Av rew ON', 'Ap highlight', ...
          'Av highlight', 'Red bar ON', 'Yellow bar ON', 'Red bar REW', 'Yellow air ON', 'Ap rew ON'};

bar(mean_pACC_cuePeriod, 'FaceColor', [0.8500 0.3250 0.0980]); % Distinct color for pACC
hold on; % Hold the current plot

% Bar plot for mean_cOFC_cuePeriod overlaid
bar(mean_cOFC_cuePeriod, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.5); % Distinct color for cOFC and make it slightly transparent

% Add a vertical line at index 40
xline(40, 'LineWidth', 2, 'Color', 'k'); % Black line with a width of 2

% Plot for Excitation
figure;
subplot(2, 1, 1); % Excitation plot
bar(1:10, [percentExcitation_cOFC; percentExcitation_pACC]', 'grouped');
title('Excitation: Significant Units','FontSize',26);
legend({'cOFC', 'pACC'}, 'Location', 'bestoutside');
set(gca, 'XTick', 1:10, 'XTickLabel', events, 'XTickLabelRotation', 45,'FontSize',18);
ylabel('Percentage (%)','FontSize',24);
xlabel('Events','FontSize',24);
box off;

% Plot for Suppression
subplot(2, 1, 2); % Suppression plot
bar(1:10, [percentSuppression_cOFC; percentSuppression_pACC]', 'grouped');
title('Suppression: Significant Units','FontSize',24);
legend({'cOFC', 'pACC'}, 'Location', 'bestoutside');
set(gca, 'XTick', 1:10, 'XTickLabel', events, 'XTickLabelRotation', 45,'FontSize',18);
ylabel('Percentage (%)','FontSize',24);
xlabel('Events','FontSize',24);
box off;

% Adjust figure properties
set(gcf, 'Color', 'w'); % Set background color to white

% Plot for Excitation
figure;
subplot(2, 1, 1); % Excitation plot
hb = bar(1:10, [percentExcitation_cOFC; percentExcitation_pACC]', 'grouped');
hb(1).FaceColor = [0 0.4470 0.7410]; % Set cOFC color
hb(1).FaceAlpha = 0.5;               % Set cOFC transparency
hb(2).FaceColor = [0.8500 0.3250 0.0980]; % Set pACC color
title('Excitation: Significant Units','FontSize',28);
legend([hb(2), hb(1)], {'pACC', 'cOFC'}, 'Location', 'bestoutside'); % Switched order in legend
set(gca, 'XTick', 1:10, 'XTickLabel', events, 'XTickLabelRotation', 45,'FontSize',32);
ylabel('Percentage (%)','FontSize',28);
xlabel('Events','FontSize',28);
box off;

% Plot for Suppression
subplot(2, 1, 2); % Suppression plot
hb2 = bar(1:10, [percentSuppression_cOFC; percentSuppression_pACC]', 'grouped');
hb2(1).FaceColor = [0 0.4470 0.7410]; % Set cOFC color
hb2(1).FaceAlpha = 0.5;               % Set cOFC transparency
hb2(2).FaceColor = [0.8500 0.3250 0.0980]; % Set pACC color
title('Inhibition: Significant Units','FontSize',28);
legend([hb2(2), hb2(1)], {'pACC', 'cOFC'}, 'Location', 'bestoutside'); % Switched order in legend
set(gca, 'XTick', 1:10, 'XTickLabel', events, 'XTickLabelRotation', 45,'FontSize',32);
ylabel('Percentage (%)','FontSize',28);
xlabel('Events','FontSize',28);
box off;

% Adjust figure properties
set(gcf, 'Color', 'w'); % Set background color to white

% Assuming matrix2 is the 2897x85 cell matrix
% Initialize an empty cell array for cOFC data
cOFC_matrix2 = {};

% Logical indexing to find rows where the 6th column is 'cOFC'
cOFC_rows = strcmp(matrix2(:, 6), 'cOFC');

% Use the logical index to extract the corresponding rows from matrix2
cOFC_matrix2 = matrix2(cOFC_rows, :);

% Filter rows based on the value in column 64 being 1 or -1
subset1 = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 64)) == 1, :); % Rows with value 1 in column 64
subset2 = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 64)) == -1, :); % Rows with value -1 in column 64

% Filter rows based on the value in column 65 being 1 or -1
subset1_air = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 65)) == 1, :); % Rows with value 1 in column 65
subset2_air = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 65)) == -1, :); % Rows with value -1 in column 65

% Filter rows based on the value in column 65 being 1 or -1
subset1_smallreward = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 80)) == 1, :); % Rows with value 1 in column 65
subset2_smallreward = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 80)) == -1, :); % Rows with value -1 in column 65

% Filter rows based on the value in column 65 being 1 or -1
subset1_bigreward = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 87)) == 1, :); % Rows with value 1 in column 65
subset2_bigreward = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 87)) == -1, :); % Rows with value -1 in column 65

% Filter rows based on the value in column 65 being 1 or -1
subset1_redcue = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 83)) == 1, :); % Rows with value 1 in column 65
subset2_redcue = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 83)) == -1, :); % Rows with value -1 in column 65

% Filter rows based on the value in column 65 being 1 or -1
subset1_yellowcue = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 84)) == 1, :); % Rows with value 1 in column 65
subset2_yellowcue = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 84)) == -1, :); % Rows with value -1 in column 65

% Filter rows based on the value in column 65 being 1 or -1
subset1_redreward = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 85)) == 1, :); % Rows with value 1 in column 65
subset2_redreward = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 85)) == -1, :); % Rows with value -1 in column 65

% Filter rows based on the value in column 65 being 1 or -1
subset1_yellowairpuff = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 86)) == 1, :); % Rows with value 1 in column 65
subset2_yellowairpuff = cOFC_matrix2(cell2mat(cOFC_matrix2(:, 86)) == -1, :); % Rows with value -1 in column 65


% Sort each subset based on the values in column 67
% Convert column 67 values to a numeric array for sorting, if not already numeric

% Sorting subset1
[~, sortIndices1] = sort(cell2mat(subset1(:, 67)), 'descend');
subset1_sorted = subset1(sortIndices1, :);

% Sorting subset2
[~, sortIndices2] = sort(cell2mat(subset2(:, 67)), 'descend');
subset2_sorted = subset2(sortIndices2, :);

% Concatenate the sorted subsets
cOFC_matrix2_sorted = [subset1_sorted; subset2_sorted];

%%%% for airpuff
% Sort each subset based on the values in column 68
% Convert column 68 values to a numeric array for sorting, if not already numeric

% Sorting subset1_air
[~, sortIndices1_air] = sort(cell2mat(subset1_air(:, 68)), 'descend');
subset1_sorted_air = subset1_air(sortIndices1_air, :);

% Sorting subset2_air
[~, sortIndices2_air] = sort(cell2mat(subset2_air(:, 68)), 'descend');
subset2_sorted_air = subset2_air(sortIndices2_air, :);

% Concatenate the sorted subsets
cOFC_matrix2_sorted_air = [subset1_sorted_air; subset2_sorted_air];
%%%%

% Sorting subset1_smallreward
[~, sortIndices1_smallreward] = sort(cell2mat(subset1_smallreward(:, 69)), 'descend');
subset1_sorted_smallreward = subset1_smallreward(sortIndices1_smallreward, :);

% Sorting subset2_smallreward
[~, sortIndices2_smallreward] = sort(cell2mat(subset2_smallreward(:, 69)), 'descend');
subset2_sorted_smallreward = subset2_smallreward(sortIndices2_smallreward, :);

% Concatenate the sorted subsets
cOFC_matrix2_sorted_smallreward = [subset1_sorted_smallreward; subset2_sorted_smallreward];
%%%%

%%%% for big reward
% Sort each subset based on the values in column 76
% Convert column 76 values to a numeric array for sorting, if not already numeric

% Sorting subset1_bigreward
[~, sortIndices1_bigreward] = sort(cell2mat(subset1_bigreward(:, 76)), 'descend');
subset1_sorted_bigreward = subset1_bigreward(sortIndices1_bigreward, :);

% Sorting subset2_bigreward
[~, sortIndices2_bigreward] = sort(cell2mat(subset2_bigreward(:, 76)), 'descend');
subset2_sorted_bigreward = subset2_bigreward(sortIndices2_bigreward, :);

% Concatenate the sorted subsets
cOFC_matrix2_sorted_bigreward = [subset1_sorted_bigreward; subset2_sorted_bigreward];
%%%%
%%%% for big reward
% Sort each subset based on the values in column 76
% Convert column 76 values to a numeric array for sorting, if not already numeric

% Sorting subset1_redcue
[~, sortIndices1_redcue] = sort(cell2mat(subset1_redcue(:, 72)), 'descend');
subset1_sorted_redcue = subset1_redcue(sortIndices1_redcue, :);

% Sorting subset2_redcue
[~, sortIndices2_redcue] = sort(cell2mat(subset2_redcue(:, 72)), 'descend');
subset2_sorted_redcue = subset2_redcue(sortIndices2_redcue, :);

% Concatenate the sorted subsets
cOFC_matrix2_sorted_redcue = [subset1_sorted_redcue; subset2_sorted_redcue];
%%%%

%%%% for big reward
% Sort each subset based on the values in column 76
% Convert column 76 values to a numeric array for sorting, if not already numeric

% Sorting subset1_yellowcue
[~, sortIndices1_yellowcue] = sort(cell2mat(subset1_yellowcue(:, 73)), 'descend');
subset1_sorted_yellowcue = subset1_yellowcue(sortIndices1_yellowcue, :);

% Sorting subset2_yellowcue
[~, sortIndices2_yellowcue] = sort(cell2mat(subset2_yellowcue(:, 73)), 'descend');
subset2_sorted_yellowcue = subset2_yellowcue(sortIndices2_yellowcue, :);

% Concatenate the sorted subsets
cOFC_matrix2_sorted_yellowcue = [subset1_sorted_yellowcue; subset2_sorted_yellowcue];
%%%%
%%%% for big reward
% Sort each subset based on the values in column 76
% Convert column 76 values to a numeric array for sorting, if not already numeric

% Sorting subset1_redreward
[~, sortIndices1_redreward] = sort(cell2mat(subset1_redreward(:, 74)), 'descend');
subset1_sorted_redreward = subset1_redreward(sortIndices1_redreward, :);

% Sorting subset2_redreward
[~, sortIndices2_redreward] = sort(cell2mat(subset2_redreward(:, 74)), 'descend');
subset2_sorted_redreward = subset2_redreward(sortIndices2_redreward, :);

% Concatenate the sorted subsets
cOFC_matrix2_sorted_redreward = [subset1_sorted_redreward; subset2_sorted_redreward];
%%%%
%%%% for big reward
% Sort each subset based on the values in column 76
% Convert column 76 values to a numeric array for sorting, if not already numeric

% Sorting subset1_yellowairpuff
[~, sortIndices1_yellowairpuff] = sort(cell2mat(subset1_yellowairpuff(:, 75)), 'descend');
subset1_sorted_yellowairpuff = subset1_yellowairpuff(sortIndices1_yellowairpuff, :);

% Sorting subset2_yellowairpuff
[~, sortIndices2_yellowairpuff] = sort(cell2mat(subset2_yellowairpuff(:, 75)), 'descend');
subset2_sorted_yellowairpuff = subset2_yellowairpuff(sortIndices2_yellowairpuff, :);

% Concatenate the sorted subsets
cOFC_matrix2_sorted_yellowairpuff = [subset1_sorted_yellowairpuff; subset2_sorted_yellowairpuff];
%%%%
onlyPlus=subset1_sorted;
% Check for NaN values in column 67 and create a logical index for rows to keep
rowsToKeep_onlyPlus = ~cellfun(@(x) isnan(x), onlyPlus(:, 67));

% Use the logical index to filter out rows with NaN values in column 67
onlyPlus_noNaN = onlyPlus(rowsToKeep_onlyPlus, :);
% Create a logical index for rows where the value in column 67 is greater than or equal to 0.05
rowsToKeep_onlyPlus2 = cell2mat(onlyPlus_noNaN(:, 67)) <= 0.05;

% Use the logical index to filter out rows
onlyPlus_noNaN = onlyPlus_noNaN(rowsToKeep_onlyPlus2, :);

%%%% for airpuff
onlyPlus_air=subset1_sorted_air;
% Check for NaN values in column 68 and create a logical index for rows to keep
rowsToKeep_onlyPlus_air = ~cellfun(@(x) isnan(x), onlyPlus_air(:, 68));

% Use the logical index to filter out rows with NaN values in column 68
onlyPlus_noNaN_air = onlyPlus_air(rowsToKeep_onlyPlus_air, :);
% Create a logical index for rows where the value in column 68 is greater than or equal to 0.05
rowsToKeep_onlyPlus2_air = cell2mat(onlyPlus_noNaN_air(:, 68)) <= 0.05;

% Use the logical index to filter out rows
onlyPlus_noNaN_air = onlyPlus_noNaN_air(rowsToKeep_onlyPlus2_air, :);

%%%%%% for small reward
onlyPlus_smallreward=subset1_sorted_smallreward;
% Check for NaN values in column 69 and create a logical index for rows to keep
rowsToKeep_onlyPlus_smallreward = ~cellfun(@(x) isnan(x), onlyPlus_smallreward(:, 69));

% Use the logical index to filter out rows with NaN values in column 69
onlyPlus_noNaN_smallreward = onlyPlus_smallreward(rowsToKeep_onlyPlus_smallreward, :);
% Create a logical index for rows where the value in column 69 is greater than or equal to 0.05
rowsToKeep_onlyPlus2_smallreward = cell2mat(onlyPlus_noNaN_smallreward(:, 69)) <= 0.05;

% Use the logical index to filter out rows
onlyPlus_noNaN_smallreward = onlyPlus_noNaN_smallreward(rowsToKeep_onlyPlus2_smallreward, :);

%%%%%% for big reward
onlyPlus_bigreward=subset1_sorted_bigreward;
% Check for NaN values in column 76 and create a logical index for rows to keep
rowsToKeep_onlyPlus_bigreward = ~cellfun(@(x) isnan(x), onlyPlus_bigreward(:, 76));

% Use the logical index to filter out rows with NaN values in column 76
onlyPlus_noNaN_bigreward = onlyPlus_bigreward(rowsToKeep_onlyPlus_bigreward, :);
% Create a logical index for rows where the value in column 76 is greater than or equal to 0.05
rowsToKeep_onlyPlus2_bigreward = cell2mat(onlyPlus_noNaN_bigreward(:, 76)) <= 0.05;

% Use the logical index to filter out rows
onlyPlus_noNaN_bigreward = onlyPlus_noNaN_bigreward(rowsToKeep_onlyPlus2_bigreward, :);

%%%%%% for big reward
onlyPlus_redcue=subset1_sorted_redcue;
% Check for NaN values in column 76 and create a logical index for rows to keep
rowsToKeep_onlyPlus_redcue = ~cellfun(@(x) isnan(x), onlyPlus_redcue(:, 72));

% Use the logical index to filter out rows with NaN values in column 76
onlyPlus_noNaN_redcue = onlyPlus_redcue(rowsToKeep_onlyPlus_redcue, :);
% Create a logical index for rows where the value in column 76 is greater than or equal to 0.05
rowsToKeep_onlyPlus2_redcue = cell2mat(onlyPlus_noNaN_redcue(:, 72)) <= 0.05;

% Use the logical index to filter out rows
onlyPlus_noNaN_redcue = onlyPlus_noNaN_redcue(rowsToKeep_onlyPlus2_redcue, :);

%%%%%% for big reward
onlyPlus_yellowcue=subset1_sorted_yellowcue;
% Check for NaN values in column 76 and create a logical index for rows to keep
rowsToKeep_onlyPlus_yellowcue = ~cellfun(@(x) isnan(x), onlyPlus_yellowcue(:, 73));

% Use the logical index to filter out rows with NaN values in column 76
onlyPlus_noNaN_yellowcue = onlyPlus_yellowcue(rowsToKeep_onlyPlus_yellowcue, :);
% Create a logical index for rows where the value in column 76 is greater than or equal to 0.05
rowsToKeep_onlyPlus2_yellowcue = cell2mat(onlyPlus_noNaN_yellowcue(:, 73)) <= 0.05;

% Use the logical index to filter out rows
onlyPlus_noNaN_yellowcue = onlyPlus_noNaN_yellowcue(rowsToKeep_onlyPlus2_yellowcue, :);

%%%%%% for big reward
onlyPlus_redreward=subset1_sorted_redreward;
% Check for NaN values in column 76 and create a logical index for rows to keep
rowsToKeep_onlyPlus_redreward = ~cellfun(@(x) isnan(x), onlyPlus_redreward(:, 74));

% Use the logical index to filter out rows with NaN values in column 76
onlyPlus_noNaN_redreward = onlyPlus_redreward(rowsToKeep_onlyPlus_redreward, :);
% Create a logical index for rows where the value in column 76 is greater than or equal to 0.05
rowsToKeep_onlyPlus2_redreward = cell2mat(onlyPlus_noNaN_redreward(:, 74)) <= 0.05;

% Use the logical index to filter out rows
onlyPlus_noNaN_redreward = onlyPlus_noNaN_redreward(rowsToKeep_onlyPlus2_redreward, :);

%%%%%% for big reward
onlyPlus_yellowairpuff=subset1_sorted_yellowairpuff;
% Check for NaN values in column 76 and create a logical index for rows to keep
rowsToKeep_onlyPlus_yellowairpuff = ~cellfun(@(x) isnan(x), onlyPlus_yellowairpuff(:, 75));

% Use the logical index to filter out rows with NaN values in column 76
onlyPlus_noNaN_yellowairpuff = onlyPlus_yellowairpuff(rowsToKeep_onlyPlus_yellowairpuff, :);
% Create a logical index for rows where the value in column 76 is greater than or equal to 0.05
rowsToKeep_onlyPlus2_yellowairpuff = cell2mat(onlyPlus_noNaN_yellowairpuff(:, 75)) <= 0.05;

% Use the logical index to filter out rows
onlyPlus_noNaN_yellowairpuff = onlyPlus_noNaN_yellowairpuff(rowsToKeep_onlyPlus2_yellowairpuff, :);

% Check for NaN values in column 67 and create a logical index for rows to keep
rowsToKeep = ~cellfun(@(x) isnan(x), cOFC_matrix2_sorted(:, 67));

% Use the logical index to filter out rows with NaN values in column 67
cOFC_matrix2_noNaN = cOFC_matrix2_sorted(rowsToKeep, :);

cOFC_matrix2_noNaN_excluded=cOFC_matrix2_noNaN;

% Assuming cOFC_matrix2_noNaN is the matrix

% Create a logical index for rows where the value in column 67 is greater than or equal to 0.05
rowsToKeep = cell2mat(cOFC_matrix2_noNaN(:, 67)) <= 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN = cOFC_matrix2_noNaN(rowsToKeep, :);


% Create a logical index for rows where the value in column 67 is greater than or equal to 0.05
rowsToKeep_excluded = cell2mat(cOFC_matrix2_noNaN_excluded(:, 67)) > 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_excluded = cOFC_matrix2_noNaN_excluded(rowsToKeep_excluded, :);

%%%%% for airpuff
% Check for NaN values in column 68 and create a logical index for rows to keep
rowsToKeep_air = ~cellfun(@(x) isnan(x), cOFC_matrix2_sorted_air(:, 68));

% Use the logical index to filter out rows with NaN values in column 68
cOFC_matrix2_noNaN_air = cOFC_matrix2_sorted_air(rowsToKeep_air, :);

cOFC_matrix2_noNaN_excluded_air=cOFC_matrix2_noNaN_air;

% Assuming cOFC_matrix2_noNaN_air is the matrix

% Create a logical index for rows where the value in column 68 is greater than or equal to 0.05
rowsToKeep_air = cell2mat(cOFC_matrix2_noNaN_air(:, 68)) <= 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_air = cOFC_matrix2_noNaN_air(rowsToKeep_air, :);


% Create a logical index for rows where the value in column 68 is greater than or equal to 0.05
rowsToKeep_excluded_air = cell2mat(cOFC_matrix2_noNaN_excluded_air(:, 68)) > 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_excluded_air = cOFC_matrix2_noNaN_excluded_air(rowsToKeep_excluded_air, :);

%%%%% for small reward

% Check for NaN values in column 69 and create a logical index for rows to keep
rowsToKeep_smallreward = ~cellfun(@(x) isnan(x), cOFC_matrix2_sorted_smallreward(:, 69));

% Use the logical index to filter out rows with NaN values in column 69
cOFC_matrix2_noNaN_smallreward = cOFC_matrix2_sorted_smallreward(rowsToKeep_smallreward, :);

cOFC_matrix2_noNaN_excluded_smallreward=cOFC_matrix2_noNaN_smallreward;

% Assuming cOFC_matrix2_noNaN_smallreward is the matrix

% Create a logical index for rows where the value in column 69 is greater than or equal to 0.05
rowsToKeep_smallreward = cell2mat(cOFC_matrix2_noNaN_smallreward(:, 69)) <= 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_smallreward = cOFC_matrix2_noNaN_smallreward(rowsToKeep_smallreward, :);


% Create a logical index for rows where the value in column 69 is greater than or equal to 0.05
rowsToKeep_excluded_smallreward = cell2mat(cOFC_matrix2_noNaN_excluded_smallreward(:, 69)) > 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_excluded_smallreward = cOFC_matrix2_noNaN_excluded_smallreward(rowsToKeep_excluded_smallreward, :);

%%%%%%%%%%%%%%%
%%%%% for big reward

% Check for NaN values in column 76 and create a logical index for rows to keep
rowsToKeep_bigreward = ~cellfun(@(x) isnan(x), cOFC_matrix2_sorted_bigreward(:, 76));

% Use the logical index to filter out rows with NaN values in column 76
cOFC_matrix2_noNaN_bigreward = cOFC_matrix2_sorted_bigreward(rowsToKeep_bigreward, :);

cOFC_matrix2_noNaN_excluded_bigreward=cOFC_matrix2_noNaN_bigreward;

% Assuming cOFC_matrix2_noNaN_bigreward is the matrix

% Create a logical index for rows where the value in column 76 is greater than or equal to 0.05
rowsToKeep_bigreward = cell2mat(cOFC_matrix2_noNaN_bigreward(:, 76)) <= 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_bigreward = cOFC_matrix2_noNaN_bigreward(rowsToKeep_bigreward, :);


% Create a logical index for rows where the value in column 76 is greater than or equal to 0.05
rowsToKeep_excluded_bigreward = cell2mat(cOFC_matrix2_noNaN_excluded_bigreward(:, 76)) > 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_excluded_bigreward = cOFC_matrix2_noNaN_excluded_bigreward(rowsToKeep_excluded_bigreward, :);

%%%%%%%%%%%%%%%
%%%%% for big reward

% Check for NaN values in column 72 and create a logical index for rows to keep
rowsToKeep_redcue = ~cellfun(@(x) isnan(x), cOFC_matrix2_sorted_redcue(:, 72));

% Use the logical index to filter out rows with NaN values in column 72
cOFC_matrix2_noNaN_redcue = cOFC_matrix2_sorted_redcue(rowsToKeep_redcue, :);

cOFC_matrix2_noNaN_excluded_redcue=cOFC_matrix2_noNaN_redcue;

% Assuming cOFC_matrix2_noNaN_redcue is the matrix

% Create a logical index for rows where the value in column 72 is greater than or equal to 0.05
rowsToKeep_redcue = cell2mat(cOFC_matrix2_noNaN_redcue(:, 72)) <= 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_redcue = cOFC_matrix2_noNaN_redcue(rowsToKeep_redcue, :);


% Create a logical index for rows where the value in column 72 is greater than or equal to 0.05
rowsToKeep_excluded_redcue = cell2mat(cOFC_matrix2_noNaN_excluded_redcue(:, 72)) > 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_excluded_redcue = cOFC_matrix2_noNaN_excluded_redcue(rowsToKeep_excluded_redcue, :);

%%%%%%%%%%%%%%%
%%%%% for big reward

% Check for NaN values in column 73 and create a logical index for rows to keep
rowsToKeep_yellowcue = ~cellfun(@(x) isnan(x), cOFC_matrix2_sorted_yellowcue(:, 73));

% Use the logical index to filter out rows with NaN values in column 73
cOFC_matrix2_noNaN_yellowcue = cOFC_matrix2_sorted_yellowcue(rowsToKeep_yellowcue, :);

cOFC_matrix2_noNaN_excluded_yellowcue=cOFC_matrix2_noNaN_yellowcue;

% Assuming cOFC_matrix2_noNaN_yellowcue is the matrix

% Create a logical index for rows where the value in column 73 is greater than or equal to 0.05
rowsToKeep_yellowcue = cell2mat(cOFC_matrix2_noNaN_yellowcue(:, 73)) <= 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_yellowcue = cOFC_matrix2_noNaN_yellowcue(rowsToKeep_yellowcue, :);


% Create a logical index for rows where the value in column 73 is greater than or equal to 0.05
rowsToKeep_excluded_yellowcue = cell2mat(cOFC_matrix2_noNaN_excluded_yellowcue(:, 73)) > 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_excluded_yellowcue = cOFC_matrix2_noNaN_excluded_yellowcue(rowsToKeep_excluded_yellowcue, :);

%%%%%%%%%%%%%%%
%%%%% for big reward

% Check for NaN values in column 74 and create a logical index for rows to keep
rowsToKeep_redreward = ~cellfun(@(x) isnan(x), cOFC_matrix2_sorted_redreward(:, 74));

% Use the logical index to filter out rows with NaN values in column 74
cOFC_matrix2_noNaN_redreward = cOFC_matrix2_sorted_redreward(rowsToKeep_redreward, :);

cOFC_matrix2_noNaN_excluded_redreward=cOFC_matrix2_noNaN_redreward;

% Assuming cOFC_matrix2_noNaN_redreward is the matrix

% Create a logical index for rows where the value in column 74 is greater than or equal to 0.05
rowsToKeep_redreward = cell2mat(cOFC_matrix2_noNaN_redreward(:, 74)) <= 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_redreward = cOFC_matrix2_noNaN_redreward(rowsToKeep_redreward, :);


% Create a logical index for rows where the value in column 74 is greater than or equal to 0.05
rowsToKeep_excluded_redreward = cell2mat(cOFC_matrix2_noNaN_excluded_redreward(:, 74)) > 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_excluded_redreward = cOFC_matrix2_noNaN_excluded_redreward(rowsToKeep_excluded_redreward, :);

%%%%%%%%%%%%%%%
%%%%% for big reward

% Check for NaN values in column 75 and create a logical index for rows to keep
rowsToKeep_yellowairpuff = ~cellfun(@(x) isnan(x), cOFC_matrix2_sorted_yellowairpuff(:, 75));

% Use the logical index to filter out rows with NaN values in column 75
cOFC_matrix2_noNaN_yellowairpuff = cOFC_matrix2_sorted_yellowairpuff(rowsToKeep_yellowairpuff, :);

cOFC_matrix2_noNaN_excluded_yellowairpuff=cOFC_matrix2_noNaN_yellowairpuff;

% Assuming cOFC_matrix2_noNaN_yellowairpuff is the matrix

% Create a logical index for rows where the value in column 75 is greater than or equal to 0.05
rowsToKeep_yellowairpuff = cell2mat(cOFC_matrix2_noNaN_yellowairpuff(:, 75)) <= 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_yellowairpuff = cOFC_matrix2_noNaN_yellowairpuff(rowsToKeep_yellowairpuff, :);


% Create a logical index for rows where the value in column 75 is greater than or equal to 0.05
rowsToKeep_excluded_yellowairpuff = cell2mat(cOFC_matrix2_noNaN_excluded_yellowairpuff(:, 75)) > 0.05;

% Use the logical index to filter out rows
cOFC_matrix2_noNaN_excluded_yellowairpuff = cOFC_matrix2_noNaN_excluded_yellowairpuff(rowsToKeep_excluded_yellowairpuff, :);

% Initialize an empty matrix for concatenation
cOFC_spikes_CuePeriod = [];

% Iterate through each row in column 63 of cOFC_matrix2_noNaN
for i = 1:size(cOFC_matrix2_noNaN, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray = cOFC_matrix2_noNaN{i, 63};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_CuePeriod = [cOFC_spikes_CuePeriod; currentArray];
end

% cOFC_spikes_CuePeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_CuePeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_CuePeriod_normalized = zeros(size(cOFC_spikes_CuePeriod));

% Normalize each row
for i = 1:size(cOFC_spikes_CuePeriod, 1)
    row = cOFC_spikes_CuePeriod(i, :);
    normalizedRow = (row - mean(row)) / std(row);
    cOFC_spikes_CuePeriod_normalized(i, :) = normalizedRow;
end

% cOFC_spikes_CuePeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_CuePeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN = all(isnan(cOFC_spikes_CuePeriod_normalized), 2);

% Find rows that are all zeros
rowsAllZero = all(cOFC_spikes_CuePeriod_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove = rowsAllNaN | rowsAllZero;

% Filter out these rows
cOFC_spikes_CuePeriod_normalized_filtered = cOFC_spikes_CuePeriod_normalized(~rowsToRemove, :);

%%%% the ones excluded
% Initialize an empty matrix for concatenation
cOFC_spikes_CuePeriod_excluded= [];

% Iterate through each row in column 63 of cOFC_matrix2_noNaN
for i = 1:size(cOFC_matrix2_noNaN_excluded, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded  = cOFC_matrix2_noNaN_excluded{i, 63};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_CuePeriod_excluded = [cOFC_spikes_CuePeriod_excluded; currentArray_excluded];
end

% cOFC_spikes_CuePeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_CuePeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_CuePeriod_excluded_normalized = zeros(size(cOFC_spikes_CuePeriod_excluded));

% Normalize each row
for i = 1:size(cOFC_spikes_CuePeriod_excluded, 1)
    row = cOFC_spikes_CuePeriod_excluded(i, :);
    normalizedRow_excluded = (row - mean(row)) / std(row);
    cOFC_spikes_CuePeriod_excluded_normalized(i, :) = normalizedRow_excluded;
end

% cOFC_spikes_CuePeriod_normalized now contains the normalized rows

% Find rows that are all NaN
rowsAllNaN = all(isnan(cOFC_spikes_CuePeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero = all(cOFC_spikes_CuePeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove = rowsAllNaN | rowsAllZero;

% Filter out these rows
cOFC_spikes_CuePeriod_excluded_normalized_filtered = cOFC_spikes_CuePeriod_excluded_normalized(~rowsToRemove, :);

%%%%%%%%%%%%%%%%%% for airpuff period
% Initialize an empty matrix for concatenation
cOFC_spikes_AirpuffPeriod = [];

% Iterate through each row in column 88 of cOFC_matrix2_noNaN_air
for i = 1:size(cOFC_matrix2_noNaN_air, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_air = cOFC_matrix2_noNaN_air{i, 88};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_AirpuffPeriod = [cOFC_spikes_AirpuffPeriod; currentArray_air];
end

% cOFC_spikes_AirpuffPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_AirpuffPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_AirpuffPeriod_normalized = zeros(size(cOFC_spikes_AirpuffPeriod));

% Normalize each row
for i = 1:size(cOFC_spikes_AirpuffPeriod, 1)
    row = cOFC_spikes_AirpuffPeriod(i, :);
    normalizedRow_air = (row - mean(row)) / std(row);
    cOFC_spikes_AirpuffPeriod_normalized(i, :) = normalizedRow_air;
end

% cOFC_spikes_AirpuffPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_AirpuffPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_air = all(isnan(cOFC_spikes_AirpuffPeriod_normalized), 2);

% Find rows that are all zeros
rowsAllZero_air = all(cOFC_spikes_AirpuffPeriod_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_air = rowsAllNaN_air | rowsAllZero_air;

% Filter out these rows
cOFC_spikes_AirpuffPeriod_normalized_filtered = cOFC_spikes_AirpuffPeriod_normalized(~rowsToRemove_air, :);

%%%% the ones excluded
% Initialize an empty matrix for concatenation
cOFC_spikes_AirpuffPeriod_excluded= [];

% Iterate through each row in column 88 of cOFC_matrix2_noNaN
for i = 1:size(cOFC_matrix2_noNaN_excluded_air, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_air  = cOFC_matrix2_noNaN_excluded_air{i, 88};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_AirpuffPeriod_excluded = [cOFC_spikes_AirpuffPeriod_excluded; currentArray_excluded_air];
end

% Initialize a matrix to hold the normalized data
cOFC_spikes_AirpuffPeriod_excluded_normalized = zeros(size(cOFC_spikes_AirpuffPeriod_excluded));

% Normalize each row
for i = 1:size(cOFC_spikes_AirpuffPeriod_excluded, 1)
    row = cOFC_spikes_AirpuffPeriod_excluded(i, :);
    normalizedRow_excluded_air = (row - mean(row)) / std(row);
    cOFC_spikes_AirpuffPeriod_excluded_normalized(i, :) = normalizedRow_excluded_air;
end

% cOFC_spikes_AirpuffPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_AirpuffPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_air = all(isnan(cOFC_spikes_AirpuffPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_air = all(cOFC_spikes_AirpuffPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_air = rowsAllNaN_air | rowsAllZero_air;

% Filter out these rows
cOFC_spikes_AirpuffPeriod_excluded_normalized_filtered = cOFC_spikes_AirpuffPeriod_excluded_normalized(~rowsToRemove_air, :);


%%%%%% for small reward
%%%%%%%%%%%%%%%%%% for airpuff period
% Initialize an empty matrix for concatenation
cOFC_spikes_SmallRewardPeriod = [];

% Iterate through each row in column 89 of cOFC_matrix2_noNaN_smallreward
for i = 1:size(cOFC_matrix2_noNaN_smallreward, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_smallreward = cOFC_matrix2_noNaN_smallreward{i, 89};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_SmallRewardPeriod = [cOFC_spikes_SmallRewardPeriod; currentArray_smallreward];
end

% cOFC_spikes_SmallRewardPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_SmallRewardPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_SmallRewardPeriod_normalized = zeros(size(cOFC_spikes_SmallRewardPeriod));

% Normalize each row
for i = 1:size(cOFC_spikes_SmallRewardPeriod, 1)
    row = cOFC_spikes_SmallRewardPeriod(i, :);
    normalizedRow_smallreward = (row - mean(row)) / std(row);
    cOFC_spikes_SmallRewardPeriod_normalized(i, :) = normalizedRow_smallreward;
end

% cOFC_spikes_SmallRewardPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_SmallRewardPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_smallreward = all(isnan(cOFC_spikes_SmallRewardPeriod_normalized), 2);

% Find rows that are all zeros
rowsAllZero_smallreward = all(cOFC_spikes_SmallRewardPeriod_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_smallreward = rowsAllNaN_smallreward | rowsAllZero_smallreward;

% Filter out these rows
cOFC_spikes_SmallRewardPeriod_normalized_filtered = cOFC_spikes_SmallRewardPeriod_normalized(~rowsToRemove_smallreward, :);



%%%% the ones excluded
% Initialize an empty matrix for concatenation
cOFC_spikes_SmallRewardPeriod_excluded= [];

% Iterate through each row in column 89 of cOFC_matrix2_noNaN
for i = 1:size(cOFC_matrix2_noNaN_excluded_smallreward, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_smallreward  = cOFC_matrix2_noNaN_excluded_smallreward{i, 89};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_SmallRewardPeriod_excluded = [cOFC_spikes_SmallRewardPeriod_excluded; currentArray_excluded_smallreward];
end

% cOFC_spikes_SmallRewardPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_SmallRewardPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_SmallRewardPeriod_excluded_normalized = zeros(size(cOFC_spikes_SmallRewardPeriod_excluded));

% Normalize each row
for i = 1:size(cOFC_spikes_SmallRewardPeriod_excluded, 1)
    row = cOFC_spikes_SmallRewardPeriod_excluded(i, :);
    normalizedRow_excluded_smallreward = (row - mean(row)) / std(row);
    cOFC_spikes_SmallRewardPeriod_excluded_normalized(i, :) = normalizedRow_excluded_smallreward;
end

% cOFC_spikes_SmallRewardPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_SmallRewardPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_smallreward = all(isnan(cOFC_spikes_SmallRewardPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_smallreward = all(cOFC_spikes_SmallRewardPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_smallreward = rowsAllNaN_smallreward | rowsAllZero_smallreward;

% Filter out these rows
cOFC_spikes_SmallRewardPeriod_excluded_normalized_filtered = cOFC_spikes_SmallRewardPeriod_excluded_normalized(~rowsToRemove_smallreward, :);



%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%% for big reward
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% for airpuff period
% Initialize an empty matrix for concatenation
cOFC_spikes_BigRewardPeriod = [];

% Iterate through each row in column 90 of cOFC_matrix2_noNaN_bigreward
for i = 1:size(cOFC_matrix2_noNaN_bigreward, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_bigreward = cOFC_matrix2_noNaN_bigreward{i, 90};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_BigRewardPeriod = [cOFC_spikes_BigRewardPeriod; currentArray_bigreward];
end

% cOFC_spikes_BigRewardPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_BigRewardPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_BigRewardPeriod_normalized = zeros(size(cOFC_spikes_BigRewardPeriod));

% Normalize each row
for i = 1:size(cOFC_spikes_BigRewardPeriod, 1)
    row = cOFC_spikes_BigRewardPeriod(i, :);
    normalizedRow_bigreward = (row - mean(row)) / std(row);
    cOFC_spikes_BigRewardPeriod_normalized(i, :) = normalizedRow_bigreward;
end

% cOFC_spikes_BigRewardPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_BigRewardPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_bigreward = all(isnan(cOFC_spikes_BigRewardPeriod_normalized), 2);

% Find rows that are all zeros
rowsAllZero_bigreward = all(cOFC_spikes_BigRewardPeriod_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_bigreward = rowsAllNaN_bigreward | rowsAllZero_bigreward;

% Filter out these rows
cOFC_spikes_BigRewardPeriod_normalized_filtered = cOFC_spikes_BigRewardPeriod_normalized(~rowsToRemove_bigreward, :);



%%%% the ones excluded
% Initialize an empty matrix for concatenation
cOFC_spikes_BigRewardPeriod_excluded= [];

% Iterate through each row in column 90 of cOFC_matrix2_noNaN
for i = 1:size(cOFC_matrix2_noNaN_excluded_bigreward, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_bigreward  = cOFC_matrix2_noNaN_excluded_bigreward{i, 90};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_BigRewardPeriod_excluded = [cOFC_spikes_BigRewardPeriod_excluded; currentArray_excluded_bigreward];
end

% cOFC_spikes_BigRewardPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_BigRewardPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_BigRewardPeriod_excluded_normalized = zeros(size(cOFC_spikes_BigRewardPeriod_excluded));

% Normalize each row
for i = 1:size(cOFC_spikes_BigRewardPeriod_excluded, 1)
    row = cOFC_spikes_BigRewardPeriod_excluded(i, :);
    normalizedRow_excluded_bigreward = (row - mean(row)) / std(row);
    cOFC_spikes_BigRewardPeriod_excluded_normalized(i, :) = normalizedRow_excluded_bigreward;
end

% cOFC_spikes_BigRewardPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_BigRewardPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_bigreward = all(isnan(cOFC_spikes_BigRewardPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_bigreward = all(cOFC_spikes_BigRewardPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_bigreward = rowsAllNaN_bigreward | rowsAllZero_bigreward;

% Filter out these rows
cOFC_spikes_BigRewardPeriod_excluded_normalized_filtered = cOFC_spikes_BigRewardPeriod_excluded_normalized(~rowsToRemove_bigreward, :);

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%% for big reward
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% for airpuff period
% Initialize an empty matrix for concatenation
cOFC_spikes_RedCuePeriod = [];

% Iterate through each row in column 91 of cOFC_matrix2_noNaN_redcue
for i = 1:size(cOFC_matrix2_noNaN_redcue, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_redcue = cOFC_matrix2_noNaN_redcue{i, 91};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_RedCuePeriod = [cOFC_spikes_RedCuePeriod; currentArray_redcue];
end

% cOFC_spikes_RedCuePeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_RedCuePeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_RedCuePeriod_normalized = zeros(size(cOFC_spikes_RedCuePeriod));

% Normalize each row
for i = 1:size(cOFC_spikes_RedCuePeriod, 1)
    row = cOFC_spikes_RedCuePeriod(i, :);
    normalizedRow_redcue = (row - mean(row)) / std(row);
    cOFC_spikes_RedCuePeriod_normalized(i, :) = normalizedRow_redcue;
end

% cOFC_spikes_RedCuePeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_RedCuePeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_redcue = all(isnan(cOFC_spikes_RedCuePeriod_normalized), 2);

% Find rows that are all zeros
rowsAllZero_redcue = all(cOFC_spikes_RedCuePeriod_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_redcue = rowsAllNaN_redcue | rowsAllZero_redcue;

% Filter out these rows
cOFC_spikes_RedCuePeriod_normalized_filtered = cOFC_spikes_RedCuePeriod_normalized(~rowsToRemove_redcue, :);



%%%% the ones excluded
% Initialize an empty matrix for concatenation
cOFC_spikes_RedCuePeriod_excluded= [];

% Iterate through each row in column 91 of cOFC_matrix2_noNaN
for i = 1:size(cOFC_matrix2_noNaN_excluded_redcue, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_redcue  = cOFC_matrix2_noNaN_excluded_redcue{i, 91};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_RedCuePeriod_excluded = [cOFC_spikes_RedCuePeriod_excluded; currentArray_excluded_redcue];
end

% cOFC_spikes_RedCuePeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_RedCuePeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_RedCuePeriod_excluded_normalized = zeros(size(cOFC_spikes_RedCuePeriod_excluded));

% Normalize each row
for i = 1:size(cOFC_spikes_RedCuePeriod_excluded, 1)
    row = cOFC_spikes_RedCuePeriod_excluded(i, :);
    normalizedRow_excluded_redcue = (row - mean(row)) / std(row);
    cOFC_spikes_RedCuePeriod_excluded_normalized(i, :) = normalizedRow_excluded_redcue;
end

% cOFC_spikes_RedCuePeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_RedCuePeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_redcue = all(isnan(cOFC_spikes_RedCuePeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_redcue = all(cOFC_spikes_RedCuePeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_redcue = rowsAllNaN_redcue | rowsAllZero_redcue;

% Filter out these rows
cOFC_spikes_RedCuePeriod_excluded_normalized_filtered = cOFC_spikes_RedCuePeriod_excluded_normalized(~rowsToRemove_redcue, :);

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%% for big reward
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% for airpuff period
% Initialize an empty matrix for concatenation
cOFC_spikes_YellowCuePeriod = [];

% Iterate through each row in column 92 of cOFC_matrix2_noNaN_yellowcue
for i = 1:size(cOFC_matrix2_noNaN_yellowcue, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_yellowcue = cOFC_matrix2_noNaN_yellowcue{i, 92};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_YellowCuePeriod = [cOFC_spikes_YellowCuePeriod; currentArray_yellowcue];
end

% cOFC_spikes_YellowCuePeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_YellowCuePeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_YellowCuePeriod_normalized = zeros(size(cOFC_spikes_YellowCuePeriod));

% Normalize each row
for i = 1:size(cOFC_spikes_YellowCuePeriod, 1)
    row = cOFC_spikes_YellowCuePeriod(i, :);
    normalizedRow_yellowcue = (row - mean(row)) / std(row);
    cOFC_spikes_YellowCuePeriod_normalized(i, :) = normalizedRow_yellowcue;
end

% cOFC_spikes_YellowCuePeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_YellowCuePeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_yellowcue = all(isnan(cOFC_spikes_YellowCuePeriod_normalized), 2);

% Find rows that are all zeros
rowsAllZero_yellowcue = all(cOFC_spikes_YellowCuePeriod_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_yellowcue = rowsAllNaN_yellowcue | rowsAllZero_yellowcue;

% Filter out these rows
cOFC_spikes_YellowCuePeriod_normalized_filtered = cOFC_spikes_YellowCuePeriod_normalized(~rowsToRemove_yellowcue, :);



%%%% the ones excluded
% Initialize an empty matrix for concatenation
cOFC_spikes_YellowCuePeriod_excluded= [];

% Iterate through each row in column 92 of cOFC_matrix2_noNaN
for i = 1:size(cOFC_matrix2_noNaN_excluded_yellowcue, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_yellowcue  = cOFC_matrix2_noNaN_excluded_yellowcue{i, 92};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_YellowCuePeriod_excluded = [cOFC_spikes_YellowCuePeriod_excluded; currentArray_excluded_yellowcue];
end

% cOFC_spikes_YellowCuePeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_YellowCuePeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_YellowCuePeriod_excluded_normalized = zeros(size(cOFC_spikes_YellowCuePeriod_excluded));

% Normalize each row
for i = 1:size(cOFC_spikes_YellowCuePeriod_excluded, 1)
    row = cOFC_spikes_YellowCuePeriod_excluded(i, :);
    normalizedRow_excluded_yellowcue = (row - mean(row)) / std(row);
    cOFC_spikes_YellowCuePeriod_excluded_normalized(i, :) = normalizedRow_excluded_yellowcue;
end

% cOFC_spikes_YellowCuePeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_YellowCuePeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_yellowcue = all(isnan(cOFC_spikes_YellowCuePeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_yellowcue = all(cOFC_spikes_YellowCuePeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_yellowcue = rowsAllNaN_yellowcue | rowsAllZero_yellowcue;

% Filter out these rows
cOFC_spikes_YellowCuePeriod_excluded_normalized_filtered = cOFC_spikes_YellowCuePeriod_excluded_normalized(~rowsToRemove_yellowcue, :);

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%% for big reward
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% for airpuff period
% Initialize an empty matrix for concatenation
cOFC_spikes_RedRewardPeriod = [];

% Iterate through each row in column 93 of cOFC_matrix2_noNaN_redreward
for i = 1:size(cOFC_matrix2_noNaN_redreward, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_redreward = cOFC_matrix2_noNaN_redreward{i, 93};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_RedRewardPeriod = [cOFC_spikes_RedRewardPeriod; currentArray_redreward];
end

% cOFC_spikes_RedRewardPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_RedRewardPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_RedRewardPeriod_normalized = zeros(size(cOFC_spikes_RedRewardPeriod));

% Normalize each row
for i = 1:size(cOFC_spikes_RedRewardPeriod, 1)
    row = cOFC_spikes_RedRewardPeriod(i, :);
    normalizedRow_redreward = (row - mean(row)) / std(row);
    cOFC_spikes_RedRewardPeriod_normalized(i, :) = normalizedRow_redreward;
end

% cOFC_spikes_RedRewardPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_RedRewardPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_redreward = all(isnan(cOFC_spikes_RedRewardPeriod_normalized), 2);

% Find rows that are all zeros
rowsAllZero_redreward = all(cOFC_spikes_RedRewardPeriod_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_redreward = rowsAllNaN_redreward | rowsAllZero_redreward;

% Filter out these rows
cOFC_spikes_RedRewardPeriod_normalized_filtered = cOFC_spikes_RedRewardPeriod_normalized(~rowsToRemove_redreward, :);



%%%% the ones excluded
% Initialize an empty matrix for concatenation
cOFC_spikes_RedRewardPeriod_excluded= [];

% Iterate through each row in column 93 of cOFC_matrix2_noNaN
for i = 1:size(cOFC_matrix2_noNaN_excluded_redreward, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_redreward  = cOFC_matrix2_noNaN_excluded_redreward{i, 93};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_RedRewardPeriod_excluded = [cOFC_spikes_RedRewardPeriod_excluded; currentArray_excluded_redreward];
end

% cOFC_spikes_RedRewardPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_RedRewardPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_RedRewardPeriod_excluded_normalized = zeros(size(cOFC_spikes_RedRewardPeriod_excluded));

% Normalize each row
for i = 1:size(cOFC_spikes_RedRewardPeriod_excluded, 1)
    row = cOFC_spikes_RedRewardPeriod_excluded(i, :);
    normalizedRow_excluded_redreward = (row - mean(row)) / std(row);
    cOFC_spikes_RedRewardPeriod_excluded_normalized(i, :) = normalizedRow_excluded_redreward;
end

% cOFC_spikes_RedRewardPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_RedRewardPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_redreward = all(isnan(cOFC_spikes_RedRewardPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_redreward = all(cOFC_spikes_RedRewardPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_redreward = rowsAllNaN_redreward | rowsAllZero_redreward;

% Filter out these rows
cOFC_spikes_RedRewardPeriod_excluded_normalized_filtered = cOFC_spikes_RedRewardPeriod_excluded_normalized(~rowsToRemove_redreward, :);

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%% for big reward
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% for airpuff period
% Initialize an empty matrix for concatenation
cOFC_spikes_YellowAirpuffPeriod = [];

% Iterate through each row in column 94 of cOFC_matrix2_noNaN_yellowairpuff
for i = 1:size(cOFC_matrix2_noNaN_yellowairpuff, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_yellowairpuff = cOFC_matrix2_noNaN_yellowairpuff{i, 94};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_YellowAirpuffPeriod = [cOFC_spikes_YellowAirpuffPeriod; currentArray_yellowairpuff];
end

% cOFC_spikes_YellowAirpuffPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_YellowAirpuffPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_YellowAirpuffPeriod_normalized = zeros(size(cOFC_spikes_YellowAirpuffPeriod));

% Normalize each row
for i = 1:size(cOFC_spikes_YellowAirpuffPeriod, 1)
    row = cOFC_spikes_YellowAirpuffPeriod(i, :);
    normalizedRow_yellowairpuff = (row - mean(row)) / std(row);
    cOFC_spikes_YellowAirpuffPeriod_normalized(i, :) = normalizedRow_yellowairpuff;
end

% cOFC_spikes_YellowAirpuffPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_YellowAirpuffPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_yellowairpuff = all(isnan(cOFC_spikes_YellowAirpuffPeriod_normalized), 2);

% Find rows that are all zeros
rowsAllZero_yellowairpuff = all(cOFC_spikes_YellowAirpuffPeriod_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_yellowairpuff = rowsAllNaN_yellowairpuff | rowsAllZero_yellowairpuff;

% Filter out these rows
cOFC_spikes_YellowAirpuffPeriod_normalized_filtered = cOFC_spikes_YellowAirpuffPeriod_normalized(~rowsToRemove_yellowairpuff, :);



%%%% the ones excluded
% Initialize an empty matrix for concatenation
cOFC_spikes_YellowAirpuffPeriod_excluded= [];

% Iterate through each row in column 94 of cOFC_matrix2_noNaN
for i = 1:size(cOFC_matrix2_noNaN_excluded_yellowairpuff, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_yellowairpuff  = cOFC_matrix2_noNaN_excluded_yellowairpuff{i, 94};
    
    % Concatenate this array to the growing matrix
    cOFC_spikes_YellowAirpuffPeriod_excluded = [cOFC_spikes_YellowAirpuffPeriod_excluded; currentArray_excluded_yellowairpuff];
end

% cOFC_spikes_YellowAirpuffPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming cOFC_spikes_YellowAirpuffPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
cOFC_spikes_YellowAirpuffPeriod_excluded_normalized = zeros(size(cOFC_spikes_YellowAirpuffPeriod_excluded));

% Normalize each row
for i = 1:size(cOFC_spikes_YellowAirpuffPeriod_excluded, 1)
    row = cOFC_spikes_YellowAirpuffPeriod_excluded(i, :);
    normalizedRow_excluded_yellowairpuff = (row - mean(row)) / std(row);
    cOFC_spikes_YellowAirpuffPeriod_excluded_normalized(i, :) = normalizedRow_excluded_yellowairpuff;
end

% cOFC_spikes_YellowAirpuffPeriod_normalized now contains the normalized rows

% Assuming cOFC_spikes_YellowAirpuffPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_yellowairpuff = all(isnan(cOFC_spikes_YellowAirpuffPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_yellowairpuff = all(cOFC_spikes_YellowAirpuffPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_yellowairpuff = rowsAllNaN_yellowairpuff | rowsAllZero_yellowairpuff;

% Filter out these rows
cOFC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered = cOFC_spikes_YellowAirpuffPeriod_excluded_normalized(~rowsToRemove_yellowairpuff, :);



% % Assuming cOFC_spikes_CuePeriod is already defined and contains the data
% 
% % Normalize each row of cOFC_spikes_CuePeriod
% normalized_cOFC_spikes_CuePeriod = cOFC_spikes_CuePeriod;
% 
% for i = 1:size(cOFC_spikes_CuePeriod, 1)
%     row = cOFC_spikes_CuePeriod(i, :);
%     normalized_cOFC_spikes_CuePeriod(i, :) = (row - min(row)) / (max(row) - min(row));
% end
% 
% % Create a heatmap of the normalized data

comeon=[cOFC_spikes_CuePeriod_excluded_normalized_filtered;cOFC_spikes_CuePeriod_normalized_filtered];
comeon_air=[cOFC_spikes_AirpuffPeriod_excluded_normalized_filtered;cOFC_spikes_AirpuffPeriod_normalized_filtered];
comeon_smallreward=[cOFC_spikes_SmallRewardPeriod_excluded_normalized_filtered;cOFC_spikes_SmallRewardPeriod_normalized_filtered];
comeon_bigreward=[cOFC_spikes_BigRewardPeriod_excluded_normalized_filtered;cOFC_spikes_BigRewardPeriod_normalized_filtered];
comeon_redcue=[cOFC_spikes_RedCuePeriod_excluded_normalized_filtered;cOFC_spikes_RedCuePeriod_normalized_filtered];
comeon_yellowcue=[cOFC_spikes_YellowCuePeriod_excluded_normalized_filtered;cOFC_spikes_YellowCuePeriod_normalized_filtered];
comeon_redreward=[cOFC_spikes_RedRewardPeriod_excluded_normalized_filtered;cOFC_spikes_RedRewardPeriod_normalized_filtered];
comeon_yellowairpuff=[cOFC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered;cOFC_spikes_YellowAirpuffPeriod_normalized_filtered];

%%%%%%%%%% for pACC
% Filter for pACC data
pACC_matrix2 = matrix2(strcmp(matrix2(:, 6), 'pACC'), :);

% Filter rows based on the value in column 64 being 1 or -1 and sort
subset1_pACC = pACC_matrix2(cell2mat(pACC_matrix2(:, 64)) == 1, :);
subset2_pACC = pACC_matrix2(cell2mat(pACC_matrix2(:, 64)) == -1, :);

% Sorting subsets based on column 67
[~, sortIndices1_pACC] = sort(cell2mat(subset1_pACC(:, 67)), 'descend');
subset1_sorted_pACC = subset1_pACC(sortIndices1_pACC, :);

[~, sortIndices2_pACC] = sort(cell2mat(subset2_pACC(:, 67)), 'descend');
subset2_sorted_pACC = subset2_pACC(sortIndices2_pACC, :);

% Concatenate the sorted subsets
pACC_matrix2_sorted = [subset1_sorted_pACC; subset2_sorted_pACC];


% Remove rows with NaN in column 67
rowsToKeep_pACC = ~cellfun(@(x) isnan(x), pACC_matrix2_sorted(:, 67));
pACC_matrix2_noNaN = pACC_matrix2_sorted(rowsToKeep_pACC, :);
pACC_matrix2_noNaN_excluded = pACC_matrix2_noNaN;

% Logical index for rows to keep in pACC_matrix2_noNaN (values in column 67 >= 0.05)
keepRows_pACC = cell2mat(pACC_matrix2_noNaN(:, 67)) <= 0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN = pACC_matrix2_noNaN(keepRows_pACC, :);



% Logical index for rows to keep in pACC_matrix2_noNaN (values in column 67 >= 0.05)
keepRows_pACC_excluded = cell2mat(pACC_matrix2_noNaN_excluded(:, 67)) >0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_excluded = pACC_matrix2_noNaN_excluded(keepRows_pACC_excluded, :);

% Concatenate data from column 63 into a new matrix
pACC_spikes_CuePeriod = [];
for i = 1:size(pACC_matrix2_noNaN, 1)
    currentArray_pACC = pACC_matrix2_noNaN{i, 63};
    pACC_spikes_CuePeriod = [pACC_spikes_CuePeriod; currentArray_pACC];
end

% Normalize each row of pACC_spikes_CuePeriod
pACC_spikes_CuePeriod_normalized = zeros(size(pACC_spikes_CuePeriod));
for i = 1:size(pACC_spikes_CuePeriod, 1)
    row = pACC_spikes_CuePeriod(i, :);
    normalizedRow = (row - mean(row)) / std(row);
    pACC_spikes_CuePeriod_normalized(i, :) = normalizedRow;
end

% Remove rows that are all NaN or all zero
rowsAllNaN_pACC = all(isnan(pACC_spikes_CuePeriod_normalized), 2);
rowsAllZero_pACC = all(pACC_spikes_CuePeriod_normalized == 0, 2);
rowsToRemove_pACC = rowsAllNaN_pACC | rowsAllZero_pACC;
pACC_spikes_CuePeriod_normalized_filtered = pACC_spikes_CuePeriod_normalized(~rowsToRemove_pACC, :);


%%%% the ones excluded
% Initialize an empty matrix for concatenation
pACC_spikes_CuePeriod_excluded= [];

% Iterate through each row in column 63 of pACC_matrix2_noNaN
for i = 1:size(pACC_matrix2_noNaN_excluded, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded  = pACC_matrix2_noNaN_excluded{i, 63};
    
    % Concatenate this array to the growing matrix
    pACC_spikes_CuePeriod_excluded = [pACC_spikes_CuePeriod_excluded; currentArray_excluded];
end

% pACC_spikes_CuePeriod now contains all the 1x80 arrays concatenated vertically

% Assuming pACC_spikes_CuePeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
pACC_spikes_CuePeriod_excluded_normalized = zeros(size(pACC_spikes_CuePeriod_excluded));

% Normalize each row
for i = 1:size(pACC_spikes_CuePeriod_excluded, 1)
    row = pACC_spikes_CuePeriod_excluded(i, :);
    normalizedRow_excluded = (row - mean(row)) / std(row);
    pACC_spikes_CuePeriod_excluded_normalized(i, :) = normalizedRow_excluded;
end

% pACC_spikes_CuePeriod_normalized now contains the normalized rows

% Assuming pACC_spikes_CuePeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN = all(isnan(pACC_spikes_CuePeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero = all(pACC_spikes_CuePeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove = rowsAllNaN | rowsAllZero;

% Filter out these rows
pACC_spikes_CuePeriod_excluded_normalized_filtered = pACC_spikes_CuePeriod_excluded_normalized(~rowsToRemove, :);


comeon1=[pACC_spikes_CuePeriod_excluded_normalized_filtered;pACC_spikes_CuePeriod_normalized_filtered];

%%%% for pACC airpuff period
% Filter for pACC at airpuff period
pACC_matrix2_air = matrix2(strcmp(matrix2(:, 6), 'pACC'), :);

% Filter rows based on the value in column 64 being 1 or -1 and sort
subset1_pACC_air = pACC_matrix2_air(cell2mat(pACC_matrix2_air(:, 65)) == 1, :);
subset2_pACC_air = pACC_matrix2_air(cell2mat(pACC_matrix2_air(:, 65)) == -1, :);

% Sorting subsets based on column 67
[~, sortIndices1_pACC_air] = sort(cell2mat(subset1_pACC_air(:, 68)), 'descend');
subset1_sorted_pACC_air = subset1_pACC_air(sortIndices1_pACC_air, :);

[~, sortIndices2_pACC_air] = sort(cell2mat(subset2_pACC_air(:, 68)), 'descend');
subset2_sorted_pACC_air = subset2_pACC_air(sortIndices2_pACC_air, :);

% Concatenate the sorted subsets
pACC_matrix2_sorted_air = [subset1_sorted_pACC_air; subset2_sorted_pACC_air];

% Remove rows with NaN in column 68
rowsToKeep_pACC_air = ~cellfun(@(x) isnan(x), pACC_matrix2_sorted_air(:, 68));
pACC_matrix2_noNaN_air = pACC_matrix2_sorted_air(rowsToKeep_pACC_air, :);
pACC_matrix2_noNaN_excluded_air = pACC_matrix2_noNaN_air;

% Logical index for rows to keep in pACC_matrix2_noNaN_air (values in column 68 >= 0.05)
keepRows_pACC_air = cell2mat(pACC_matrix2_noNaN_air(:, 68)) <= 0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_air = pACC_matrix2_noNaN_air(keepRows_pACC_air, :);

% Logical index for rows to keep in pACC_matrix2_noNaN (values in column 68 >= 0.05)
keepRows_pACC_excluded_air = cell2mat(pACC_matrix2_noNaN_excluded_air(:, 68)) >0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_excluded_air = pACC_matrix2_noNaN_excluded_air(keepRows_pACC_excluded_air, :);

% Concatenate data from column 63 into a new matrix
pACC_spikes_AirpuffPeriod = [];
for i = 1:size(pACC_matrix2_noNaN_air, 1)
    currentArray_pACC_air = pACC_matrix2_noNaN_air{i, 88};
    pACC_spikes_AirpuffPeriod = [pACC_spikes_AirpuffPeriod; currentArray_pACC_air];
end

% Normalize each row of pACC_spikes_AirpuffPeriod
pACC_spikes_AirpuffPeriod_normalized = zeros(size(pACC_spikes_AirpuffPeriod));
for i = 1:size(pACC_spikes_AirpuffPeriod, 1)
    row = pACC_spikes_AirpuffPeriod(i, :);
    normalizedRow_air = (row - mean(row)) / std(row);
    pACC_spikes_AirpuffPeriod_normalized(i, :) = normalizedRow_air;
end

% Remove rows that are all NaN or all zero
rowsAllNaN_pACC_air = all(isnan(pACC_spikes_AirpuffPeriod_normalized), 2);
rowsAllZero_pACC_air = all(pACC_spikes_AirpuffPeriod_normalized == 0, 2);
rowsToRemove_pACC_air = rowsAllNaN_pACC_air | rowsAllZero_pACC_air;
pACC_spikes_AirpuffPeriod_normalized_filtered = pACC_spikes_AirpuffPeriod_normalized(~rowsToRemove_pACC_air, :);

%%%% the ones excluded
% Initialize an empty matrix for concatenation
pACC_spikes_AirpuffPeriod_excluded= [];

% Iterate through each row in column 88 of pACC_matrix2_noNaN_air
for i = 1:size(pACC_matrix2_noNaN_excluded_air, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_air  = pACC_matrix2_noNaN_excluded_air{i, 88};
    
    % Concatenate this array to the growing matrix
    pACC_spikes_AirpuffPeriod_excluded = [pACC_spikes_AirpuffPeriod_excluded; currentArray_excluded_air];
end

% pACC_spikes_AirpuffPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming pACC_spikes_AirpuffPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
pACC_spikes_AirpuffPeriod_excluded_normalized = zeros(size(pACC_spikes_AirpuffPeriod_excluded));

% Normalize each row
for i = 1:size(pACC_spikes_AirpuffPeriod_excluded, 1)
    row = pACC_spikes_AirpuffPeriod_excluded(i, :);
    normalizedRow_excluded_air = (row - mean(row)) / std(row);
    pACC_spikes_AirpuffPeriod_excluded_normalized(i, :) = normalizedRow_excluded_air;
end

% pACC_spikes_AirpuffPeriod_normalized now contains the normalized rows

% Assuming pACC_spikes_AirpuffPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_air = all(isnan(pACC_spikes_AirpuffPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_air = all(pACC_spikes_AirpuffPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_air = rowsAllNaN_air | rowsAllZero_air;

% Filter out these rows
pACC_spikes_AirpuffPeriod_excluded_normalized_filtered = pACC_spikes_AirpuffPeriod_excluded_normalized(~rowsToRemove_air, :);

comeon1_air=[pACC_spikes_AirpuffPeriod_excluded_normalized_filtered;pACC_spikes_AirpuffPeriod_normalized_filtered];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  small reward
%%%% for pACC airpuff period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter for pACC at airpuff period
pACC_matrix2_smallreward = matrix2(strcmp(matrix2(:, 6), 'pACC'), :);

% Filter rows based on the value in column 64 being 1 or -1 and sort
subset1_pACC_smallreward = pACC_matrix2_smallreward(cell2mat(pACC_matrix2_smallreward(:, 80)) == 1, :);
subset2_pACC_smallreward = pACC_matrix2_smallreward(cell2mat(pACC_matrix2_smallreward(:, 80)) == -1, :);

% Sorting subsets based on column 69
[~, sortIndices1_pACC_smallreward] = sort(cell2mat(subset1_pACC_smallreward(:, 69)), 'descend');
subset1_sorted_pACC_smallreward = subset1_pACC_smallreward(sortIndices1_pACC_smallreward, :);

[~, sortIndices2_pACC_smallreward] = sort(cell2mat(subset2_pACC_smallreward(:, 69)), 'descend');
subset2_sorted_pACC_smallreward = subset2_pACC_smallreward(sortIndices2_pACC_smallreward, :);

% Concatenate the sorted subsets
pACC_matrix2_sorted_smallreward = [subset1_sorted_pACC_smallreward; subset2_sorted_pACC_smallreward];


% Remove rows with NaN in column 69
rowsToKeep_pACC_smallreward = ~cellfun(@(x) isnan(x), pACC_matrix2_sorted_smallreward(:, 69));
pACC_matrix2_noNaN_smallreward = pACC_matrix2_sorted_smallreward(rowsToKeep_pACC_smallreward, :);
pACC_matrix2_noNaN_excluded_smallreward = pACC_matrix2_noNaN_smallreward;

% Logical index for rows to keep in pACC_matrix2_noNaN_smallreward (values in column 69 >= 0.05)
keepRows_pACC_smallreward = cell2mat(pACC_matrix2_noNaN_smallreward(:, 69)) <= 0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_smallreward = pACC_matrix2_noNaN_smallreward(keepRows_pACC_smallreward, :);

% Logical index for rows to keep in pACC_matrix2_noNaN (values in column 69 >= 0.05)
keepRows_pACC_excluded_smallreward = cell2mat(pACC_matrix2_noNaN_excluded_smallreward(:, 69)) >0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_excluded_smallreward = pACC_matrix2_noNaN_excluded_smallreward(keepRows_pACC_excluded_smallreward, :);

% Concatenate data from column 89 into a new matrix
pACC_spikes_SmallRewardPeriod = [];
for i = 1:size(pACC_matrix2_noNaN_smallreward, 1)
    currentArray_pACC_smallreward = pACC_matrix2_noNaN_smallreward{i, 89};
    pACC_spikes_SmallRewardPeriod = [pACC_spikes_SmallRewardPeriod; currentArray_pACC_smallreward];
end

% Normalize each row of pACC_spikes_SmallRewardPeriod
pACC_spikes_SmallRewardPeriod_normalized = zeros(size(pACC_spikes_SmallRewardPeriod));
for i = 1:size(pACC_spikes_SmallRewardPeriod, 1)
    row = pACC_spikes_SmallRewardPeriod(i, :);
    normalizedRow_smallreward = (row - mean(row)) / std(row);
    pACC_spikes_SmallRewardPeriod_normalized(i, :) = normalizedRow_smallreward;
end

% Remove rows that are all NaN or all zero
rowsAllNaN_pACC_smallreward = all(isnan(pACC_spikes_SmallRewardPeriod_normalized), 2);
rowsAllZero_pACC_smallreward = all(pACC_spikes_SmallRewardPeriod_normalized == 0, 2);
rowsToRemove_pACC_smallreward = rowsAllNaN_pACC_smallreward | rowsAllZero_pACC_smallreward;
pACC_spikes_SmallRewardPeriod_normalized_filtered = pACC_spikes_SmallRewardPeriod_normalized(~rowsToRemove_pACC_smallreward, :);

%%%% the ones excluded
% Initialize an empty matrix for concatenation
pACC_spikes_SmallRewardPeriod_excluded= [];

% Iterate through each row in column 89 of pACC_matrix2_noNaN_smallreward
for i = 1:size(pACC_matrix2_noNaN_excluded_smallreward, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_smallreward  = pACC_matrix2_noNaN_excluded_smallreward{i, 89};
    
    % Concatenate this array to the growing matrix
    pACC_spikes_SmallRewardPeriod_excluded = [pACC_spikes_SmallRewardPeriod_excluded; currentArray_excluded_smallreward];
end

% pACC_spikes_SmallRewardPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming pACC_spikes_SmallRewardPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
pACC_spikes_SmallRewardPeriod_excluded_normalized = zeros(size(pACC_spikes_SmallRewardPeriod_excluded));

% Normalize each row
for i = 1:size(pACC_spikes_SmallRewardPeriod_excluded, 1)
    row = pACC_spikes_SmallRewardPeriod_excluded(i, :);
    normalizedRow_excluded_smallreward = (row - mean(row)) / std(row);
    pACC_spikes_SmallRewardPeriod_excluded_normalized(i, :) = normalizedRow_excluded_smallreward;
end

% pACC_spikes_SmallRewardPeriod_normalized now contains the normalized rows

% Assuming pACC_spikes_SmallRewardPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_smallreward = all(isnan(pACC_spikes_SmallRewardPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_smallreward = all(pACC_spikes_SmallRewardPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_smallreward = rowsAllNaN_smallreward | rowsAllZero_smallreward;

% Filter out these rows
pACC_spikes_SmallRewardPeriod_excluded_normalized_filtered = pACC_spikes_SmallRewardPeriod_excluded_normalized(~rowsToRemove_smallreward, :);


comeon1_smallreward=[pACC_spikes_SmallRewardPeriod_excluded_normalized_filtered;pACC_spikes_SmallRewardPeriod_normalized_filtered];

%%%%%%%%%%%%%%%%%
% Filter for pACC at airpuff period
pACC_matrix2_bigreward = matrix2(strcmp(matrix2(:, 6), 'pACC'), :);

% Filter rows based on the value in column 64 being 1 or -1 and sort
subset1_pACC_bigreward = pACC_matrix2_bigreward(cell2mat(pACC_matrix2_bigreward(:, 87)) == 1, :);
subset2_pACC_bigreward = pACC_matrix2_bigreward(cell2mat(pACC_matrix2_bigreward(:, 87)) == -1, :);

% Sorting subsets based on column 76
[~, sortIndices1_pACC_bigreward] = sort(cell2mat(subset1_pACC_bigreward(:, 76)), 'descend');
subset1_sorted_pACC_bigreward = subset1_pACC_bigreward(sortIndices1_pACC_bigreward, :);

[~, sortIndices2_pACC_bigreward] = sort(cell2mat(subset2_pACC_bigreward(:, 76)), 'descend');
subset2_sorted_pACC_bigreward = subset2_pACC_bigreward(sortIndices2_pACC_bigreward, :);

% Concatenate the sorted subsets
pACC_matrix2_sorted_bigreward = [subset1_sorted_pACC_bigreward; subset2_sorted_pACC_bigreward];


% Remove rows with NaN in column 76
rowsToKeep_pACC_bigreward = ~cellfun(@(x) isnan(x), pACC_matrix2_sorted_bigreward(:, 76));
pACC_matrix2_noNaN_bigreward = pACC_matrix2_sorted_bigreward(rowsToKeep_pACC_bigreward, :);
pACC_matrix2_noNaN_excluded_bigreward = pACC_matrix2_noNaN_bigreward;

% Logical index for rows to keep in pACC_matrix2_noNaN_bigreward (values in column 76 >= 0.05)
keepRows_pACC_bigreward = cell2mat(pACC_matrix2_noNaN_bigreward(:, 76)) <= 0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_bigreward = pACC_matrix2_noNaN_bigreward(keepRows_pACC_bigreward, :);

% Logical index for rows to keep in pACC_matrix2_noNaN (values in column 76 >= 0.05)
keepRows_pACC_excluded_bigreward = cell2mat(pACC_matrix2_noNaN_excluded_bigreward(:, 76)) >0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_excluded_bigreward = pACC_matrix2_noNaN_excluded_bigreward(keepRows_pACC_excluded_bigreward, :);

% Concatenate data from column 90 into a new matrix
pACC_spikes_BigRewardPeriod = [];
for i = 1:size(pACC_matrix2_noNaN_bigreward, 1)
    currentArray_pACC_bigreward = pACC_matrix2_noNaN_bigreward{i, 90};
    pACC_spikes_BigRewardPeriod = [pACC_spikes_BigRewardPeriod; currentArray_pACC_bigreward];
end

% Normalize each row of pACC_spikes_BigRewardPeriod
pACC_spikes_BigRewardPeriod_normalized = zeros(size(pACC_spikes_BigRewardPeriod));
for i = 1:size(pACC_spikes_BigRewardPeriod, 1)
    row = pACC_spikes_BigRewardPeriod(i, :);
    normalizedRow_bigreward = (row - mean(row)) / std(row);
    pACC_spikes_BigRewardPeriod_normalized(i, :) = normalizedRow_bigreward;
end

% Remove rows that are all NaN or all zero
rowsAllNaN_pACC_bigreward = all(isnan(pACC_spikes_BigRewardPeriod_normalized), 2);
rowsAllZero_pACC_bigreward = all(pACC_spikes_BigRewardPeriod_normalized == 0, 2);
rowsToRemove_pACC_bigreward = rowsAllNaN_pACC_bigreward | rowsAllZero_pACC_bigreward;
pACC_spikes_BigRewardPeriod_normalized_filtered = pACC_spikes_BigRewardPeriod_normalized(~rowsToRemove_pACC_bigreward, :);

%%%% the ones excluded
% Initialize an empty matrix for concatenation
pACC_spikes_BigRewardPeriod_excluded= [];

% Iterate through each row in column 90 of pACC_matrix2_noNaN_bigreward
for i = 1:size(pACC_matrix2_noNaN_excluded_bigreward, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_bigreward  = pACC_matrix2_noNaN_excluded_bigreward{i, 90};
    
    % Concatenate this array to the growing matrix
    pACC_spikes_BigRewardPeriod_excluded = [pACC_spikes_BigRewardPeriod_excluded; currentArray_excluded_bigreward];
end

% pACC_spikes_BigRewardPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming pACC_spikes_BigRewardPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
pACC_spikes_BigRewardPeriod_excluded_normalized = zeros(size(pACC_spikes_BigRewardPeriod_excluded));

% Normalize each row
for i = 1:size(pACC_spikes_BigRewardPeriod_excluded, 1)
    row = pACC_spikes_BigRewardPeriod_excluded(i, :);
    normalizedRow_excluded_bigreward = (row - mean(row)) / std(row);
    pACC_spikes_BigRewardPeriod_excluded_normalized(i, :) = normalizedRow_excluded_bigreward;
end

% pACC_spikes_BigRewardPeriod_normalized now contains the normalized rows

% Assuming pACC_spikes_BigRewardPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_bigreward = all(isnan(pACC_spikes_BigRewardPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_bigreward = all(pACC_spikes_BigRewardPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_bigreward = rowsAllNaN_bigreward | rowsAllZero_bigreward;

% Filter out these rows
pACC_spikes_BigRewardPeriod_excluded_normalized_filtered = pACC_spikes_BigRewardPeriod_excluded_normalized(~rowsToRemove_bigreward, :);

comeon1_bigreward=[pACC_spikes_BigRewardPeriod_excluded_normalized_filtered;pACC_spikes_BigRewardPeriod_normalized_filtered];

% Filter for pACC at airpuff period
pACC_matrix2_redcue = matrix2(strcmp(matrix2(:, 6), 'pACC'), :);

% Filter rows based on the value in column 64 being 1 or -1 and sort
subset1_pACC_redcue = pACC_matrix2_redcue(cell2mat(pACC_matrix2_redcue(:, 83)) == 1, :);
subset2_pACC_redcue = pACC_matrix2_redcue(cell2mat(pACC_matrix2_redcue(:, 83)) == -1, :);

% Sorting subsets based on column 72
[~, sortIndices1_pACC_redcue] = sort(cell2mat(subset1_pACC_redcue(:, 72)), 'descend');
subset1_sorted_pACC_redcue = subset1_pACC_redcue(sortIndices1_pACC_redcue, :);

[~, sortIndices2_pACC_redcue] = sort(cell2mat(subset2_pACC_redcue(:, 72)), 'descend');
subset2_sorted_pACC_redcue = subset2_pACC_redcue(sortIndices2_pACC_redcue, :);

% Concatenate the sorted subsets
pACC_matrix2_sorted_redcue = [subset1_sorted_pACC_redcue; subset2_sorted_pACC_redcue];


% Remove rows with NaN in column 72
rowsToKeep_pACC_redcue = ~cellfun(@(x) isnan(x), pACC_matrix2_sorted_redcue(:, 72));
pACC_matrix2_noNaN_redcue = pACC_matrix2_sorted_redcue(rowsToKeep_pACC_redcue, :);
pACC_matrix2_noNaN_excluded_redcue = pACC_matrix2_noNaN_redcue;

% Logical index for rows to keep in pACC_matrix2_noNaN_redcue (values in column 72 >= 0.05)
keepRows_pACC_redcue = cell2mat(pACC_matrix2_noNaN_redcue(:, 72)) <= 0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_redcue = pACC_matrix2_noNaN_redcue(keepRows_pACC_redcue, :);

% Logical index for rows to keep in pACC_matrix2_noNaN (values in column 72 >= 0.05)
keepRows_pACC_excluded_redcue = cell2mat(pACC_matrix2_noNaN_excluded_redcue(:, 72)) >0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_excluded_redcue = pACC_matrix2_noNaN_excluded_redcue(keepRows_pACC_excluded_redcue, :);

% Concatenate data from column 91 into a new matrix
pACC_spikes_RedCuePeriod = [];
for i = 1:size(pACC_matrix2_noNaN_redcue, 1)
    currentArray_pACC_redcue = pACC_matrix2_noNaN_redcue{i, 91};
    pACC_spikes_RedCuePeriod = [pACC_spikes_RedCuePeriod; currentArray_pACC_redcue];
end

% Normalize each row of pACC_spikes_RedCuePeriod
pACC_spikes_RedCuePeriod_normalized = zeros(size(pACC_spikes_RedCuePeriod));
for i = 1:size(pACC_spikes_RedCuePeriod, 1)
    row = pACC_spikes_RedCuePeriod(i, :);
    normalizedRow_redcue = (row - mean(row)) / std(row);
    pACC_spikes_RedCuePeriod_normalized(i, :) = normalizedRow_redcue;
end

% Remove rows that are all NaN or all zero
rowsAllNaN_pACC_redcue = all(isnan(pACC_spikes_RedCuePeriod_normalized), 2);
rowsAllZero_pACC_redcue = all(pACC_spikes_RedCuePeriod_normalized == 0, 2);
rowsToRemove_pACC_redcue = rowsAllNaN_pACC_redcue | rowsAllZero_pACC_redcue;
pACC_spikes_RedCuePeriod_normalized_filtered = pACC_spikes_RedCuePeriod_normalized(~rowsToRemove_pACC_redcue, :);

%%%% the ones excluded
% Initialize an empty matrix for concatenation
pACC_spikes_RedCuePeriod_excluded= [];

% Iterate through each row in column 91 of pACC_matrix2_noNaN_redcue
for i = 1:size(pACC_matrix2_noNaN_excluded_redcue, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_redcue  = pACC_matrix2_noNaN_excluded_redcue{i, 91};
    
    % Concatenate this array to the growing matrix
    pACC_spikes_RedCuePeriod_excluded = [pACC_spikes_RedCuePeriod_excluded; currentArray_excluded_redcue];
end

% pACC_spikes_RedCuePeriod now contains all the 1x80 arrays concatenated vertically

% Assuming pACC_spikes_RedCuePeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
pACC_spikes_RedCuePeriod_excluded_normalized = zeros(size(pACC_spikes_RedCuePeriod_excluded));

% Normalize each row
for i = 1:size(pACC_spikes_RedCuePeriod_excluded, 1)
    row = pACC_spikes_RedCuePeriod_excluded(i, :);
    normalizedRow_excluded_redcue = (row - mean(row)) / std(row);
    pACC_spikes_RedCuePeriod_excluded_normalized(i, :) = normalizedRow_excluded_redcue;
end

% pACC_spikes_RedCuePeriod_normalized now contains the normalized rows

% Assuming pACC_spikes_RedCuePeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_redcue = all(isnan(pACC_spikes_RedCuePeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_redcue = all(pACC_spikes_RedCuePeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_redcue = rowsAllNaN_redcue | rowsAllZero_redcue;

% Filter out these rows
pACC_spikes_RedCuePeriod_excluded_normalized_filtered = pACC_spikes_RedCuePeriod_excluded_normalized(~rowsToRemove_redcue, :);


comeon1_redcue=[pACC_spikes_RedCuePeriod_excluded_normalized_filtered;pACC_spikes_RedCuePeriod_normalized_filtered];

% Filter for pACC at airpuff period
pACC_matrix2_yellowcue = matrix2(strcmp(matrix2(:, 6), 'pACC'), :);

% Filter rows based on the value in column 64 being 1 or -1 and sort
subset1_pACC_yellowcue = pACC_matrix2_yellowcue(cell2mat(pACC_matrix2_yellowcue(:, 84)) == 1, :);
subset2_pACC_yellowcue = pACC_matrix2_yellowcue(cell2mat(pACC_matrix2_yellowcue(:, 84)) == -1, :);

% Sorting subsets based on column 73
[~, sortIndices1_pACC_yellowcue] = sort(cell2mat(subset1_pACC_yellowcue(:, 73)), 'descend');
subset1_sorted_pACC_yellowcue = subset1_pACC_yellowcue(sortIndices1_pACC_yellowcue, :);

[~, sortIndices2_pACC_yellowcue] = sort(cell2mat(subset2_pACC_yellowcue(:, 73)), 'descend');
subset2_sorted_pACC_yellowcue = subset2_pACC_yellowcue(sortIndices2_pACC_yellowcue, :);

% Concatenate the sorted subsets
pACC_matrix2_sorted_yellowcue = [subset1_sorted_pACC_yellowcue; subset2_sorted_pACC_yellowcue];


% Remove rows with NaN in column 73
rowsToKeep_pACC_yellowcue = ~cellfun(@(x) isnan(x), pACC_matrix2_sorted_yellowcue(:, 73));
pACC_matrix2_noNaN_yellowcue = pACC_matrix2_sorted_yellowcue(rowsToKeep_pACC_yellowcue, :);
pACC_matrix2_noNaN_excluded_yellowcue = pACC_matrix2_noNaN_yellowcue;




% Logical index for rows to keep in pACC_matrix2_noNaN_yellowcue (values in column 73 >= 0.05)
keepRows_pACC_yellowcue = cell2mat(pACC_matrix2_noNaN_yellowcue(:, 73)) <= 0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_yellowcue = pACC_matrix2_noNaN_yellowcue(keepRows_pACC_yellowcue, :);



% Logical index for rows to keep in pACC_matrix2_noNaN (values in column 73 >= 0.05)
keepRows_pACC_excluded_yellowcue = cell2mat(pACC_matrix2_noNaN_excluded_yellowcue(:, 73)) >0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_excluded_yellowcue = pACC_matrix2_noNaN_excluded_yellowcue(keepRows_pACC_excluded_yellowcue, :);

% Concatenate data from column 92 into a new matrix
pACC_spikes_YellowCuePeriod = [];
for i = 1:size(pACC_matrix2_noNaN_yellowcue, 1)
    currentArray_pACC_yellowcue = pACC_matrix2_noNaN_yellowcue{i, 92};
    pACC_spikes_YellowCuePeriod = [pACC_spikes_YellowCuePeriod; currentArray_pACC_yellowcue];
end

% Normalize each row of pACC_spikes_YellowCuePeriod
pACC_spikes_YellowCuePeriod_normalized = zeros(size(pACC_spikes_YellowCuePeriod));
for i = 1:size(pACC_spikes_YellowCuePeriod, 1)
    row = pACC_spikes_YellowCuePeriod(i, :);
    normalizedRow_yellowcue = (row - mean(row)) / std(row);
    pACC_spikes_YellowCuePeriod_normalized(i, :) = normalizedRow_yellowcue;
end

% Remove rows that are all NaN or all zero
rowsAllNaN_pACC_yellowcue = all(isnan(pACC_spikes_YellowCuePeriod_normalized), 2);
rowsAllZero_pACC_yellowcue = all(pACC_spikes_YellowCuePeriod_normalized == 0, 2);
rowsToRemove_pACC_yellowcue = rowsAllNaN_pACC_yellowcue | rowsAllZero_pACC_yellowcue;
pACC_spikes_YellowCuePeriod_normalized_filtered = pACC_spikes_YellowCuePeriod_normalized(~rowsToRemove_pACC_yellowcue, :);

%%%% the ones excluded
% Initialize an empty matrix for concatenation
pACC_spikes_YellowCuePeriod_excluded= [];

% Iterate through each row in column 92 of pACC_matrix2_noNaN_yellowcue
for i = 1:size(pACC_matrix2_noNaN_excluded_yellowcue, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_yellowcue  = pACC_matrix2_noNaN_excluded_yellowcue{i, 92};
    
    % Concatenate this array to the growing matrix
    pACC_spikes_YellowCuePeriod_excluded = [pACC_spikes_YellowCuePeriod_excluded; currentArray_excluded_yellowcue];
end

% pACC_spikes_YellowCuePeriod now contains all the 1x80 arrays concatenated vertically

% Assuming pACC_spikes_YellowCuePeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
pACC_spikes_YellowCuePeriod_excluded_normalized = zeros(size(pACC_spikes_YellowCuePeriod_excluded));

% Normalize each row
for i = 1:size(pACC_spikes_YellowCuePeriod_excluded, 1)
    row = pACC_spikes_YellowCuePeriod_excluded(i, :);
    normalizedRow_excluded_yellowcue = (row - mean(row)) / std(row);
    pACC_spikes_YellowCuePeriod_excluded_normalized(i, :) = normalizedRow_excluded_yellowcue;
end

% pACC_spikes_YellowCuePeriod_normalized now contains the normalized rows

% Assuming pACC_spikes_YellowCuePeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_yellowcue = all(isnan(pACC_spikes_YellowCuePeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_yellowcue = all(pACC_spikes_YellowCuePeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_yellowcue = rowsAllNaN_yellowcue | rowsAllZero_yellowcue;

% Filter out these rows
pACC_spikes_YellowCuePeriod_excluded_normalized_filtered = pACC_spikes_YellowCuePeriod_excluded_normalized(~rowsToRemove_yellowcue, :);


comeon1_yellowcue=[pACC_spikes_YellowCuePeriod_excluded_normalized_filtered;pACC_spikes_YellowCuePeriod_normalized_filtered];

% Filter for pACC at airpuff period
pACC_matrix2_redreward = matrix2(strcmp(matrix2(:, 6), 'pACC'), :);

% Filter rows based on the value in column 64 being 1 or -1 and sort
subset1_pACC_redreward = pACC_matrix2_redreward(cell2mat(pACC_matrix2_redreward(:, 85)) == 1, :);
subset2_pACC_redreward = pACC_matrix2_redreward(cell2mat(pACC_matrix2_redreward(:, 85)) == -1, :);

% Sorting subsets based on column 74
[~, sortIndices1_pACC_redreward] = sort(cell2mat(subset1_pACC_redreward(:, 74)), 'descend');
subset1_sorted_pACC_redreward = subset1_pACC_redreward(sortIndices1_pACC_redreward, :);

[~, sortIndices2_pACC_redreward] = sort(cell2mat(subset2_pACC_redreward(:, 74)), 'descend');
subset2_sorted_pACC_redreward = subset2_pACC_redreward(sortIndices2_pACC_redreward, :);

% Concatenate the sorted subsets
pACC_matrix2_sorted_redreward = [subset1_sorted_pACC_redreward; subset2_sorted_pACC_redreward];


% Remove rows with NaN in column 74
rowsToKeep_pACC_redreward = ~cellfun(@(x) isnan(x), pACC_matrix2_sorted_redreward(:, 74));
pACC_matrix2_noNaN_redreward = pACC_matrix2_sorted_redreward(rowsToKeep_pACC_redreward, :);
pACC_matrix2_noNaN_excluded_redreward = pACC_matrix2_noNaN_redreward;




% Logical index for rows to keep in pACC_matrix2_noNaN_redreward (values in column 74 >= 0.05)
keepRows_pACC_redreward = cell2mat(pACC_matrix2_noNaN_redreward(:, 74)) <= 0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_redreward = pACC_matrix2_noNaN_redreward(keepRows_pACC_redreward, :);



% Logical index for rows to keep in pACC_matrix2_noNaN (values in column 74 >= 0.05)
keepRows_pACC_excluded_redreward = cell2mat(pACC_matrix2_noNaN_excluded_redreward(:, 74)) >0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_excluded_redreward = pACC_matrix2_noNaN_excluded_redreward(keepRows_pACC_excluded_redreward, :);


% Concatenate data from column 93 into a new matrix
pACC_spikes_RedRewardPeriod = [];
for i = 1:size(pACC_matrix2_noNaN_redreward, 1)
    currentArray_pACC_redreward = pACC_matrix2_noNaN_redreward{i, 93};
    pACC_spikes_RedRewardPeriod = [pACC_spikes_RedRewardPeriod; currentArray_pACC_redreward];
end

% Normalize each row of pACC_spikes_RedRewardPeriod
pACC_spikes_RedRewardPeriod_normalized = zeros(size(pACC_spikes_RedRewardPeriod));
for i = 1:size(pACC_spikes_RedRewardPeriod, 1)
    row = pACC_spikes_RedRewardPeriod(i, :);
    normalizedRow_redreward = (row - mean(row)) / std(row);
    pACC_spikes_RedRewardPeriod_normalized(i, :) = normalizedRow_redreward;
end

% Remove rows that are all NaN or all zero
rowsAllNaN_pACC_redreward = all(isnan(pACC_spikes_RedRewardPeriod_normalized), 2);
rowsAllZero_pACC_redreward = all(pACC_spikes_RedRewardPeriod_normalized == 0, 2);
rowsToRemove_pACC_redreward = rowsAllNaN_pACC_redreward | rowsAllZero_pACC_redreward;
pACC_spikes_RedRewardPeriod_normalized_filtered = pACC_spikes_RedRewardPeriod_normalized(~rowsToRemove_pACC_redreward, :);

%%%% the ones excluded
% Initialize an empty matrix for concatenation
pACC_spikes_RedRewardPeriod_excluded= [];

% Iterate through each row in column 93 of pACC_matrix2_noNaN_redreward
for i = 1:size(pACC_matrix2_noNaN_excluded_redreward, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_redreward  = pACC_matrix2_noNaN_excluded_redreward{i, 93};
    
    % Concatenate this array to the growing matrix
    pACC_spikes_RedRewardPeriod_excluded = [pACC_spikes_RedRewardPeriod_excluded; currentArray_excluded_redreward];
end

% pACC_spikes_RedRewardPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming pACC_spikes_RedRewardPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
pACC_spikes_RedRewardPeriod_excluded_normalized = zeros(size(pACC_spikes_RedRewardPeriod_excluded));

% Normalize each row
for i = 1:size(pACC_spikes_RedRewardPeriod_excluded, 1)
    row = pACC_spikes_RedRewardPeriod_excluded(i, :);
    normalizedRow_excluded_redreward = (row - mean(row)) / std(row);
    pACC_spikes_RedRewardPeriod_excluded_normalized(i, :) = normalizedRow_excluded_redreward;
end

% pACC_spikes_RedRewardPeriod_normalized now contains the normalized rows

% Assuming pACC_spikes_RedRewardPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_redreward = all(isnan(pACC_spikes_RedRewardPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_redreward = all(pACC_spikes_RedRewardPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_redreward = rowsAllNaN_redreward | rowsAllZero_redreward;

% Filter out these rows
pACC_spikes_RedRewardPeriod_excluded_normalized_filtered = pACC_spikes_RedRewardPeriod_excluded_normalized(~rowsToRemove_redreward, :);


comeon1_redreward=[pACC_spikes_RedRewardPeriod_excluded_normalized_filtered;pACC_spikes_RedRewardPeriod_normalized_filtered];

% Filter for pACC at airpuff period
pACC_matrix2_yellowairpuff = matrix2(strcmp(matrix2(:, 6), 'pACC'), :);

% Filter rows based on the value in column 64 being 1 or -1 and sort
subset1_pACC_yellowairpuff = pACC_matrix2_yellowairpuff(cell2mat(pACC_matrix2_yellowairpuff(:, 86)) == 1, :);
subset2_pACC_yellowairpuff = pACC_matrix2_yellowairpuff(cell2mat(pACC_matrix2_yellowairpuff(:, 86)) == -1, :);

% Sorting subsets based on column 75
[~, sortIndices1_pACC_yellowairpuff] = sort(cell2mat(subset1_pACC_yellowairpuff(:, 75)), 'descend');
subset1_sorted_pACC_yellowairpuff = subset1_pACC_yellowairpuff(sortIndices1_pACC_yellowairpuff, :);

[~, sortIndices2_pACC_yellowairpuff] = sort(cell2mat(subset2_pACC_yellowairpuff(:, 75)), 'descend');
subset2_sorted_pACC_yellowairpuff = subset2_pACC_yellowairpuff(sortIndices2_pACC_yellowairpuff, :);

% Concatenate the sorted subsets
pACC_matrix2_sorted_yellowairpuff = [subset1_sorted_pACC_yellowairpuff; subset2_sorted_pACC_yellowairpuff];

% Remove rows with NaN in column 75
rowsToKeep_pACC_yellowairpuff = ~cellfun(@(x) isnan(x), pACC_matrix2_sorted_yellowairpuff(:, 75));
pACC_matrix2_noNaN_yellowairpuff = pACC_matrix2_sorted_yellowairpuff(rowsToKeep_pACC_yellowairpuff, :);
pACC_matrix2_noNaN_excluded_yellowairpuff = pACC_matrix2_noNaN_yellowairpuff;

% Logical index for rows to keep in pACC_matrix2_noNaN_yellowairpuff (values in column 75 >= 0.05)
keepRows_pACC_yellowairpuff = cell2mat(pACC_matrix2_noNaN_yellowairpuff(:, 75)) <= 0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_yellowairpuff = pACC_matrix2_noNaN_yellowairpuff(keepRows_pACC_yellowairpuff, :);

% Logical index for rows to keep in pACC_matrix2_noNaN (values in column 75 >= 0.05)
keepRows_pACC_excluded_yellowairpuff = cell2mat(pACC_matrix2_noNaN_excluded_yellowairpuff(:, 75)) >0.05;

% Use the logical index to filter out rows
pACC_matrix2_noNaN_excluded_yellowairpuff = pACC_matrix2_noNaN_excluded_yellowairpuff(keepRows_pACC_excluded_yellowairpuff, :);

% Concatenate data from column 94 into a new matrix
pACC_spikes_YellowAirpuffPeriod = [];
for i = 1:size(pACC_matrix2_noNaN_yellowairpuff, 1)
    currentArray_pACC_yellowairpuff = pACC_matrix2_noNaN_yellowairpuff{i, 94};
    pACC_spikes_YellowAirpuffPeriod = [pACC_spikes_YellowAirpuffPeriod; currentArray_pACC_yellowairpuff];
end

% Normalize each row of pACC_spikes_YellowAirpuffPeriod
pACC_spikes_YellowAirpuffPeriod_normalized = zeros(size(pACC_spikes_YellowAirpuffPeriod));
for i = 1:size(pACC_spikes_YellowAirpuffPeriod, 1)
    row = pACC_spikes_YellowAirpuffPeriod(i, :);
    normalizedRow_yellowairpuff = (row - mean(row)) / std(row);
    pACC_spikes_YellowAirpuffPeriod_normalized(i, :) = normalizedRow_yellowairpuff;
end

% Remove rows that are all NaN or all zero
rowsAllNaN_pACC_yellowairpuff = all(isnan(pACC_spikes_YellowAirpuffPeriod_normalized), 2);
rowsAllZero_pACC_yellowairpuff = all(pACC_spikes_YellowAirpuffPeriod_normalized == 0, 2);
rowsToRemove_pACC_yellowairpuff = rowsAllNaN_pACC_yellowairpuff | rowsAllZero_pACC_yellowairpuff;
pACC_spikes_YellowAirpuffPeriod_normalized_filtered = pACC_spikes_YellowAirpuffPeriod_normalized(~rowsToRemove_pACC_yellowairpuff, :);

%%%% the ones excluded
% Initialize an empty matrix for concatenation
pACC_spikes_YellowAirpuffPeriod_excluded= [];

% Iterate through each row in column 94 of pACC_matrix2_noNaN_yellowairpuff
for i = 1:size(pACC_matrix2_noNaN_excluded_yellowairpuff, 1)
    % Extract the 1x80 double array from the current row's cell
    currentArray_excluded_yellowairpuff  = pACC_matrix2_noNaN_excluded_yellowairpuff{i, 94};
    
    % Concatenate this array to the growing matrix
    pACC_spikes_YellowAirpuffPeriod_excluded = [pACC_spikes_YellowAirpuffPeriod_excluded; currentArray_excluded_yellowairpuff];
end

% pACC_spikes_YellowAirpuffPeriod now contains all the 1x80 arrays concatenated vertically

% Assuming pACC_spikes_YellowAirpuffPeriod is the matrix where each row needs to be normalized

% Initialize a matrix to hold the normalized data
pACC_spikes_YellowAirpuffPeriod_excluded_normalized = zeros(size(pACC_spikes_YellowAirpuffPeriod_excluded));

% Normalize each row
for i = 1:size(pACC_spikes_YellowAirpuffPeriod_excluded, 1)
    row = pACC_spikes_YellowAirpuffPeriod_excluded(i, :);
    normalizedRow_excluded_yellowairpuff = (row - mean(row)) / std(row);
    pACC_spikes_YellowAirpuffPeriod_excluded_normalized(i, :) = normalizedRow_excluded_yellowairpuff;
end

% pACC_spikes_YellowAirpuffPeriod_normalized now contains the normalized rows
% Assuming pACC_spikes_YellowAirpuffPeriod_normalized is the normalized matrix

% Find rows that are all NaN
rowsAllNaN_yellowairpuff = all(isnan(pACC_spikes_YellowAirpuffPeriod_excluded_normalized), 2);

% Find rows that are all zeros
rowsAllZero_yellowairpuff = all(pACC_spikes_YellowAirpuffPeriod_excluded_normalized == 0, 2);

% Combine the conditions to identify rows to remove (either all NaN or all zero)
rowsToRemove_yellowairpuff = rowsAllNaN_yellowairpuff | rowsAllZero_yellowairpuff;

% Filter out these rows
pACC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered = pACC_spikes_YellowAirpuffPeriod_excluded_normalized(~rowsToRemove_yellowairpuff, :);

sigma=1;
comeon1_yellowairpuff=[pACC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered;pACC_spikes_YellowAirpuffPeriod_normalized_filtered];
comeon1_yellowairpuff_smoothed = imgaussfilt(comeon1_yellowairpuff, sigma);

% Extract the values from columns 64 and 67
valuesCol64 = cell2mat(pACC_matrix2_sorted(:, 64));
valuesCol67 = cell2mat(pACC_matrix2_sorted(:, 67));

% Define conditions for excitation, inhibition, and unresponsive pACC
excitation_pACC_rows = find(valuesCol64 == 1 & valuesCol67 <= 0.05);
inhibition_pACC_rows = find(valuesCol64 == -1 & valuesCol67 <= 0.05);
condition1 = (valuesCol64 == 1 & valuesCol67 <= 0.05);
condition2 = (valuesCol64 == -1 & valuesCol67 <= 0.05);

% Combine conditions to identify unresponsive rows
unresponsive_pACC_rows = ~(condition1 | condition2);
totalUnresponsive_pACC = sum(unresponsive_pACC_rows);

% Total rows in pACC_matrix2_sorted
totalRows = size(pACC_matrix2_sorted, 1);

% Calculate percentages for excitation, inhibition, and unresponsive
percentageExcitation_pACC = (length(excitation_pACC_rows) / totalRows) * 100;
percentageInhibition_pACC = (length(inhibition_pACC_rows) / totalRows) * 100;
percentageUnresponsive_pACC = (totalUnresponsive_pACC / totalRows) * 100;

% Display the counts and percentages
disp(['Number of excitation_pACC rows: ', num2str(length(excitation_pACC_rows))]);
disp(['Number of inhibition_pACC rows: ', num2str(length(inhibition_pACC_rows))]);
disp(['Number of unresponsive_pACC rows: ', num2str(totalUnresponsive_pACC)]);

disp(['Percentage of excitation_pACC rows: ', num2str(percentageExcitation_pACC), '%']);
disp(['Percentage of inhibition_pACC rows: ', num2str(percentageInhibition_pACC), '%']);
disp(['Percentage of unresponsive_pACC rows: ', num2str(percentageUnresponsive_pACC), '%']);

%%%%%%%%%%%%% for airpuff
% Extract the values from columns 64 and 67
valuesCol65_air = cell2mat(pACC_matrix2_sorted_air(:, 65));
valuesCol68_air = cell2mat(pACC_matrix2_sorted_air(:, 68));

% Define conditions for excitation, inhibition, and unresponsive pACC
excitation_pACC_rows_air = find(valuesCol65_air == 1 & valuesCol68_air <= 0.05);
inhibition_pACC_rows_air = find(valuesCol65_air == -1 & valuesCol68_air <= 0.05);
condition1_air = (valuesCol65_air == 1 & valuesCol68_air <= 0.05);
condition2_air = (valuesCol65_air == -1 & valuesCol68_air <= 0.05);

% Combine conditions to identify unresponsive rows
unresponsive_pACC_rows_air = ~(condition1_air | condition2_air);
totalUnresponsive_pACC_air = sum(unresponsive_pACC_rows_air);

% Total rows in pACC_matrix2_sorted
totalRows_air = size(pACC_matrix2_sorted_air, 1);

% Calculate percentages for excitation, inhibition, and unresponsive
percentageExcitation_pACC_air = (length(excitation_pACC_rows_air) / totalRows_air) * 100;
percentageInhibition_pACC_air = (length(inhibition_pACC_rows_air) / totalRows_air) * 100;
percentageUnresponsive_pACC_air = (totalUnresponsive_pACC_air / totalRows_air) * 100;

% Display the counts and percentages
disp(['Number of excitation_pACC rows at Approach Outcome: ', num2str(length(excitation_pACC_rows_air))]);
disp(['Number of inhibition_pACC rows at Approach Outcome: ', num2str(length(inhibition_pACC_rows_air))]);
disp(['Number of unresponsive_pACC rows at Approach Outcome: ', num2str(totalUnresponsive_pACC_air)]);

disp(['Percentage of excitation_pACC rows at Approach Outcome: ', num2str(percentageExcitation_pACC_air), '%']);
disp(['Percentage of inhibition_pACC rows at Approach Outcome: ', num2str(percentageInhibition_pACC_air), '%']);
disp(['Percentage of unresponsive_pACC rows at Approach Outcome: ', num2str(percentageUnresponsive_pACC_air), '%']);

%%%%%%%%%%%%% for small reward
% Extract the values from columns 80 and 69
valuesCol80_smallreward = cell2mat(pACC_matrix2_sorted_smallreward(:, 80));
valuesCol69_smallreward = cell2mat(pACC_matrix2_sorted_smallreward(:, 69));

% Define conditions for excitation, inhibition, and unresponsive pACC
excitation_pACC_rows_smallreward = find(valuesCol80_smallreward == 1 & valuesCol69_smallreward <= 0.05);
inhibition_pACC_rows_smallreward = find(valuesCol80_smallreward == -1 & valuesCol69_smallreward <= 0.05);
condition1_smallreward = (valuesCol80_smallreward == 1 & valuesCol69_smallreward <= 0.05);
condition2_smallreward = (valuesCol80_smallreward == -1 & valuesCol69_smallreward <= 0.05);

% Combine conditions to identify unresponsive rows
unresponsive_pACC_rows_smallreward = ~(condition1_smallreward | condition2_smallreward);
totalUnresponsive_pACC_smallreward = sum(unresponsive_pACC_rows_smallreward);

% Total rows in pACC_matrix2_sorted
totalRows_smallreward = size(pACC_matrix2_sorted_smallreward, 1);

% Calculate percentages for excitation, inhibition, and unresponsive
percentageExcitation_pACC_smallreward = (length(excitation_pACC_rows_smallreward) / totalRows_smallreward) * 100;
percentageInhibition_pACC_smallreward = (length(inhibition_pACC_rows_smallreward) / totalRows_smallreward) * 100;
percentageUnresponsive_pACC_smallreward = (totalUnresponsive_pACC_smallreward / totalRows_smallreward) * 100;

% Display the counts and percentages
disp(['Number of excitation_pACC rows at Avoidance Outcome: ', num2str(length(excitation_pACC_rows_smallreward))]);
disp(['Number of inhibition_pACC rows at Avoidance Outcome: ', num2str(length(inhibition_pACC_rows_smallreward))]);
disp(['Number of unresponsive_pACC rows at Avoidance Outcome: ', num2str(totalUnresponsive_pACC_smallreward)]);

disp(['Percentage of excitation_pACC rows at Avoidance Outcome: ', num2str(percentageExcitation_pACC_smallreward), '%']);
disp(['Percentage of inhibition_pACC rows at Avoidance Outcome: ', num2str(percentageInhibition_pACC_smallreward), '%']);
disp(['Percentage of unresponsive_pACC rows at Avoidance Outcome: ', num2str(percentageUnresponsive_pACC_smallreward), '%']);

%%%%%%%%%%%%% for small reward
% Extract the values from columns 87 and 76
valuesCol87_bigreward = cell2mat(pACC_matrix2_sorted_bigreward(:, 87));
valuesCol76_bigreward = cell2mat(pACC_matrix2_sorted_bigreward(:, 76));

% Define conditions for excitation, inhibition, and unresponsive pACC
excitation_pACC_rows_bigreward = find(valuesCol87_bigreward == 1 & valuesCol76_bigreward <= 0.05);
inhibition_pACC_rows_bigreward = find(valuesCol87_bigreward == -1 & valuesCol76_bigreward <= 0.05);
condition1_bigreward = (valuesCol87_bigreward == 1 & valuesCol76_bigreward <= 0.05);
condition2_bigreward = (valuesCol87_bigreward == -1 & valuesCol76_bigreward <= 0.05);

% Combine conditions to identify unresponsive rows
unresponsive_pACC_rows_bigreward = ~(condition1_bigreward | condition2_bigreward);
totalUnresponsive_pACC_bigreward = sum(unresponsive_pACC_rows_bigreward);

% Total rows in pACC_matrix2_sorted
totalRows_bigreward = size(pACC_matrix2_sorted_bigreward, 1);

% Calculate percentages for excitation, inhibition, and unresponsive
percentageExcitation_pACC_bigreward = (length(excitation_pACC_rows_bigreward) / totalRows_bigreward) * 100;
percentageInhibition_pACC_bigreward = (length(inhibition_pACC_rows_bigreward) / totalRows_bigreward) * 100;
percentageUnresponsive_pACC_bigreward = (totalUnresponsive_pACC_bigreward / totalRows_bigreward) * 100;

% Display the counts and percentages
disp(['Number of excitation_pACC rows at Approach Reward Outcome: ', num2str(length(excitation_pACC_rows_bigreward))]);
disp(['Number of inhibition_pACC rows at Approach Reward Outcome: ', num2str(length(inhibition_pACC_rows_bigreward))]);
disp(['Number of unresponsive_pACC rows at Approach Reward Outcome: ', num2str(totalUnresponsive_pACC_bigreward)]);

disp(['Percentage of excitation_pACC rows at Approach Reward Outcome: ', num2str(percentageExcitation_pACC_bigreward), '%']);
disp(['Percentage of inhibition_pACC rows at Approach Reward Outcome: ', num2str(percentageInhibition_pACC_bigreward), '%']);
disp(['Percentage of unresponsive_pACC rows at Approach Reward Outcome: ', num2str(percentageUnresponsive_pACC_bigreward), '%']);

% Extract the values from columns 83 and 72
valuesCol83_redcue = cell2mat(pACC_matrix2_sorted_redcue(:, 83));
valuesCol72_redcue = cell2mat(pACC_matrix2_sorted_redcue(:, 72));

% Define conditions for excitation, inhibition, and unresponsive pACC
excitation_pACC_rows_redcue = find(valuesCol83_redcue == 1 & valuesCol72_redcue <= 0.05);
inhibition_pACC_rows_redcue = find(valuesCol83_redcue == -1 & valuesCol72_redcue <= 0.05);
condition1_redcue = (valuesCol83_redcue == 1 & valuesCol72_redcue <= 0.05);
condition2_redcue = (valuesCol83_redcue == -1 & valuesCol72_redcue <= 0.05);

% Combine conditions to identify unresponsive rows
unresponsive_pACC_rows_redcue = ~(condition1_redcue | condition2_redcue);
totalUnresponsive_pACC_redcue = sum(unresponsive_pACC_rows_redcue);

% Total rows in pACC_matrix2_sorted
totalRows_redcue = size(pACC_matrix2_sorted_redcue, 1);

% Calculate percentages for excitation, inhibition, and unresponsive
percentageExcitation_pACC_redcue = (length(excitation_pACC_rows_redcue) / totalRows_redcue) * 100;
percentageInhibition_pACC_redcue = (length(inhibition_pACC_rows_redcue) / totalRows_redcue) * 100;
percentageUnresponsive_pACC_redcue = (totalUnresponsive_pACC_redcue / totalRows_redcue) * 100;

% Display the counts and percentages
disp(['Number of excitation_pACC rows at Red Cue: ', num2str(length(excitation_pACC_rows_redcue))]);
disp(['Number of inhibition_pACC rows at Red Cue: ', num2str(length(inhibition_pACC_rows_redcue))]);
disp(['Number of unresponsive_pACC rows at Red Cue: ', num2str(totalUnresponsive_pACC_redcue)]);

disp(['Percentage of excitation_pACC rows at Red Cue: ', num2str(percentageExcitation_pACC_redcue), '%']);
disp(['Percentage of inhibition_pACC rows at Red Cue: ', num2str(percentageInhibition_pACC_redcue), '%']);
disp(['Percentage of unresponsive_pACC rows at Red Cue: ', num2str(percentageUnresponsive_pACC_redcue), '%']);

% Extract the values from columns 84 and 73
valuesCol84_yellowcue = cell2mat(pACC_matrix2_sorted_yellowcue(:, 84));
valuesCol73_yellowcue = cell2mat(pACC_matrix2_sorted_yellowcue(:, 73));

% Define conditions for excitation, inhibition, and unresponsive pACC
excitation_pACC_rows_yellowcue = find(valuesCol84_yellowcue == 1 & valuesCol73_yellowcue <= 0.05);
inhibition_pACC_rows_yellowcue = find(valuesCol84_yellowcue == -1 & valuesCol73_yellowcue <= 0.05);
condition1_yellowcue = (valuesCol84_yellowcue == 1 & valuesCol73_yellowcue <= 0.05);
condition2_yellowcue = (valuesCol84_yellowcue == -1 & valuesCol73_yellowcue <= 0.05);

% Combine conditions to identify unresponsive rows
unresponsive_pACC_rows_yellowcue = ~(condition1_yellowcue | condition2_yellowcue);
totalUnresponsive_pACC_yellowcue = sum(unresponsive_pACC_rows_yellowcue);

% Total rows in pACC_matrix2_sorted
totalRows_yellowcue = size(pACC_matrix2_sorted_yellowcue, 1);

% Calculate percentages for excitation, inhibition, and unresponsive
percentageExcitation_pACC_yellowcue = (length(excitation_pACC_rows_yellowcue) / totalRows_yellowcue) * 100;
percentageInhibition_pACC_yellowcue = (length(inhibition_pACC_rows_yellowcue) / totalRows_yellowcue) * 100;
percentageUnresponsive_pACC_yellowcue = (totalUnresponsive_pACC_yellowcue / totalRows_yellowcue) * 100;

% Display the counts and percentages
disp(['Number of excitation_pACC rows at Yellow Cue: ', num2str(length(excitation_pACC_rows_yellowcue))]);
disp(['Number of inhibition_pACC rows at Yellow Cue: ', num2str(length(inhibition_pACC_rows_yellowcue))]);
disp(['Number of unresponsive_pACC rows at Yellow Cue: ', num2str(totalUnresponsive_pACC_yellowcue)]);

disp(['Percentage of excitation_pACC rows at Yellow Cue: ', num2str(percentageExcitation_pACC_yellowcue), '%']);
disp(['Percentage of inhibition_pACC rows at Yellow Cue: ', num2str(percentageInhibition_pACC_yellowcue), '%']);
disp(['Percentage of unresponsive_pACC rows at Yellow Cue: ', num2str(percentageUnresponsive_pACC_yellowcue), '%']);

% Extract the values from columns 85 and 74
valuesCol85_redreward = cell2mat(pACC_matrix2_sorted_redreward(:, 85));
valuesCol74_redreward = cell2mat(pACC_matrix2_sorted_redreward(:, 74));

% Define conditions for excitation, inhibition, and unresponsive pACC
excitation_pACC_rows_redreward = find(valuesCol85_redreward == 1 & valuesCol74_redreward <= 0.05);
inhibition_pACC_rows_redreward = find(valuesCol85_redreward == -1 & valuesCol74_redreward <= 0.05);
condition1_redreward = (valuesCol85_redreward == 1 & valuesCol74_redreward <= 0.05);
condition2_redreward = (valuesCol85_redreward == -1 & valuesCol74_redreward <= 0.05);

% Combine conditions to identify unresponsive rows
unresponsive_pACC_rows_redreward = ~(condition1_redreward | condition2_redreward);
totalUnresponsive_pACC_redreward = sum(unresponsive_pACC_rows_redreward);

% Total rows in pACC_matrix2_sorted
totalRows_redreward = size(pACC_matrix2_sorted_redreward, 1);

% Calculate percentages for excitation, inhibition, and unresponsive
percentageExcitation_pACC_redreward = (length(excitation_pACC_rows_redreward) / totalRows_redreward) * 100;
percentageInhibition_pACC_redreward = (length(inhibition_pACC_rows_redreward) / totalRows_redreward) * 100;
percentageUnresponsive_pACC_redreward = (totalUnresponsive_pACC_redreward / totalRows_redreward) * 100;

% Display the counts and percentages
disp(['Number of excitation_pACC rows at Red Reward: ', num2str(length(excitation_pACC_rows_redreward))]);
disp(['Number of inhibition_pACC rows at Red Reward: ', num2str(length(inhibition_pACC_rows_redreward))]);
disp(['Number of unresponsive_pACC rows at Red Reward: ', num2str(totalUnresponsive_pACC_redreward)]);

disp(['Percentage of excitation_pACC rows at Red Reward: ', num2str(percentageExcitation_pACC_redreward), '%']);
disp(['Percentage of inhibition_pACC rows at Red Reward: ', num2str(percentageInhibition_pACC_redreward), '%']);
disp(['Percentage of unresponsive_pACC rows at Red Reward: ', num2str(percentageUnresponsive_pACC_redreward), '%']);

% Extract the values from columns 86 and 75
valuesCol86_yellowairpuff = cell2mat(pACC_matrix2_sorted_yellowairpuff(:, 86));
valuesCol75_yellowairpuff = cell2mat(pACC_matrix2_sorted_yellowairpuff(:, 75));

% Define conditions for excitation, inhibition, and unresponsive pACC
excitation_pACC_rows_yellowairpuff = find(valuesCol86_yellowairpuff == 1 & valuesCol75_yellowairpuff <= 0.05);
inhibition_pACC_rows_yellowairpuff = find(valuesCol86_yellowairpuff == -1 & valuesCol75_yellowairpuff <= 0.05);
condition1_yellowairpuff = (valuesCol86_yellowairpuff == 1 & valuesCol75_yellowairpuff <= 0.05);
condition2_yellowairpuff = (valuesCol86_yellowairpuff == -1 & valuesCol75_yellowairpuff <= 0.05);

% Combine conditions to identify unresponsive rows
unresponsive_pACC_rows_yellowairpuff = ~(condition1_yellowairpuff | condition2_yellowairpuff);
totalUnresponsive_pACC_yellowairpuff = sum(unresponsive_pACC_rows_yellowairpuff);

% Total rows in pACC_matrix2_sorted
totalRows_yellowairpuff = size(pACC_matrix2_sorted_yellowairpuff, 1);

% Calculate percentages for excitation, inhibition, and unresponsive
percentageExcitation_pACC_yellowairpuff = (length(excitation_pACC_rows_yellowairpuff) / totalRows_yellowairpuff) * 100;
percentageInhibition_pACC_yellowairpuff = (length(inhibition_pACC_rows_yellowairpuff) / totalRows_yellowairpuff) * 100;
percentageUnresponsive_pACC_yellowairpuff = (totalUnresponsive_pACC_yellowairpuff / totalRows_yellowairpuff) * 100;

% Display the counts and percentages
disp(['Number of excitation_pACC rows at Yellow Airpuff: ', num2str(length(excitation_pACC_rows_yellowairpuff))]);
disp(['Number of inhibition_pACC rows at Yellow Airpuff: ', num2str(length(inhibition_pACC_rows_yellowairpuff))]);
disp(['Number of unresponsive_pACC rows at Yellow Airpuff: ', num2str(totalUnresponsive_pACC_yellowairpuff)]);

disp(['Percentage of excitation_pACC rows at Yellow Airpuff: ', num2str(percentageExcitation_pACC_yellowairpuff), '%']);
disp(['Percentage of inhibition_pACC rows at Yellow Airpuff: ', num2str(percentageInhibition_pACC_yellowairpuff), '%']);
disp(['Percentage of unresponsive_pACC rows at Yellow Airpuff: ', num2str(percentageUnresponsive_pACC_yellowairpuff), '%']);

predictorNames={'Reward','Aversion','Eutility','Choice','Reward*Choice', 'Aversion*Choice','Conflict', 'Reaction Time'};
numPredictors = length(predictorNames); % Number of predictors

sigma = 1; % Standard deviation for Gaussian kernel (adjust as needed)
comeon_smoothed = imgaussfilt(comeon, sigma);

% Apply Gaussian smoothing to the comeon matrix
sigma = 1; % Standard deviation for Gaussian kernel (adjust as needed)
comeon_smoothed_air = imgaussfilt(comeon_air, sigma);

%%%% for small reward
% Apply Gaussian smoothing to the comeon matrix
sigma = 1; % Standard deviation for Gaussian kernel (adjust as needed)
comeon_smoothed_smallreward = imgaussfilt(comeon_smallreward, sigma);

%%%% for big reward
% Apply Gaussian smoothing to the comeon matrix
sigma = 1; % Standard deviation for Gaussian kernel (adjust as needed)
comeon_smoothed_bigreward = imgaussfilt(comeon_bigreward, sigma);

% Apply Gaussian smoothing to the comeon1 matrix
sigma = 1; % Standard deviation for Gaussian kernel (adjust as needed)
comeon1_smoothed = imgaussfilt(comeon1, sigma);

% Apply Gaussian smoothing to the comeon1_air matrix
sigma = 1; % Standard deviation for Gaussian kernel (adjust as needed)
comeon1_smoothed_air = imgaussfilt(comeon1_air, sigma);

% Apply Gaussian smoothing to the comeon1_smallreward matrix
sigma = 1; % Standard deviation for Gaussian kernel (adjust as needed)
comeon1_smoothed_smallreward = imgaussfilt(comeon1_smallreward, sigma);

% Apply Gaussian smoothing to the comeon1_bigreward matrix
sigma = 1; % Standard deviation for Gaussian kernel (adjust as needed)
comeon1_smoothed_bigreward = imgaussfilt(comeon1_bigreward, sigma);

% Extract the values from columns 64 and 67 for cOFC
valuesCol64_cOFC = cell2mat(cOFC_matrix2_sorted(:, 64));
valuesCol67_cOFC = cell2mat(cOFC_matrix2_sorted(:, 67));

% Define conditions for excitation, inhibition, and unresponsive cOFC
excitation_cOFC_rows = find(valuesCol64_cOFC == 1 & valuesCol67_cOFC <= 0.05);
inhibition_cOFC_rows = find(valuesCol64_cOFC == -1 & valuesCol67_cOFC <= 0.05);
condition1_cOFC = (valuesCol64_cOFC == 1 & valuesCol67_cOFC <= 0.05);
condition2_cOFC = (valuesCol64_cOFC == -1 & valuesCol67_cOFC <= 0.05);

% Combine conditions to identify unresponsive rows for cOFC
unresponsive_cOFC_rows = ~(condition1_cOFC | condition2_cOFC);
totalUnresponsive_cOFC = sum(unresponsive_cOFC_rows);

% Total rows in cOFC_matrix2_sorted
totalRows_cOFC = size(cOFC_matrix2_sorted, 1);

% Calculate percentages for excitation, inhibition, and unresponsive for cOFC
percentageExcitation_cOFC = (length(excitation_cOFC_rows) / totalRows_cOFC) * 100;
percentageInhibition_cOFC = (length(inhibition_cOFC_rows) / totalRows_cOFC) * 100;
percentageUnresponsive_cOFC = (totalUnresponsive_cOFC / totalRows_cOFC) * 100;

% Display the counts and percentages for cOFC
disp(['Number of excitation_cOFC rows: ', num2str(length(excitation_cOFC_rows))]);
disp(['Number of inhibition_cOFC rows: ', num2str(length(inhibition_cOFC_rows))]);
disp(['Number of unresponsive_cOFC rows: ', num2str(totalUnresponsive_cOFC)]);

disp(['Percentage of excitation_cOFC rows: ', num2str(percentageExcitation_cOFC), '%']);
disp(['Percentage of inhibition_cOFC rows: ', num2str(percentageInhibition_cOFC), '%']);
disp(['Percentage of unresponsive_cOFC rows: ', num2str(percentageUnresponsive_cOFC), '%']);

% Extract the values from columns 64 and 67 for cOFC
valuesCol65_cOFC_air = cell2mat(cOFC_matrix2_sorted_air(:, 65));
valuesCol68_cOFC_air = cell2mat(cOFC_matrix2_sorted_air(:, 68));

% Define conditions for excitation, inhibition, and unresponsive cOFC
excitation_cOFC_rows_air = find(valuesCol65_cOFC_air == 1 & valuesCol68_cOFC_air <= 0.05);
inhibition_cOFC_rows_air = find(valuesCol65_cOFC_air == -1 & valuesCol68_cOFC_air <= 0.05);
condition1_cOFC_air = (valuesCol65_cOFC_air == 1 & valuesCol68_cOFC_air <= 0.05);
condition2_cOFC_air = (valuesCol65_cOFC_air == -1 & valuesCol68_cOFC_air <= 0.05);

% Combine conditions to identify unresponsive rows for cOFC
unresponsive_cOFC_rows_air = ~(condition1_cOFC_air | condition2_cOFC_air);
totalUnresponsive_cOFC_air = sum(unresponsive_cOFC_rows_air);

% Total rows in cOFC_matrix2_sorted
totalRows_cOFC_air = size(cOFC_matrix2_sorted_air, 1);

% Calculate percentages for excitation, inhibition, and unresponsive for cOFC
percentageExcitation_cOFC_air = (length(excitation_cOFC_rows_air) / totalRows_cOFC_air) * 100;
percentageInhibition_cOFC_air = (length(inhibition_cOFC_rows_air) / totalRows_cOFC_air) * 100;
percentageUnresponsive_cOFC_air = (totalUnresponsive_cOFC_air / totalRows_cOFC_air) * 100;

% Display the counts and percentages for cOFC
disp(['Number of excitation_cOFC rows at Approach Outcome: ', num2str(length(excitation_cOFC_rows_air))]);
disp(['Number of inhibition_cOFC rows at Approach Outcome: ', num2str(length(inhibition_cOFC_rows_air))]);
disp(['Number of unresponsive_cOFC rows at Approach Outcome: ', num2str(totalUnresponsive_cOFC_air)]);

disp(['Percentage of excitation_cOFC rows at Approach Outcome: ', num2str(percentageExcitation_cOFC_air), '%']);
disp(['Percentage of inhibition_cOFC rows at Approach Outcome: ', num2str(percentageInhibition_cOFC_air), '%']);
disp(['Percentage of unresponsive_cOFC rows at Approach Outcome: ', num2str(percentageUnresponsive_cOFC_air), '%']);

% Extract the values from columns 80 and 69 for cOFC
valuesCol80_cOFC_smallreward = cell2mat(cOFC_matrix2_sorted_smallreward(:, 80));
valuesCol69_cOFC_smallreward = cell2mat(cOFC_matrix2_sorted_smallreward(:, 69));

% Define conditions for excitation, inhibition, and unresponsive cOFC
excitation_cOFC_rows_smallreward = find(valuesCol80_cOFC_smallreward == 1 & valuesCol69_cOFC_smallreward <= 0.05);
inhibition_cOFC_rows_smallreward = find(valuesCol80_cOFC_smallreward == -1 & valuesCol69_cOFC_smallreward <= 0.05);
condition1_cOFC_smallreward = (valuesCol80_cOFC_smallreward == 1 & valuesCol69_cOFC_smallreward <= 0.05);
condition2_cOFC_smallreward = (valuesCol80_cOFC_smallreward == -1 & valuesCol69_cOFC_smallreward <= 0.05);

% Combine conditions to identify unresponsive rows for cOFC
unresponsive_cOFC_rows_smallreward = ~(condition1_cOFC_smallreward | condition2_cOFC_smallreward);
totalUnresponsive_cOFC_smallreward = sum(unresponsive_cOFC_rows_smallreward);

% Total rows in cOFC_matrix2_sorted
totalRows_cOFC_smallreward = size(cOFC_matrix2_sorted_smallreward, 1);

% Calculate percentages for excitation, inhibition, and unresponsive for cOFC
percentageExcitation_cOFC_smallreward = (length(excitation_cOFC_rows_smallreward) / totalRows_cOFC_smallreward) * 100;
percentageInhibition_cOFC_smallreward = (length(inhibition_cOFC_rows_smallreward) / totalRows_cOFC_smallreward) * 100;
percentageUnresponsive_cOFC_smallreward = (totalUnresponsive_cOFC_smallreward / totalRows_cOFC_smallreward) * 100;

% Display the counts and percentages for cOFC
disp(['Number of excitation_cOFC rows at Avoidance Outcome: ', num2str(length(excitation_cOFC_rows_smallreward))]);
disp(['Number of inhibition_cOFC rows at Avoidance Outcome: ', num2str(length(inhibition_cOFC_rows_smallreward))]);
disp(['Number of unresponsive_cOFC rows at Avoidance Outcome: ', num2str(totalUnresponsive_cOFC_smallreward)]);

disp(['Percentage of excitation_cOFC rows at Avoidance Outcome: ', num2str(percentageExcitation_cOFC_smallreward), '%']);
disp(['Percentage of inhibition_cOFC rows at Avoidance Outcome: ', num2str(percentageInhibition_cOFC_smallreward), '%']);
disp(['Percentage of unresponsive_cOFC rows at Avoidance Outcome: ', num2str(percentageUnresponsive_cOFC_smallreward), '%']);

% Extract the values from columns 87 and 76 for cOFC
valuesCol87_cOFC_bigreward = cell2mat(cOFC_matrix2_sorted_bigreward(:, 87));
valuesCol76_cOFC_bigreward = cell2mat(cOFC_matrix2_sorted_bigreward(:, 76));

% Define conditions for excitation, inhibition, and unresponsive cOFC
excitation_cOFC_rows_bigreward = find(valuesCol87_cOFC_bigreward == 1 & valuesCol76_cOFC_bigreward <= 0.05);
inhibition_cOFC_rows_bigreward = find(valuesCol87_cOFC_bigreward == -1 & valuesCol76_cOFC_bigreward <= 0.05);
condition1_cOFC_bigreward = (valuesCol87_cOFC_bigreward == 1 & valuesCol76_cOFC_bigreward <= 0.05);
condition2_cOFC_bigreward = (valuesCol87_cOFC_bigreward == -1 & valuesCol76_cOFC_bigreward <= 0.05);

% Combine conditions to identify unresponsive rows for cOFC
unresponsive_cOFC_rows_bigreward = ~(condition1_cOFC_bigreward | condition2_cOFC_bigreward);
totalUnresponsive_cOFC_bigreward = sum(unresponsive_cOFC_rows_bigreward);

% Total rows in cOFC_matrix2_sorted
totalRows_cOFC_bigreward = size(cOFC_matrix2_sorted_bigreward, 1);

% Calculate percentages for excitation, inhibition, and unresponsive for cOFC
percentageExcitation_cOFC_bigreward = (length(excitation_cOFC_rows_bigreward) / totalRows_cOFC_bigreward) * 100;
percentageInhibition_cOFC_bigreward = (length(inhibition_cOFC_rows_bigreward) / totalRows_cOFC_bigreward) * 100;
percentageUnresponsive_cOFC_bigreward = (totalUnresponsive_cOFC_bigreward / totalRows_cOFC_bigreward) * 100;

% Display the counts and percentages for cOFC
disp(['Number of excitation_cOFC rows at Approach Reward Outcome: ', num2str(length(excitation_cOFC_rows_bigreward))]);
disp(['Number of inhibition_cOFC rows at Approach Reward Outcome: ', num2str(length(inhibition_cOFC_rows_bigreward))]);
disp(['Number of unresponsive_cOFC rows at Approach Reward Outcome: ', num2str(totalUnresponsive_cOFC_bigreward)]);

disp(['Percentage of excitation_cOFC rows at Approach Reward Outcome: ', num2str(percentageExcitation_cOFC_bigreward), '%']);
disp(['Percentage of inhibition_cOFC rows at Approach Reward Outcome: ', num2str(percentageInhibition_cOFC_bigreward), '%']);
disp(['Percentage of unresponsive_cOFC rows at Approach Reward Outcome: ', num2str(percentageUnresponsive_cOFC_bigreward), '%']);

% Extract the values from columns 83 and 72 for cOFC
valuesCol83_cOFC_redcue = cell2mat(cOFC_matrix2_sorted_redcue(:, 83));
valuesCol72_cOFC_redcue = cell2mat(cOFC_matrix2_sorted_redcue(:, 72));

% Define conditions for excitation, inhibition, and unresponsive cOFC
excitation_cOFC_rows_redcue = find(valuesCol83_cOFC_redcue == 1 & valuesCol72_cOFC_redcue <= 0.05);
inhibition_cOFC_rows_redcue = find(valuesCol83_cOFC_redcue == -1 & valuesCol72_cOFC_redcue <= 0.05);
condition1_cOFC_redcue = (valuesCol83_cOFC_redcue == 1 & valuesCol72_cOFC_redcue <= 0.05);
condition2_cOFC_redcue = (valuesCol83_cOFC_redcue == -1 & valuesCol72_cOFC_redcue <= 0.05);

% Combine conditions to identify unresponsive rows for cOFC
unresponsive_cOFC_rows_redcue = ~(condition1_cOFC_redcue | condition2_cOFC_redcue);
totalUnresponsive_cOFC_redcue = sum(unresponsive_cOFC_rows_redcue);

% Total rows in cOFC_matrix2_sorted
totalRows_cOFC_redcue = size(cOFC_matrix2_sorted_redcue, 1);

% Calculate percentages for excitation, inhibition, and unresponsive for cOFC
percentageExcitation_cOFC_redcue = (length(excitation_cOFC_rows_redcue) / totalRows_cOFC_redcue) * 100;
percentageInhibition_cOFC_redcue = (length(inhibition_cOFC_rows_redcue) / totalRows_cOFC_redcue) * 100;
percentageUnresponsive_cOFC_redcue = (totalUnresponsive_cOFC_redcue / totalRows_cOFC_redcue) * 100;

% Display the counts and percentages for cOFC
disp(['Number of excitation_cOFC rows at Red Cue: ', num2str(length(excitation_cOFC_rows_redcue))]);
disp(['Number of inhibition_cOFC rows at Red Cue: ', num2str(length(inhibition_cOFC_rows_redcue))]);
disp(['Number of unresponsive_cOFC rows at Red Cue: ', num2str(totalUnresponsive_cOFC_redcue)]);

disp(['Percentage of excitation_cOFC rows at Red Cue: ', num2str(percentageExcitation_cOFC_redcue), '%']);
disp(['Percentage of inhibition_cOFC rows at Red Cue: ', num2str(percentageInhibition_cOFC_redcue), '%']);
disp(['Percentage of unresponsive_cOFC rows at Red Cue: ', num2str(percentageUnresponsive_cOFC_redcue), '%']);

% Extract the values from columns 84 and 73 for cOFC
valuesCol84_cOFC_yellowcue = cell2mat(cOFC_matrix2_sorted_yellowcue(:, 84));
valuesCol73_cOFC_yellowcue = cell2mat(cOFC_matrix2_sorted_yellowcue(:, 73));

% Define conditions for excitation, inhibition, and unresponsive cOFC
excitation_cOFC_rows_yellowcue = find(valuesCol84_cOFC_yellowcue == 1 & valuesCol73_cOFC_yellowcue <= 0.05);
inhibition_cOFC_rows_yellowcue = find(valuesCol84_cOFC_yellowcue == -1 & valuesCol73_cOFC_yellowcue <= 0.05);
condition1_cOFC_yellowcue = (valuesCol84_cOFC_yellowcue == 1 & valuesCol73_cOFC_yellowcue <= 0.05);
condition2_cOFC_yellowcue = (valuesCol84_cOFC_yellowcue == -1 & valuesCol73_cOFC_yellowcue <= 0.05);

% Combine conditions to identify unresponsive rows for cOFC
unresponsive_cOFC_rows_yellowcue = ~(condition1_cOFC_yellowcue | condition2_cOFC_yellowcue);
totalUnresponsive_cOFC_yellowcue = sum(unresponsive_cOFC_rows_yellowcue);

% Total rows in cOFC_matrix2_sorted
totalRows_cOFC_yellowcue = size(cOFC_matrix2_sorted_yellowcue, 1);

% Calculate percentages for excitation, inhibition, and unresponsive for cOFC
percentageExcitation_cOFC_yellowcue = (length(excitation_cOFC_rows_yellowcue) / totalRows_cOFC_yellowcue) * 100;
percentageInhibition_cOFC_yellowcue = (length(inhibition_cOFC_rows_yellowcue) / totalRows_cOFC_yellowcue) * 100;
percentageUnresponsive_cOFC_yellowcue = (totalUnresponsive_cOFC_yellowcue / totalRows_cOFC_yellowcue) * 100;

% Display the counts and percentages for cOFC
disp(['Number of excitation_cOFC rows at Yellow Cue: ', num2str(length(excitation_cOFC_rows_yellowcue))]);
disp(['Number of inhibition_cOFC rows at Yellow Cue: ', num2str(length(inhibition_cOFC_rows_yellowcue))]);
disp(['Number of unresponsive_cOFC rows at Yellow Cue: ', num2str(totalUnresponsive_cOFC_yellowcue)]);

disp(['Percentage of excitation_cOFC rows at Yellow Cue: ', num2str(percentageExcitation_cOFC_yellowcue), '%']);
disp(['Percentage of inhibition_cOFC rows at Yellow Cue: ', num2str(percentageInhibition_cOFC_yellowcue), '%']);
disp(['Percentage of unresponsive_cOFC rows at Yellow Cue: ', num2str(percentageUnresponsive_cOFC_yellowcue), '%']);

% Extract the values from columns 85 and 74 for cOFC
valuesCol85_cOFC_redreward = cell2mat(cOFC_matrix2_sorted_redreward(:, 85));
valuesCol74_cOFC_redreward = cell2mat(cOFC_matrix2_sorted_redreward(:, 74));

% Define conditions for excitation, inhibition, and unresponsive cOFC
excitation_cOFC_rows_redreward = find(valuesCol85_cOFC_redreward == 1 & valuesCol74_cOFC_redreward <= 0.05);
inhibition_cOFC_rows_redreward = find(valuesCol85_cOFC_redreward == -1 & valuesCol74_cOFC_redreward <= 0.05);
condition1_cOFC_redreward = (valuesCol85_cOFC_redreward == 1 & valuesCol74_cOFC_redreward <= 0.05);
condition2_cOFC_redreward = (valuesCol85_cOFC_redreward == -1 & valuesCol74_cOFC_redreward <= 0.05);

% Combine conditions to identify unresponsive rows for cOFC
unresponsive_cOFC_rows_redreward = ~(condition1_cOFC_redreward | condition2_cOFC_redreward);
totalUnresponsive_cOFC_redreward = sum(unresponsive_cOFC_rows_redreward);

% Total rows in cOFC_matrix2_sorted
totalRows_cOFC_redreward = size(cOFC_matrix2_sorted_redreward, 1);

% Calculate percentages for excitation, inhibition, and unresponsive for cOFC
percentageExcitation_cOFC_redreward = (length(excitation_cOFC_rows_redreward) / totalRows_cOFC_redreward) * 100;
percentageInhibition_cOFC_redreward = (length(inhibition_cOFC_rows_redreward) / totalRows_cOFC_redreward) * 100;
percentageUnresponsive_cOFC_redreward = (totalUnresponsive_cOFC_redreward / totalRows_cOFC_redreward) * 100;

% Display the counts and percentages for cOFC
disp(['Number of excitation_cOFC rows at Red Reward: ', num2str(length(excitation_cOFC_rows_redreward))]);
disp(['Number of inhibition_cOFC rows at Red Reward: ', num2str(length(inhibition_cOFC_rows_redreward))]);
disp(['Number of unresponsive_cOFC rows at Red Reward: ', num2str(totalUnresponsive_cOFC_redreward)]);

disp(['Percentage of excitation_cOFC rows at Red Reward: ', num2str(percentageExcitation_cOFC_redreward), '%']);
disp(['Percentage of inhibition_cOFC rows at Red Reward: ', num2str(percentageInhibition_cOFC_redreward), '%']);
disp(['Percentage of unresponsive_cOFC rows at Red Reward: ', num2str(percentageUnresponsive_cOFC_redreward), '%']);

% Extract the values from columns 86 and 75 for cOFC
valuesCol86_cOFC_yellowairpuff = cell2mat(cOFC_matrix2_sorted_yellowairpuff(:, 86));
valuesCol75_cOFC_yellowairpuff = cell2mat(cOFC_matrix2_sorted_yellowairpuff(:, 75));

% Define conditions for excitation, inhibition, and unresponsive cOFC
excitation_cOFC_rows_yellowairpuff = find(valuesCol86_cOFC_yellowairpuff == 1 & valuesCol75_cOFC_yellowairpuff <= 0.05);
inhibition_cOFC_rows_yellowairpuff = find(valuesCol86_cOFC_yellowairpuff == -1 & valuesCol75_cOFC_yellowairpuff <= 0.05);
condition1_cOFC_yellowairpuff = (valuesCol86_cOFC_yellowairpuff == 1 & valuesCol75_cOFC_yellowairpuff <= 0.05);
condition2_cOFC_yellowairpuff = (valuesCol86_cOFC_yellowairpuff == -1 & valuesCol75_cOFC_yellowairpuff <= 0.05);

% Combine conditions to identify unresponsive rows for cOFC
unresponsive_cOFC_rows_yellowairpuff = ~(condition1_cOFC_yellowairpuff | condition2_cOFC_yellowairpuff);
totalUnresponsive_cOFC_yellowairpuff = sum(unresponsive_cOFC_rows_yellowairpuff);

% Total rows in cOFC_matrix2_sorted
totalRows_cOFC_yellowairpuff = size(cOFC_matrix2_sorted_yellowairpuff, 1);

% Calculate percentages for excitation, inhibition, and unresponsive for cOFC
percentageExcitation_cOFC_yellowairpuff = (length(excitation_cOFC_rows_yellowairpuff) / totalRows_cOFC_yellowairpuff) * 100;
percentageInhibition_cOFC_yellowairpuff = (length(inhibition_cOFC_rows_yellowairpuff) / totalRows_cOFC_yellowairpuff) * 100;
percentageUnresponsive_cOFC_yellowairpuff = (totalUnresponsive_cOFC_yellowairpuff / totalRows_cOFC_yellowairpuff) * 100;

% Display the counts and percentages for cOFC
disp(['Number of excitation_cOFC rows at Red Reward: ', num2str(length(excitation_cOFC_rows_yellowairpuff))]);
disp(['Number of inhibition_cOFC rows at Red Reward: ', num2str(length(inhibition_cOFC_rows_yellowairpuff))]);
disp(['Number of unresponsive_cOFC rows at Red Reward: ', num2str(totalUnresponsive_cOFC_yellowairpuff)]);

disp(['Percentage of excitation_cOFC rows at Red Reward: ', num2str(percentageExcitation_cOFC_yellowairpuff), '%']);
disp(['Percentage of inhibition_cOFC rows at Red Reward: ', num2str(percentageInhibition_cOFC_yellowairpuff), '%']);
disp(['Percentage of unresponsive_cOFC rows at Red Reward: ', num2str(percentageUnresponsive_cOFC_yellowairpuff), '%']);



% Create a new figure with specified dimensions for better visualization
figure('Position', [100, 100, 1800, 1200]); % Adjust size as needed

% === Subplot 1: pACC - Cue Period ===
subplot(2, 4, 1);
imagesc(comeon1_smoothed); % Assuming comeon1_smoothed is the data matrix for pACC - Cue Period
caxis([-3, 3]);
title('pACC: Cue Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Draw vertical and horizontal lines as needed with hold on/off appropriately
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Add the horizontal dotted lines based on certain lengths and conditions
line(xlim, [length(pACC_spikes_CuePeriod_excluded_normalized_filtered) length(pACC_spikes_CuePeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(pACC_spikes_CuePeriod_excluded_normalized_filtered)+length(excitation_pACC_rows) length(pACC_spikes_CuePeriod_excluded_normalized_filtered)+length(excitation_pACC_rows)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line
hold off;

% === Subplot 2: pACC - Approach Outcome Period ===
ax2 = subplot(2, 4, 2);
imagesc(comeon1_smoothed_air); % Assuming comeon1_smoothed_air is another data matrix for pACC - Approach Outcome Period
caxis([-3, 3]);
title('Approach Airpuff Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Draw vertical and horizontal lines as needed with hold on/off appropriately
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Add the horizontal dotted lines based on certain lengths and conditions
line(xlim, [length(pACC_spikes_AirpuffPeriod_excluded_normalized_filtered) length(pACC_spikes_AirpuffPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(pACC_spikes_AirpuffPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_air) length(pACC_spikes_AirpuffPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_air)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line
hold off;
% Repeat the plotting commands for the second plot here

% === Subplot 3: pACC - Approach Reward Period ===
subplot(2, 4, 3);
imagesc(comeon1_smoothed_bigreward); % Assuming this is another data matrix for pACC - Avoidance Outcome Period
caxis([-3, 3]);
title('Approach Reward Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Draw vertical and horizontal lines as needed with hold on/off appropriately
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Add the horizontal dotted lines based on certain lengths and conditions
line(xlim, [length(pACC_spikes_BigRewardPeriod_excluded_normalized_filtered) length(pACC_spikes_BigRewardPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(pACC_spikes_BigRewardPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_bigreward) length(pACC_spikes_BigRewardPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_bigreward)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line
hold off;


% === Subplot 3: pACC - Avoidance Outcome Period ===
subplot(2, 4, 4);
imagesc(comeon1_smoothed_smallreward); % Assuming this is another data matrix for pACC - Avoidance Outcome Period
caxis([-3, 3]);
title('Avoidance Outcome Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Draw vertical and horizontal lines as needed with hold on/off appropriately
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Add the horizontal dotted lines based on certain lengths and conditions
line(xlim, [length(pACC_spikes_SmallRewardPeriod_excluded_normalized_filtered) length(pACC_spikes_SmallRewardPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(pACC_spikes_SmallRewardPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_smallreward) length(pACC_spikes_SmallRewardPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_smallreward)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line
hold off;
% Repeat the plotting commands for the third plot here

% === Subplot 4: cOFC - Cue Period ===
subplot(2, 4, 5);
imagesc(comeon_smoothed); % Assuming comeon_smoothed is the data matrix for cOFC - Cue Period
caxis([-3, 3]);
title('cOFC: Cue Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
hold on;

% Assuming cOFC_spikes_CuePeriod_excluded_normalized_filtered and onlyPlus_noNaN are predefined
% Add horizontal dotted lines at specified positions based on the length calculations
line(xlim, [length(cOFC_spikes_CuePeriod_excluded_normalized_filtered) length(cOFC_spikes_CuePeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line for excluded normalized filtered
line(xlim, [length(cOFC_spikes_CuePeriod_excluded_normalized_filtered)+length(onlyPlus_noNaN) length(cOFC_spikes_CuePeriod_excluded_normalized_filtered)+length(onlyPlus_noNaN)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line for onlyPlus_noNaN

hold off; % Release the plot hold
% === Subplot 5: cOFC - Approach Outcome Period ===
subplot(2, 4, 6);
imagesc(comeon_smoothed_air); % Assuming comeon_smoothed_air is another data matrix for cOFC - Approach Outcome Period
caxis([-3,3]);
title('Approach Airpuff Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
hold on;
% Assuming cOFC_spikes_CuePeriod_excluded_normalized_filtered and onlyPlus_noNaN are predefined
% Add horizontal dotted lines at specified positions based on the length calculations
line(xlim, [length(cOFC_spikes_AirpuffPeriod_excluded_normalized_filtered) length(cOFC_spikes_AirpuffPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line for excluded normalized filtered
line(xlim, [length(cOFC_spikes_AirpuffPeriod_excluded_normalized_filtered)+length(onlyPlus_noNaN_air) length(cOFC_spikes_AirpuffPeriod_excluded_normalized_filtered)+length(onlyPlus_noNaN_air)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line for onlyPlus_noNaN_air

hold off; % Release the plot hold

% === Subplot 6: cOFC - Approach Reward Outcome Period ===
subplot(2, 4, 7);
imagesc(comeon_smoothed_bigreward); % Assuming comeon_smoothed_bigreward is data matrix for cOFC - Avoidance Outcome Period
caxis([-3, 3]);
title('Approach Reward Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
hold on;
line(xlim, [length(cOFC_spikes_BigRewardPeriod_excluded_normalized_filtered) length(cOFC_spikes_BigRewardPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line for excluded normalized filtered
line(xlim, [length(cOFC_spikes_BigRewardPeriod_excluded_normalized_filtered)+length(onlyPlus_noNaN_bigreward) length(cOFC_spikes_BigRewardPeriod_excluded_normalized_filtered)+length(onlyPlus_noNaN_bigreward)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line for onlyPlus_noNaN_bigreward
hold off

% === Subplot 6: cOFC - Avoidance Outcome Period ===
subplot(2, 4, 8);
imagesc(comeon_smoothed_smallreward); % Assuming comeon_smoothed_smallreward is data matrix for cOFC - Avoidance Outcome Period
caxis([-3, 3]);
title('Avoidance Outcome Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
hold on;
line(xlim, [length(cOFC_spikes_SmallRewardPeriod_excluded_normalized_filtered) length(cOFC_spikes_SmallRewardPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line for excluded normalized filtered
line(xlim, [length(cOFC_spikes_SmallRewardPeriod_excluded_normalized_filtered)+length(onlyPlus_noNaN_smallreward) length(cOFC_spikes_SmallRewardPeriod_excluded_normalized_filtered)+length(onlyPlus_noNaN_smallreward)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line for onlyPlus_noNaN_smallreward
hold off

% Make sure the figure uses the 'jet' colormap globally if desired
colormap(jet);

% Adjust figure background color if needed
set(gcf, 'Color', 'w'); 
comeon1_redcue=[pACC_spikes_RedCuePeriod_excluded_normalized_filtered;pACC_spikes_RedCuePeriod_normalized_filtered];
comeon1_yellowcue=[pACC_spikes_YellowCuePeriod_excluded_normalized_filtered;pACC_spikes_YellowCuePeriod_normalized_filtered];
comeon1_redreward=[pACC_spikes_RedRewardPeriod_excluded_normalized_filtered;pACC_spikes_RedRewardPeriod_normalized_filtered];
comeon1_yellowairpuff=[pACC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered;pACC_spikes_YellowAirpuffPeriod_normalized_filtered];
comeon_redcue=[cOFC_spikes_RedCuePeriod_excluded_normalized_filtered;cOFC_spikes_RedCuePeriod_normalized_filtered];
comeon_yellowcue=[cOFC_spikes_YellowCuePeriod_excluded_normalized_filtered;cOFC_spikes_YellowCuePeriod_normalized_filtered];
comeon_redreward=[cOFC_spikes_RedRewardPeriod_excluded_normalized_filtered;cOFC_spikes_RedRewardPeriod_normalized_filtered];
comeon_yellowairpuff=[cOFC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered;cOFC_spikes_YellowAirpuffPeriod_normalized_filtered];

sigma = 1; % Standard deviation for Gaussian kernel (adjust as needed)
comeon1_smoothed_redcue = imgaussfilt(comeon1_redcue, sigma);
comeon1_smoothed_yellowcue = imgaussfilt(comeon1_yellowcue, sigma);
comeon1_smoothed_redreward = imgaussfilt(comeon1_redreward, sigma);
comeon1_smoothed_yellowairpuff = imgaussfilt(comeon1_yellowairpuff, sigma);
comeon_smoothed_redcue = imgaussfilt(comeon_redcue, sigma);
comeon_smoothed_yellowcue = imgaussfilt(comeon_yellowcue, sigma);
comeon_smoothed_redreward = imgaussfilt(comeon_redreward, sigma);
comeon_smoothed_yellowairpuff = imgaussfilt(comeon_yellowairpuff, sigma);

% Create a new figure with specified dimensions for better visualization
figure('Position', [100, 100, 1800, 1200]); % Adjust size as needed

% === Subplot 1: pACC - Cue Period ===
subplot(2, 4, 1);
imagesc(comeon1_smoothed_redcue); % Assuming comeon1_redcue is the data matrix for pACC - Cue Period
caxis([-3, 3]);
title('pACC: Reward Cue Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Draw vertical and horizontal lines as needed with hold on/off appropriately
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Add the horizontal dotted lines based on certain lengths and conditions
line(xlim, [length(pACC_spikes_RedCuePeriod_excluded_normalized_filtered) length(pACC_spikes_RedCuePeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(pACC_spikes_RedCuePeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_redcue) length(pACC_spikes_RedCuePeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_redcue)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line
hold off;

% === Subplot 2: pACC - Approach Outcome Period ===
ax2 = subplot(2, 4, 2);
imagesc(comeon1_smoothed_yellowcue); % Assuming comeon1_smoothed_yellowcue is another data matrix for pACC - Approach Outcome Period
caxis([-3, 3]);
title('Airpuff Cue Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Draw vertical and horizontal lines as needed with hold on/off appropriately
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Add the horizontal dotted lines based on certain lengths and conditions
line(xlim, [length(pACC_spikes_YellowCuePeriod_excluded_normalized_filtered) length(pACC_spikes_YellowCuePeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(pACC_spikes_YellowCuePeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_yellowcue) length(pACC_spikes_YellowCuePeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_yellowcue)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line
hold off;
% Repeat the plotting commands for the second plot here

% === Subplot 3: pACC - Approach Reward Period ===
subplot(2, 4, 3);
imagesc(comeon1_smoothed_redreward); % Assuming this is another data matrix for pACC - Pavlovian Reward Period
caxis([-3, 3]);
title('Pavlovian Reward Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Draw vertical and horizontal lines as needed with hold on/off appropriately
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Add the horizontal dotted lines based on certain lengths and conditions
line(xlim, [length(pACC_spikes_RedRewardPeriod_excluded_normalized_filtered) length(pACC_spikes_RedRewardPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(pACC_spikes_RedRewardPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_redreward) length(pACC_spikes_RedRewardPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_redreward)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line
hold off;


% === Subplot 3: pACC - Pavlovian Airpuff Period ===
subplot(2, 4, 4);
imagesc(comeon1_smoothed_yellowairpuff); % Assuming this is another data matrix for pACC - Pavlovian Airpuff Period
caxis([-3, 3]);
title('Pavlovian Airpuff Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Draw vertical and horizontal lines as needed with hold on/off appropriately
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
% Add the horizontal dotted lines based on certain lengths and conditions
line(xlim, [length(pACC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered) length(pACC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(pACC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_yellowairpuff) length(pACC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered)+length(excitation_pACC_rows_yellowairpuff)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line
hold off;
% Repeat the plotting commands for the third plot here

% === Subplot 4: cOFC - Cue Period ===
subplot(2, 4, 5);
imagesc(comeon_smoothed_redcue); % Assuming comeon_smoothed is the data matrix for cOFC - Cue Period
caxis([-3, 3]);
title('cOFC: Reward Cue Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
hold on;

% Assuming cOFC_spikes_RedCuePeriod_excluded_normalized_filtered and onlyPlus_noNaN are predefined
% Add horizontal dotted lines at specified positions based on the length calculations
line(xlim, [length(cOFC_spikes_RedCuePeriod_excluded_normalized_filtered) length(cOFC_spikes_RedCuePeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(cOFC_spikes_RedCuePeriod_excluded_normalized_filtered)+length(excitation_cOFC_rows_redcue) length(cOFC_spikes_RedCuePeriod_excluded_normalized_filtered)+length(excitation_cOFC_rows_redcue)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line

hold off; % Release the plot hold
% === Subplot 5: cOFC - Approach Outcome Period ===
subplot(2, 4, 6);
imagesc(comeon_smoothed_yellowcue); % Assuming comeon_smoothed_yellowcue is another data matrix for cOFC - Approach Outcome Period
caxis([-3,3]);
title('Airpuff Cue Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
hold on;
% Assuming cOFC_spikes_RedCuePeriod_excluded_normalized_filtered and onlyPlus_noNaN are predefined
% Add horizontal dotted lines at specified positions based on the length calculations
line(xlim, [length(cOFC_spikes_YellowCuePeriod_excluded_normalized_filtered) length(cOFC_spikes_YellowCuePeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line
line(xlim, [length(cOFC_spikes_YellowCuePeriod_excluded_normalized_filtered)+length(excitation_cOFC_rows_yellowcue) length(cOFC_spikes_YellowCuePeriod_excluded_normalized_filtered)+length(excitation_cOFC_rows_yellowcue)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line

hold off; % Release the plot hold

% === Subplot 6: cOFC - Approach Reward Outcome Period ===
subplot(2, 4, 7);
imagesc(comeon_smoothed_redreward); % Assuming comeon_smoothed_redreward is data matrix for cOFC - Pavlovian Reward Period
caxis([-3, 3]);
title('Pavlovian Reward Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
hold on;
line(xlim, [length(cOFC_spikes_RedRewardPeriod_excluded_normalized_filtered) length(cOFC_spikes_RedRewardPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line for excluded normalized filtered
line(xlim, [length(cOFC_spikes_RedRewardPeriod_excluded_normalized_filtered)+length(excitation_cOFC_rows_redreward) length(cOFC_spikes_RedRewardPeriod_excluded_normalized_filtered)+length(excitation_cOFC_rows_redreward)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line for onlyPlus_noNaN_redreward
hold off

% === Subplot 6: cOFC - Avoidance Outcome Period ===
subplot(2, 4, 8);
imagesc(comeon_smoothed_yellowairpuff); % Assuming comeon_smoothed_yellowairpuff is data matrix for cOFC - Avoidance Outcome Period
caxis([-3, 3]);
title('Pavlovian Airpuff Period', 'FontSize', 18);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Units #', 'FontSize', 12);
set(gca, 'XTick', linspace(1, 80, 9), 'XTickLabels', linspace(-2, 2, 9), 'FontSize', 18);
colormap(jet);
hold on;
line([40.5 40.5], ylim, 'Color', 'k', 'LineWidth', 2); % Vertical black line at the zero point
hold on;
line(xlim, [length(cOFC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered) length(cOFC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % First dotted line for excluded normalized filtered
line(xlim, [length(cOFC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered)+length(excitation_cOFC_rows_yellowairpuff) length(cOFC_spikes_YellowAirpuffPeriod_excluded_normalized_filtered)+length(excitation_cOFC_rows_yellowairpuff)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2); % Second dotted line for onlyPlus_noNaN_yellowairpuff
hold off

% Make sure the figure uses the 'jet' colormap globally if desired
colormap(jet);

% Adjust figure background color if needed
set(gcf, 'Color', 'w'); 

% Initialize a container for counting occurrences of predictor combinations in pACC
predictorCombinationFrequencies_pACC = containers.Map('KeyType', 'char', 'ValueType', 'double');

% Iterate over the dataset to count occurrences
for row = 1:size(filtered_manipulation, 1)
    if strcmp(filtered_manipulation{row, 6}, 'pACC')  % Check if the row is for cOFC
        predictorsList = filtered_manipulation{row, 62}; % Retrieve the predictors for this row
        if isempty(predictorsList) || all(cellfun(@isempty, predictorsList))
            % Skip cases with no predictors
            continue;
        else
            % Create a key for the current combination of predictors (sorted to ensure consistency)
            combinationKey = strjoin(sort(predictorsList), '+');
            % Increment the frequency count for this combination
            if isKey(predictorCombinationFrequencies_pACC, combinationKey)
                predictorCombinationFrequencies_pACC(combinationKey) = predictorCombinationFrequencies_pACC(combinationKey) + 1;
            else
                predictorCombinationFrequencies_pACC(combinationKey) = 1;
            end
        end
    end
end

% Prepare data for plotting
combinationNames = keys(predictorCombinationFrequencies_pACC);
combinationCounts = cell2mat(values(predictorCombinationFrequencies_pACC));

% Sort combinations by their counts
[sortedCounts_pACC, sortIdx] = sort(combinationCounts, 'descend');
sortedCombinations_pACC = combinationNames(sortIdx);

% Find and exclude the index of combinations with no predictors, if present
noPredictorIndex = find(strcmp(sortedCombinations_pACC, ''), 1);
if ~isempty(noPredictorIndex)
    sortedCounts_pACC(noPredictorIndex) = [];  % Remove count for no predictor combination
    sortedCombinations_pACC(noPredictorIndex) = [];  % Remove label for no predictor combination
end

fig = figure;

% Set figure size
fig.Position = [100, 100, 1000, 600]; 

% Bar plot
ax1 = subplot('Position', [0.1 0.4 0.8 0.5]);
bar(ax1, sortedCounts_pACC, 'FaceColor', 'b');
ylabel('# of Units');
title('Classification of units recorded in pACC');
set(ax1, 'XTick', 1:length(sortedCombinations_pACC), 'XTickLabel', {}, 'XLim', [0.5, length(sortedCombinations_pACC)+0.5]);
set(ax1, 'Box', 'off', 'Color', 'none');

% Add count labels above bars
for i = 1:length(sortedCounts_pACC)
    text(ax1, i, sortedCounts_pACC(i), num2str(sortedCounts_pACC(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% Matrix plot
ax2 = subplot('Position', [0.1 0.19 0.8 0.2]);
hold(ax2, 'on');

% Assuming numPredictors is defined and equal to 8
numPredictors = 8; % Total number of predictors, assuming x1 to x8

% Draw the grid matrix manually based on 'sortedCombinations_cOFC'
for col = 1:length(sortedCombinations_pACC)
    comboParts = split(sortedCombinations_pACC{col}, '+');
    for partIndex = 1:length(comboParts)
        predictorIndex = str2double(extractAfter(comboParts(partIndex), 'x')); % Convert 'x1' to 1, 'x2' to 2, etc.
        % Calculate the Y position to keep x1 at the top and x8 at the bottom
        yPosition = numPredictors - predictorIndex + 0.5;
        % Fill the corresponding cell in the matrix
        rectangle(ax2, 'Position', [col-0.5, yPosition, 1, 1], 'FaceColor', 'k', 'EdgeColor', 'none');
    end
end

% Set matrix plot properties
set(ax2, 'XLim', [0.5, length(sortedCombinations_pACC)+0.5], 'YLim', [0.5, numPredictors+0.5]);
set(ax2, 'YTick', 1:numPredictors, 'YTickLabel', flip(predictorNames), 'TickLength', [0 0]);
set(ax2, 'XTick', [], 'XTickLabel', {});
set(ax2, 'Box', 'off', 'Color', 'none');

% Draw grid lines
for i = 1:numPredictors+1
    line(ax2, [0.5, length(sortedCombinations_pACC)+0.5], [i-0.5, i-0.5], 'Color', 'k');
end
for col = 1:length(sortedCombinations_pACC)+1
    line(ax2, [col-0.5, col-0.5], [0.5, numPredictors+0.5], 'Color', 'k');
end

hold(ax2, 'off');

% Set size change callback to adjust the figure when resized
set(fig, 'SizeChangedFcn', @(src, evnt) resizeFig(src, ax1, ax2));

matchingRows = [];

% Iterate over each row in filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Check if column 6 contains 'pACC'
    if strcmp(filtered_manipulation{i, 6}, 'pACC')
        % Check the number of predictors in column 62
        numPredictors = size(filtered_manipulation{i, 62}, 1);
        
        % If there is exactly one predictor
        if numPredictors == 1
            % This row meets the criteria, add it to matchingRows
            matchingRows = [matchingRows; i];
        end
    end
end

% Display the number of matching rows
disp(['Number of matching rows with exactly one predictor and "pACC": ', num2str(length(matchingRows))]);

% Plotting
figure;
bar(sortedCounts_pACC, 'FaceColor', 'b');  % Color the bars blue
set(gca, 'XTick', 1:length(sortedCombinations_pACC), 'XTickLabel', sortedCombinations_pACC, 'XTickLabelRotation', 45);
ylabel('# of Units');
title('Classification of units recorded in pACC');
set(gcf, 'Color', 'w');  % Set background color to white
box off;

% Optionally, add counts on top of each bar for clarity
for i = 1:length(sortedCounts_pACC)
    text(i, sortedCounts_pACC(i), num2str(sortedCounts_pACC(i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom');
end



% Adding descriptions for x1, x2, ..., x8 with actual meanings
descriptions = {
    'x1: Reward', 
    'x2: Aversion', 
    'x3: Eutility', 
    'x4: Choice', 
    'x5: Reward*Choice', 
    'x6: Aversion*Choice', 
    'x7: Conflict', 
    'x8: Reaction Time'
};
descText = strjoin(descriptions, '\n');

% Place the description text on the top right of the plot
dim = [0.65 0.6 0.3 0.3]; % Normalized units for textbox position and size
annotation('textbox', dim, 'String', descText, 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'none');
% Extracting the names next to x1:, x2:, etc., for the predictors names

% Initialize a container for counting occurrences of predictor combinations in cOFC
predictorCombinationFrequencies_cOFC = containers.Map('KeyType', 'char', 'ValueType', 'double');

% Iterate over the dataset to count occurrences
for row = 1:size(filtered_manipulation, 1)
    if strcmp(filtered_manipulation{row, 6}, 'cOFC')  % Check if the row is for cOFC
        predictorsList = filtered_manipulation{row, 62}; % Retrieve the predictors for this row
        if isempty(predictorsList) || all(cellfun(@isempty, predictorsList))
            % Skip cases with no predictors
            continue;
        else
            % Create a key for the current combination of predictors (sorted to ensure consistency)
            combinationKey = strjoin(sort(predictorsList), '+');
            % Increment the frequency count for this combination
            if isKey(predictorCombinationFrequencies_cOFC, combinationKey)
                predictorCombinationFrequencies_cOFC(combinationKey) = predictorCombinationFrequencies_cOFC(combinationKey) + 1;
            else
                predictorCombinationFrequencies_cOFC(combinationKey) = 1;
            end
        end
    end
end

% Prepare data for plotting
combinationNames = keys(predictorCombinationFrequencies_cOFC);
combinationCounts = cell2mat(values(predictorCombinationFrequencies_cOFC));

% Sort combinations by their counts
[sortedCounts, sortIdx] = sort(combinationCounts, 'descend');
sortedCombinations_cOFC = combinationNames(sortIdx);

% Find and exclude the index of combinations with no predictors, if present
noPredictorIndex = find(strcmp(sortedCombinations_cOFC, ''), 1);
if ~isempty(noPredictorIndex)
    sortedCounts(noPredictorIndex) = [];  % Remove count for no predictor combination
    sortedCombinations_cOFC(noPredictorIndex) = [];  % Remove label for no predictor combination
end

matchingRows = [];

% Iterate over each row in filtered_manipulation
for i = 1:size(filtered_manipulation, 1)
    % Check if column 6 contains 'cOFC'
    if strcmp(filtered_manipulation{i, 6}, 'cOFC')
        % Check the number of predictors in column 62
        numPredictors = size(filtered_manipulation{i, 62}, 1);
        
        % If there is exactly one predictor
        if numPredictors == 1
            % This row meets the criteria, add it to matchingRows
            matchingRows = [matchingRows; i];
        end
    end
end

% Display the number of matching rows
disp(['Number of matching rows with exactly one predictor and "cOFC": ', num2str(length(matchingRows))]);

%%%% with matrix
% Predictor names for 8 predictors
predictorNames = {'Reward', 'Aversion', 'Eutility', 'Choice', 'Reward*Choice', 'Aversion*Choice', 'Conflict', 'Reaction Time'};
numPredictors = length(predictorNames);

% Assuming 'sortedCounts' contains the frequency for each combination
% Assuming 'sortedCombinations_cOFC' contains strings representing the predictors for each bar

% Create the figure
fig = figure;

% Set figure size
fig.Position = [100, 100, 1000, 600]; 

% Bar plot
ax1 = subplot('Position', [0.1 0.4 0.8 0.5]);
bar(ax1, sortedCounts, 'FaceColor', 'b');
ylabel('# of Units');
title('Classification of units recorded in cOFC');
set(ax1, 'XTick', 1:length(sortedCombinations_cOFC), 'XTickLabel', {}, 'XLim', [0.5, length(sortedCombinations_cOFC)+0.5]);
set(ax1, 'Box', 'off', 'Color', 'none');

% Add count labels above bars
for i = 1:length(sortedCounts)
    text(ax1, i, sortedCounts(i), num2str(sortedCounts(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% Matrix plot
ax2 = subplot('Position', [0.1 0.19 0.8 0.2]);
hold(ax2, 'on');

% Assuming numPredictors is defined and equal to 8
numPredictors = 8; % Total number of predictors, assuming x1 to x8

% Draw the grid matrix manually based on 'sortedCombinations_cOFC'
for col = 1:length(sortedCombinations_cOFC)
    comboParts = split(sortedCombinations_cOFC{col}, '+');
    for partIndex = 1:length(comboParts)
        predictorIndex = str2double(extractAfter(comboParts(partIndex), 'x')); % Convert 'x1' to 1, 'x2' to 2, etc.
        % Calculate the Y position to keep x1 at the top and x8 at the bottom
        yPosition = numPredictors - predictorIndex + 0.5;
        % Fill the corresponding cell in the matrix
        rectangle(ax2, 'Position', [col-0.5, yPosition, 1, 1], 'FaceColor', 'k', 'EdgeColor', 'none');
    end
end

% Set matrix plot properties
set(ax2, 'XLim', [0.5, length(sortedCombinations_cOFC)+0.5], 'YLim', [0.5, numPredictors+0.5]);
set(ax2, 'YTick', 1:numPredictors, 'YTickLabel', flip(predictorNames), 'TickLength', [0 0]);
set(ax2, 'XTick', [], 'XTickLabel', {});
set(ax2, 'Box', 'off', 'Color', 'none');

% Draw grid lines
for i = 1:numPredictors+1
    line(ax2, [0.5, length(sortedCombinations_cOFC)+0.5], [i-0.5, i-0.5], 'Color', 'k');
end
for col = 1:length(sortedCombinations_cOFC)+1
    line(ax2, [col-0.5, col-0.5], [0.5, numPredictors+0.5], 'Color', 'k');
end

hold(ax2, 'off');

% Set size change callback to adjust the figure when resized
set(fig, 'SizeChangedFcn', @(src, evnt) resizeFig(src, ax1, ax2));

% Set size change callback to adjust the figure when resized
set(fig, 'SizeChangedFcn', @(src, evnt) resizeFig(src, ax1, ax2));

% Assuming extendedMatrix2_clean is the cleaned data matrix
% Define the base column for 'rew30co_all' which is column 95
baseColumn = 95;

% Define the list of category suffixes
categories = {'rew30co', 'rew60co', 'rew90co', 'rew120co', 'rew150co', 'rew180co', ...
              'air30co', 'air60co', 'air90co', 'air120co', 'air150co', 'air180co', ...
              'rew30rd', 'rew60rd', 'rew90rd', 'rew120rd', 'rew150rd', 'rew180rd', ...
              'air30rd', 'air60rd', 'air90rd', 'air120rd', 'air150rd', 'air180rd'};

% Define regions
regions = {'cofc', 'pacc'};
region_data_means = struct();

for reg = regions
    % Find rows containing current region in column 6
    contains_region = cellfun(@(x) contains(lower(strtrim(x)), reg{1}), extendedMatrix2_clean(:, 6));
    
    % Initialize a structure to hold the means for current region
    region_means = struct();

    % Loop over each category to extract and concatenate data
    for i = 1:length(categories)
        currentColumn = baseColumn + (i - 1);  % Calculate current column
        cells_from_current_col = extendedMatrix2_clean(contains_region, currentColumn);

        % Filter out any empty cells before concatenation
        non_empty_cells = cells_from_current_col(~cellfun(@isempty, cells_from_current_col));

        % Concatenate non-empty cells into a matrix and calculate the mean
        if ~isempty(non_empty_cells)
            category_data = vertcat(non_empty_cells{:});
            region_means.(categories{i}) = mean(category_data, 1);  % Calculate mean across rows
        else
            disp(['No data to process for ' categories{i} ' in region ' reg{1} '.']);
        end
    end
    
    % Store region means in the main structure
    region_data_means.(reg{1}) = region_means;
end

% Define the groups as per the categories
groups = {
    {'rew30co', 'rew60co', 'rew90co', 'rew120co', 'rew150co', 'rew180co'},
    {'air30co', 'air60co', 'air90co', 'air120co', 'air150co', 'air180co'},
    {'rew30rd', 'rew60rd', 'rew90rd', 'rew120rd', 'rew150rd', 'rew180rd'},
    {'air30rd', 'air60rd', 'air90rd', 'air120rd', 'air150rd', 'air180rd'}
};
groupNames = {'Reward Cue', 'Airpuff Cue', 'Reward Delivery', 'Airpuff Delivery'};

% Assuming extendedMatrix2_clean is the cleaned data matrix
% Define the base column for 'rew30co_all' which is column 95
baseColumn = 95;

% Define the list of category suffixes
categories = {'rew30co', 'rew60co', 'rew90co', 'rew120co', 'rew150co', 'rew180co', ...
              'air30co', 'air60co', 'air90co', 'air120co', 'air150co', 'air180co', ...
              'rew30rd', 'rew60rd', 'rew90rd', 'rew120rd', 'rew150rd', 'rew180rd', ...
              'air30rd', 'air60rd', 'air90rd', 'air120rd', 'air150rd', 'air180rd'};

% Define regions
regions = {'cofc', 'pacc'};
region_data_means = struct();

for reg = regions
    % Find rows containing current region in column 6
    contains_region = cellfun(@(x) contains(lower(strtrim(x)), reg{1}), extendedMatrix2_clean(:, 6));
    
    % Initialize a structure to hold the means for current region
    region_means = struct();

    % Loop over each category to extract and concatenate data
    for i = 1:length(categories)
        currentColumn = baseColumn + (i - 1);  % Calculate current column
        cells_from_current_col = extendedMatrix2_clean(contains_region, currentColumn);

        % Filter out any empty cells before concatenation
        non_empty_cells = cells_from_current_col(~cellfun(@isempty, cells_from_current_col));

        % Concatenate non-empty cells into a matrix and calculate the mean
        if ~isempty(non_empty_cells)
            category_data = vertcat(non_empty_cells{:});
            category_data_mean = mean(category_data, 1);  % Calculate mean across rows
            
            % Adjustment: Subtract the mean of columns 21 to 40 from all data points
            if size(category_data, 1) >= 40  % Ensure there are at least 40 columns
                mean_21_to_40 = mean(category_data(:, 21:40), 'all');
                adjusted_mean = category_data_mean - mean_21_to_40;
            else
                adjusted_mean = category_data_mean - mean(category_data_mean);  % Use overall mean if fewer than 40 columns
            end

            region_means.(categories{i}) = adjusted_mean;
        else
            disp(['No data to process for ' categories{i} ' in region ' reg{1} '.']);
        end
    end
    
    % Store region means in the main structure
    region_data_means.(reg{1}) = region_means;
end

% Define the groups as per the categories
groups = {
    {'rew30co', 'rew60co', 'rew90co', 'rew120co', 'rew150co', 'rew180co'},
    {'air30co', 'air60co', 'air90co', 'air120co', 'air150co', 'air180co'},
    {'rew30rd', 'rew60rd', 'rew90rd', 'rew120rd', 'rew150rd', 'rew180rd'},
    {'air30rd', 'air60rd', 'air90rd', 'air120rd', 'air150rd', 'air180rd'}
};

groupNames = {'Reward Cue', 'Airpuff Cue', 'Reward Delivery', 'Airpuff Delivery'};

% Create a figure with 8 subplots
figure;
for i = 1:length(groups)
    for reg = regions
        subplot(4, 2, 2*(i-1) + find(strcmp(regions, reg)))
        hold on;  % Hold on to add multiple lines to the same plot
        for j = 1:length(groups{i})
            category = groups{i}{j};
            if isfield(region_data_means.(reg{1}), category)
                plot(region_data_means.(reg{1}).(category));  % Plot the adjusted mean values
            else
                disp(['Means for ' category ' in ' reg{1} ' not calculated or not available.']);
            end
        end
        % Add a red vertical line at index 41
        xline(41, 'r', 'LineWidth', 2);  % Draw a red vertical line
        hold off;
        title([reg{1} ' - ' groupNames{i}]);  % Title for each group
        legend(groups{i}, 'Location', 'bestoutside');  % Add a legend
        xlabel('Index');  % Assuming each index corresponds to some measurement index
        ylabel('Adjusted Mean Value');  % Adjusted mean values on the Y-axis
    end
end

% Assuming extendedMatrix2_clean is the cleaned data matrix
% Define the base column for 'rew30co_all' which is column 95
baseColumn = 95;

% Define the list of category suffixes
categories = {'rew30co', 'rew60co', 'rew90co', 'rew120co', 'rew150co', 'rew180co', ...
              'air30co', 'air60co', 'air90co', 'air120co', 'air150co', 'air180co', ...
              'rew30rd', 'rew60rd', 'rew90rd', 'rew120rd', 'rew150rd', 'rew180rd', ...
              'air30rd', 'air60rd', 'air90rd', 'air120rd', 'air150rd', 'air180rd'};

% Define regions
regions = {'cofc', 'pacc'};
region_data_means = struct();

for reg = regions
    % Find rows containing current region in column 6
    contains_region = cellfun(@(x) contains(lower(strtrim(x)), reg{1}), extendedMatrix2_clean(:, 6));
    
    % Initialize a structure to hold the means for current region
    region_means = struct();

    % Loop over each category to extract and concatenate data
    for i = 1:length(categories)
        currentColumn = baseColumn + (i - 1);  % Calculate current column
        cells_from_current_col = extendedMatrix2_clean(contains_region, currentColumn);

        % Filter out any empty cells before concatenation
        non_empty_cells = cells_from_current_col(~cellfun(@isempty, cells_from_current_col));

        % Concatenate non-empty cells into a matrix and calculate the mean
        if ~isempty(non_empty_cells)
            category_data = vertcat(non_empty_cells{:});
            category_data_mean = mean(category_data, 1);  % Calculate mean across rows
            
            % Adjustment: Subtract the mean of columns 21 to 40 from all data points
            if size(category_data, 1) >= 40  % Ensure there are at least 40 columns
                mean_21_to_40 = mean(category_data(:, 21:40), 'all');
                adjusted_mean = category_data_mean - mean_21_to_40;
            else
                adjusted_mean = category_data_mean - mean(category_data_mean);  % Use overall mean if fewer than 40 columns
            end

            region_means.(categories{i}) = adjusted_mean;
        else
            disp(['No data to process for ' categories{i} ' in region ' reg{1} '.']);
        end
    end
    
    % Store region means in the main structure
    region_data_means.(reg{1}) = region_means;
end

% Define the groups as per the categories
groups = {
    {'rew30co', 'rew60co', 'rew90co', 'rew120co', 'rew150co', 'rew180co'},
    {'air30co', 'air60co', 'air90co', 'air120co', 'air150co', 'air180co'},
    {'rew30rd', 'rew60rd', 'rew90rd', 'rew120rd', 'rew150rd', 'rew180rd'},
    {'air30rd', 'air60rd', 'air90rd', 'air120rd', 'air150rd', 'air180rd'}
};

groupNames = {'Reward Cue', 'Airpuff Cue', 'Reward Delivery', 'Airpuff Delivery'};

% Create a figure with 8 subplots
figure;
for i = 1:length(groups)
    for reg = regions
        subplot(4, 2, 2*(i-1) + find(strcmp(regions, reg)))
        hold on;  % Hold on to add multiple lines to the same plot
        for j = 1:length(groups{i})
            category = groups{i}{j};
            if isfield(region_data_means.(reg{1}), category)
                plot(region_data_means.(reg{1}).(category));  % Plot the adjusted mean values
            else
                disp(['Means for ' category ' in ' reg{1} ' not calculated or not available.']);
            end
        end
        % Add a red vertical line at index 41
        xline(41, 'r', 'LineWidth', 2);  % Draw a red vertical line
        hold off;
        set(gca, 'FontSize', 14);  % Change font size for axes
        title([reg{1} ' - ' groupNames{i}], 'FontSize', 20);  % Title with increased font size
        legend(groups{i}, 'Location', 'bestoutside', 'FontSize', 12);  % Legend with specific font size
        xlabel('Time (ms)', 'FontSize', 14);  % X-label with increased font size
        ylabel('Normalized Mean Value (Hz)', 'FontSize', 14);  % Y-label with increased font size
    end
end

% Create a figure with 8 subplots
figure;
for i = 1:length(groups)
    for reg = regions
        subplot(4, 2, 2*(i-1) + find(strcmp(regions, reg)))
        hold on;  % Hold on to add multiple lines to the same plot
        for j = 1:length(groups{i})
            category = groups{i}{j};
            if isfield(region_data_means.(reg{1}), category)
                plot(region_data_means.(reg{1}).(category));  % Plot the mean values
            else
                disp(['Means for ' category ' in ' reg{1} ' not calculated or not available.']);
            end
        end
        hold off;
        title([reg{1} ' - ' groupNames{i}]);  % Title for each group
        legend(groups{i}, 'Location', 'bestoutside');  % Add a legend
        xlabel('Index');  % Assuming each index corresponds to some measurement index
        ylabel('Mean Value');  % Mean values on the Y-axis
    end
end

% Assuming the matrix is named extendedMatrix2_clean
numRows = size(extendedMatrix2_clean, 1);
numbers = 1:8; % Define the range of numbers to check

% Initialize counters for each number for each condition
positive_cOFC_counts = zeros(size(numbers));
negative_cOFC_counts = zeros(size(numbers));
positive_pACC_counts = zeros(size(numbers));
negative_pACC_counts = zeros(size(numbers));

% Loop through each row of the matrix
for i = 1:numRows
    region = extendedMatrix2_clean{i, 6};
    number = extendedMatrix2_clean{i, 119};
    value = extendedMatrix2_clean{i, 120};

    % Check if the number is within the range we are interested in
    if ismember(number, numbers)
        % Check the conditions for cOFC
        if strcmp(region, 'cOFC')
            if value > 0
                positive_cOFC_counts(number) = positive_cOFC_counts(number) + 1;
            elseif value < 0
                negative_cOFC_counts(number) = negative_cOFC_counts(number) + 1;
            end
        end

        % Check the conditions for pACC
        if strcmp(region, 'pACC')
            if value > 0
                positive_pACC_counts(number) = positive_pACC_counts(number) + 1;
            elseif value < 0
                negative_pACC_counts(number) = negative_pACC_counts(number) + 1;
            end
        end
    end
end

% Output or save the counts
disp('Counts for positive_cOFC:');
disp(array2table(positive_cOFC_counts', 'RowNames', cellstr(num2str(numbers'))));

disp('Counts for negative_cOFC:');
disp(array2table(negative_cOFC_counts', 'RowNames', cellstr(num2str(numbers'))));

disp('Counts for positive_pACC:');
disp(array2table(positive_pACC_counts', 'RowNames', cellstr(num2str(numbers'))));

disp('Counts for negative_pACC:');
disp(array2table(negative_pACC_counts', 'RowNames', cellstr(num2str(numbers'))));

% Create the figure
fig = figure;

% Set figure size
fig.Position = [100, 100, 1000, 600]; 

% Initialize data for two types of bars: stacked and white
stackedBarData = zeros(length(sortedCombinations_cOFC), 2); % For positive and negative counts
whiteBarData = NaN(length(sortedCombinations_cOFC), 1); % For combinations

% Define colors and properties
positiveColor = 'r'; % Color for positive counts
negativeColor = 'b'; % Color for negative counts
combinationColor = [0 0.5 0];%[1 1 1]; % White color for combinations

% Analyze combinations and assign data to bar types
for i = 1:length(sortedCombinations_cOFC)
    if contains(sortedCombinations_cOFC{i}, '+') % Check if it is a combination
        whiteBarData(i) = sum(sortedCounts(i)); % Use total count for white bars
    else % Single predictors
        predictorIndex = str2double(extractAfter(sortedCombinations_cOFC{i}, 'x'));
        if predictorIndex <= length(positive_cOFC_counts) && predictorIndex <= length(negative_cOFC_counts)
            stackedBarData(i, 1) = positive_cOFC_counts(predictorIndex); % Red part
            stackedBarData(i, 2) = negative_cOFC_counts(predictorIndex); % Blue part
        end
    end
end

% Subplot for bar plots
ax1 = subplot('Position', [0.1 0.4 0.8 0.5]);
hold(ax1, 'on');

% Plot white bars for combinations
hb1 = bar(ax1, whiteBarData, 'FaceColor', combinationColor, 'EdgeColor', 'k');

% Plot stacked bars for single predictors
hb2 = bar(ax1, stackedBarData, 'stacked');
set(hb2(1), 'FaceColor', positiveColor); % Red for positive counts
set(hb2(2), 'FaceColor', negativeColor); % Blue for negative counts

% Additional plot properties
ylabel('# of Units');
title('Classification of units recorded in cOFC', 'FontSize', 22);
set(ax1, 'XTick', 1:length(sortedCombinations_cOFC), 'XTickLabel', {}, 'XLim', [0.5, length(sortedCombinations_cOFC)+0.5]);
set(ax1, 'Box', 'off', 'Color', 'none');

% Add count labels above bars
for i = 1:length(sortedCombinations_cOFC)
    totalHeight = nansum([whiteBarData(i), sum(stackedBarData(i, :))]);
    if ~isnan(totalHeight)
        text(ax1, i, totalHeight, num2str(totalHeight), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

hold(ax1, 'off');

% Matrix plot
ax2 = subplot('Position', [0.1 0.19 0.8 0.2]);
hold(ax2, 'on');

% Assuming numPredictors is defined and equal to 8
numPredictors = 8; % Total number of predictors, assuming x1 to x8

% Draw the grid matrix manually based on 'sortedCombinations_cOFC'
for col = 1:length(sortedCombinations_cOFC)
    comboParts = split(sortedCombinations_cOFC{col}, '+');
    for partIndex = 1:length(comboParts)
        predictorIndex = str2double(extractAfter(comboParts(partIndex), 'x')); % Convert 'x1' to 1, 'x2' to 2, etc.
        % Calculate the Y position to keep x1 at the top and x8 at the bottom
        yPosition = numPredictors - predictorIndex + 0.5;
        % Fill the corresponding cell in the matrix
        rectangle(ax2, 'Position', [col-0.5, yPosition, 1, 1], 'FaceColor', 'k', 'EdgeColor', 'none');
    end
end

% Set matrix plot properties
set(ax2, 'XLim', [0.5, length(sortedCombinations_cOFC)+0.5], 'YLim', [0.5, numPredictors+0.5]);
set(ax2, 'YTick', 1:numPredictors, 'YTickLabel', flip(predictorNames), 'TickLength', [0 0]);
set(ax2, 'XTick', [], 'XTickLabel', {});
set(ax2, 'Box', 'off', 'Color', 'none');

% Draw grid lines
for i = 1:numPredictors+1
    line(ax2, [0.5, length(sortedCombinations_cOFC)+0.5], [i-0.5, i-0.5], 'Color', 'k');
end
for col = 1:length(sortedCombinations_cOFC)+1
    line(ax2, [col-0.5, col-0.5], [0.5, numPredictors+0.5], 'Color', 'k');
end

hold(ax2, 'off');

% Create the figure
fig = figure;

% Set figure size
fig.Position = [100, 100, 1000, 600]; 

% Initialize data for two types of bars: stacked and white
stackedBarData = zeros(length(sortedCombinations_pACC), 2); % For positive and negative counts
whiteBarData = NaN(length(sortedCombinations_pACC), 1); % For combinations

% Define colors and properties
positiveColor = 'r'; % Color for positive counts
negativeColor = 'b'; % Color for negative counts
combinationColor = [0 0.5 0];%[1 1 1]; % White color for combinations

% Analyze combinations and assign data to bar types
for i = 1:length(sortedCombinations_pACC)
    if contains(sortedCombinations_pACC{i}, '+') % Check if it is a combination
        whiteBarData(i) = sum(sortedCounts_pACC(i)); % Use total count for white bars
    else % Single predictors
        predictorIndex = str2double(extractAfter(sortedCombinations_pACC{i}, 'x'));
        if predictorIndex <= length(positive_pACC_counts) && predictorIndex <= length(negative_pACC_counts)
            stackedBarData(i, 1) = positive_pACC_counts(predictorIndex); % Red part
            stackedBarData(i, 2) = negative_pACC_counts(predictorIndex); % Blue part
        end
    end
end

% Subplot for bar plots
ax1 = subplot('Position', [0.1 0.4 0.8 0.5]);
hold(ax1, 'on');

% Plot white bars for combinations
hb1 = bar(ax1, whiteBarData, 'FaceColor', combinationColor, 'EdgeColor', 'k');

% Plot stacked bars for single predictors
hb2 = bar(ax1, stackedBarData, 'stacked');
set(hb2(1), 'FaceColor', positiveColor); % Red for positive counts
set(hb2(2), 'FaceColor', negativeColor); % Blue for negative counts

% Additional plot properties
ylabel('# of Units');
title('Classification of units recorded in pACC', 'FontSize', 22);
set(ax1, 'XTick', 1:length(sortedCombinations_pACC), 'XTickLabel', {}, 'XLim', [0.5, length(sortedCombinations_pACC)+0.5]);
set(ax1, 'Box', 'off', 'Color', 'none');

% Add count labels above bars
for i = 1:length(sortedCombinations_pACC)
    totalHeight = nansum([whiteBarData(i), sum(stackedBarData(i, :))]);
    if ~isnan(totalHeight)
        text(ax1, i, totalHeight, num2str(totalHeight), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

hold(ax1, 'off');

% Matrix plot
ax2 = subplot('Position', [0.1 0.19 0.8 0.2]);
hold(ax2, 'on');

% Assuming numPredictors is defined and equal to 8
numPredictors = 8; % Total number of predictors, assuming x1 to x8

% Draw the grid matrix manually based on 'sortedCombinations_pACC'
for col = 1:length(sortedCombinations_pACC)
    comboParts = split(sortedCombinations_pACC{col}, '+');
    for partIndex = 1:length(comboParts)
        predictorIndex = str2double(extractAfter(comboParts(partIndex), 'x')); % Convert 'x1' to 1, 'x2' to 2, etc.
        % Calculate the Y position to keep x1 at the top and x8 at the bottom
        yPosition = numPredictors - predictorIndex + 0.5;
        % Fill the corresponding cell in the matrix
        rectangle(ax2, 'Position', [col-0.5, yPosition, 1, 1], 'FaceColor', 'k', 'EdgeColor', 'none');
    end
end

% Set matrix plot properties
set(ax2, 'XLim', [0.5, length(sortedCombinations_pACC)+0.5], 'YLim', [0.5, numPredictors+0.5]);
set(ax2, 'YTick', 1:numPredictors, 'YTickLabel', flip(predictorNames), 'TickLength', [0 0]);
set(ax2, 'XTick', [], 'XTickLabel', {});
set(ax2, 'Box', 'off', 'Color', 'none');

% Draw grid lines
for i = 1:numPredictors+1
    line(ax2, [0.5, length(sortedCombinations_pACC)+0.5], [i-0.5, i-0.5], 'Color', 'k');
end
for col = 1:length(sortedCombinations_pACC)+1
    line(ax2, [col-0.5, col-0.5], [0.5, numPredictors+0.5], 'Color', 'k');
end

hold(ax2, 'off');

figureBar = figure;

% Initialize storage for bar values across all groups and regions
allBarValues = [];

for i = 1:length(groups)
    for reg = regions
        % Select subplot for bar plots
        subplot(4, 2, 2*(i-1) + find(strcmp(regions, reg)), 'Parent', figureBar);
        hold on;  % Hold on to add bars to the same plot
        barValues = [];  % Initialize array to store mean values for bar plot
        
        for j = 1:length(groups{i})
            category = groups{i}{j};
            if isfield(region_data_means.(reg{1}), category)
                % Retrieve the data for the current category
                data = region_data_means.(reg{1}).(category);
                % Ensure there is enough data to cover up to the 60th index
                if length(data) >= 60
                    % Extract data from indices 21 to 40
                    mean1 = mean(data(21:40));
                    % Extract data from indices 41 to 60
                    mean2 = mean(data(41:60));
                    % Calculate the difference of the means
                    mean_difference = mean2 - mean1;
                    % Store the mean difference for the bar plot
                    barValues(end+1) = mean_difference;
                else
                    % If not enough data, handle it (e.g., set to NaN or skip)
                    barValues(end+1) = NaN;  % Assign NaN for missing data
                end
            else
                disp(['Means for ' category ' in ' reg{1} ' not calculated or not available.']);
                barValues(end+1) = NaN;  % Assign NaN for missing data
            end
        end
        
        % Plot bar chart for the current group and region
        bar(barValues, 'FaceColor', 'b');
        
        % Fit and plot a trend line over the bar values
        x = 1:length(barValues);  % X coordinates for the bars
        validIndices = ~isnan(barValues);  % Exclude NaN values for fitting
        if sum(validIndices) > 2  % Ensure there are enough points to fit
            [p, S] = polyfit(x(validIndices), barValues(validIndices), 1);  % Fit a linear polynomial, obtain residuals
            yFit = polyval(p, x);  % Evaluate the polynomial
            plot(x, yFit, 'r-', 'LineWidth', 2);  % Draw the line in red over the bars

            % Calculate the standard error of the slope and the t-statistic
            x_bar = mean(x(validIndices));
            n = sum(validIndices);
            se = sqrt(S.normr^2 / (n-2) / sum((x(validIndices) - x_bar).^2));
            t_stat = p(1) / se;  % p(1) is the slope
            p_value = 2 * (1 - tcdf(abs(t_stat), n - 2));  % Two-tailed p-value

            % Display significance with group name
            sigMsg = sprintf('%s - %s: Slope t-stat = %.2f, p-value = %.3f', reg{1}, groupNames{i}, t_stat, p_value);
            disp(sigMsg);
        end
        
        hold off;
        
        title([reg{1} ' - ' groupNames{i} ' Change in Mean Values'], 'FontSize', 18);
        ylabel('Change in Mean Value', 'FontSize', 14);
        set(gca, 'XTickLabel', groups{i}, 'XTick', 1:numel(groups{i}), 'FontSize', 14);
        xtickangle(45);  % Angle the x-axis labels for better readability

        % Optionally, store all bar values for global operations or comparisons
        allBarValues = [allBarValues; barValues];  % Append current bar values
    end
end

%%%% SECOND WAY
% Create the figure for bar plots
figureBar = figure;

% Initialize storage for bar values across all groups and regions
allBarValues = [];

for i = 1:length(groups)
    for reg = regions
        % Select subplot for bar plots
        subplot(4, 2, 2*(i-1) + find(strcmp(regions, reg)), 'Parent', figureBar);
        hold on;  % Hold on to add bars to the same plot
        barValues = [];  % Initialize array to store mean differences for bar plot
        
        for j = 1:length(groups{i})
            category = groups{i}{j};
            if isfield(region_data_means.(reg{1}), category)
                % Retrieve the data for the current category
                data = region_data_means.(reg{1}).(category);
                % Ensure there is enough data to cover up to the 60th index
                if length(data) >= 60
                    % Extract data from indices 21 to 40
                    mean1 = mean(data(21:40));
                    % Extract data from indices 41 to 60
                    mean2 = mean(data(41:60));
                    % Calculate the difference of the means
                    mean_difference = mean2 - mean1;
                    % Store the mean difference for the bar plot
                    barValues(end+1) = mean_difference;
                else
                    % If not enough data, handle it (e.g., set to NaN or skip)
                    barValues(end+1) = NaN;  % Assign NaN for missing data
                end
            else
                disp(['Means for ' category ' in ' reg{1} ' not calculated or not available.']);
                barValues(end+1) = NaN;  % Assign NaN for missing data
            end
        end
        
        % Plot bar chart for the current group and region
        b = bar(barValues, 'FaceColor', 'b');
        
        % Fit linear model to barValues to test significance using t-statistic
        x = 1:length(barValues);  % X coordinates for the bars
        validIndices = ~isnan(barValues);  % Exclude NaN values for fitting
        
        if sum(validIndices) > 2  % Ensure there are enough points to fit
            % Create a table for linear model
            tbl = table(x(validIndices)', barValues(validIndices)', 'VariableNames', {'X', 'Y'});
            lm = fitlm(tbl, 'Y ~ X');
            
            % Extract t-statistic and p-value for the slope
            t_stat = lm.Coefficients.tStat(2); % t-statistic for the slope
            p_value = lm.Coefficients.pValue(2); % p-value for the t-statistic
            
            % Plot the fitted trend line
            yFit = lm.Fitted;
            plot(x(validIndices), yFit, 'r-', 'LineWidth', 2);  % Draw the line in red over the bars
            
            % Determine significance markers based on t-statistic p-value
            if p_value < 0.001
                sigMarker = '**';
            elseif p_value < 0.05
                sigMarker = '*';
            else
                sigMarker = '';
            end
            
            % Display significance message in the terminal using t-statistic
            sigMsg = sprintf('%s - %s: t-stat = %.2f, p-value = %.3f', reg{1}, groupNames{i}, t_stat, p_value);
            disp(sigMsg);
            
            % Add significance marker to the plot
            if ~isempty(sigMarker)
                text(mean(x(validIndices)), max(barValues(validIndices)) + 5, sigMarker, ...
                    'FontSize', 18, 'Color', 'k', 'HorizontalAlignment', 'center');
            end
        else
            disp(['Not enough valid bar values to perform trend analysis for ' reg{1} ' - ' groupNames{i}]);
        end
        
        hold off;
        
        % Set plot titles and labels
        title([reg{1} ' - ' groupNames{i} ' Change in Mean Values'], 'FontSize', 18);
        ylabel('Change in Mean Value', 'FontSize', 14);
        set(gca, 'XTickLabel', groups{i}, 'XTick', 1:numel(groups{i}), 'FontSize', 14);
        xtickangle(45);  % Angle the x-axis labels for better readability

        % Optionally, store all bar values for global operations or comparisons
        allBarValues = [allBarValues; barValues];  % Append current bar values
    end
end

%%%%
% Create a figure for line plots
figure;

for i = 1:length(groups)
    for reg = regions
        % Select the subplot for line plots
        subplot(4, 2, 2*(i-1) + find(strcmp(regions, reg)));
        hold on;  % Hold on to add multiple lines to the same plot
        for j = 1:length(groups{i})
            category = groups{i}{j};
            if isfield(region_data_means.(reg{1}), category)
                % Retrieve the data for the current category
                data = region_data_means.(reg{1}).(category);
                % Check if data length is at least 5
                if length(data) >= 5
                    % Calculate the average of the first 5 data points
                    baseline_avg = mean(data(1:5));
                    % Normalize the data by subtracting the baseline average
                    normalized_data = data - baseline_avg;
                else
                    normalized_data = data; % If not enough data, skip normalization
                end
                % Plot the normalized data
                plot(normalized_data);
            else
                disp(['Means for ' category ' in ' reg{1} ' not calculated or not available.']);
            end
        end
        % Draw a vertical line at index 40
        xline(40, 'r-', 'LineWidth', 1.5);  % Red line with a line width of 1.5
        hold off;
        title([reg{1} ' - ' groupNames{i}],'FontSize', 18);  % Title for each group
        legend(groups{i}, 'Location', 'bestoutside');  % Add a legend
        xlabel('Time (ms)','Fontsize',14);  % Assuming each index corresponds to some measurement index
        ylabel('Normalized Mean Value (Hz)','Fontsize',14);  % Baseline-normalized mean values on the Y-axis
    end
end

% Create the contingency table
observed = [significant_cOFC; significant_pACC];

% Sum the observed frequencies for each event across both brain regions
total_cOFC = sum(cOFC);
total_pACC = sum(pACC);
total_events = sum(observed, 1);
total_overall = sum(total_events);

% Calculate the expected frequencies
expected_cOFC = total_cOFC * total_events / total_overall;
expected_pACC = total_pACC * total_events / total_overall;
expected = [expected_cOFC; expected_pACC];

% Perform the Chi-Square Test of Independence
chi2_stat = sum((observed(:) - expected(:)).^2 ./ expected(:));
df = (size(observed, 1) - 1) * (size(observed, 2) - 1);
p_value = 1 - chi2cdf(chi2_stat, df);

% Display the results
disp(['Chi-Square Statistic: ', num2str(chi2_stat)]);
disp(['P-value: ', num2str(p_value)]);
disp(['Degrees of Freedom: ', num2str(df)]);

if p_value < 0.05
    disp('There is a significant difference between the brain regions across task events.');
else
    disp('There is no significant difference between the brain regions across task events.');
end

% Calculate percentages
percent_cOFC = (cOFC / total_cOFC) * 100;
percent_pACC = (pACC / total_pACC) * 100;

% Perform Chi-Square test for each event
chi2_stat = zeros(1, length(cOFC));
p_value = zeros(1, length(cOFC));
significant = zeros(1, length(cOFC));

for i = 1:length(cOFC)
    % Create the contingency table
    observed = [cOFC(i), total_cOFC - cOFC(i); pACC(i), total_pACC - pACC(i)];
    
    % Perform Chi-Square test
    [h, p, stats] = chi2test(observed);
    chi2_stat(i) = stats.chi2stat;
    p_value(i) = p;
    significant(i) = p < 0.05;
end

% Display results
fprintf('Event\tChi-Square Statistic\tP-value\tSignificant\tPercentage cOFC\tPercentage pACC\n');
for i = 1:length(cOFC)
    fprintf('%d\t%.4f\t\t\t%.4f\t%d\t%.2f%%\t\t%.2f%%\n', i, chi2_stat(i), p_value(i), significant(i), percent_significant_cOFC(i), percent_significant_pACC(i));
end

% Specific event of interest
event_index = 3;
fprintf('\nFor Event %d (Approach Airpuff Onset):\n', event_index);
fprintf('Chi-Square Statistic = %.4f\n', chi2_stat(event_index));
fprintf('P-value = %.4f\n', p_value(event_index));
fprintf('Significant = %d\n', significant(event_index));
fprintf('Percentage of cOFC neurons = %.2f%%\n', percent_cOFC(event_index));
fprintf('Percentage of pACC neurons = %.2f%%\n', percent_pACC(event_index));

% Chi-Square test function
function [h, p, stats] = chi2test(observed)
    % Calculate expected frequencies
    total = sum(observed(:));
    expected = sum(observed, 2) * sum(observed, 1) / total;

    % Perform the chi-square test
    chi2stat = sum((observed - expected).^2 ./ expected, 'all');
    df = (size(observed, 1) - 1) * (size(observed, 2) - 1);
    p = 1 - chi2cdf(chi2stat, df);
    
    % Check if the result is significant
    h = p < 0.05;
    
    % Return stats
    stats.chi2stat = chi2stat;
    stats.df = df;
    stats.p = p;
end

% Function to adjust the figure when resized
function resizeFig(src, ax1, ax2)
    src.Units = 'pixels';
    pos = src.Position;
    matrixHeight = pos(4) * 0.2;
    barHeight = pos(4) * 0.6;
    ax1.Position = [0.1, matrixHeight + (pos(4) * 0.1), 0.8, barHeight];
    ax2.Position = [0.1, 0.05 * pos(4), 0.8, matrixHeight];
end


% Function to calculate ratio and avoid division by zero
function ratio = calculateRatio(count1, count_minus1)
    if count_minus1 ~= 0
        ratio = count1 / count_minus1;
    else
        ratio = NaN; % Avoid division by zero
    end
end
