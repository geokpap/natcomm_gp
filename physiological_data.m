clc;
clear;
close all;

currentFigure = 1;
subplotIndex = 1;

load('/Users/georgiospapageorgiou/Dropbox (MIT)/Analyses/PhysiologicalData/Victoria/ApAv/nonclean_ApAv_withLick.mat');
load('/Users/georgiospapageorgiou/Dropbox (MIT)/Analyses/PhysiologicalData/Victoria/AllApAv_aggroTrialPupil_20240909_nodupes.mat');


% Define the rows to be removed
rowsToRemove = [17, 55, 124, 125, 126, 130, 140, 143, 147];

% Filter out rows that are beyond the current number of rows in matrixData
rowsToRemove = rowsToRemove(rowsToRemove <= size(matrixData, 1));

% Remove the specified rows from matrixData
matrixData(rowsToRemove, :) = [];
% matrixData=matrixData(1:50,:);

% Extract the date part of the string for each cell in the second column
Alldates1 = cellfun(@(x) x(1:10), matrixData(:, 2), 'UniformOutput', false);

% Now extractedDates is a cell array with each element being a 'YYYY-MM-DD' string
% Loop through each session in matrixData
for i = 1:size(matrixData, 1)
    % Get the current session's data
    currentSessionData = matrixData{i, 3};

    % Find the rows where all elements are zero
    rowsWithAllZeros = all(currentSessionData == 0, 2);
    
    % Remove those rows
    currentSessionData(rowsWithAllZeros, :) = [];
    
    % Place the filtered data back into matrixData
    matrixData{i, 3} = currentSessionData;
end

for i = 1 : size(matrixData,1)
    Red_bar_perc=matrixData{i,3}(:,1)./2;
    Yellow_bar_perc=matrixData{i,3}(:,2)./2;
    Chosen_option=matrixData{i,3}(:,3);


X0=[Red_bar_perc,Yellow_bar_perc];
y=Chosen_option;

[b,dev,stats]=glmfit(X0,y,'binomial','link','logit');

betaValues{1,i}=b;
    x=0:100;
    plot_x = [0:100];
    plot_y=(-b(1)-b(2).*x)./b(3);
    plot_y_1{1,i}=plot_y;
        
Correct_green=matrixData{i,3}(:,1:3);
Correct_green(Correct_green(:,3)~=1,:)=[];
Correct_green=(Correct_green./200).*100;
Correct_green_1{1,i}=Correct_green;

Correct_green_HR=matrixData{i,3}(:,3:7);
Correct_green_HR(Correct_green_HR(:,1)~=1,:)=[];
Correct_green_HR=Correct_green_HR(:,2);
Correct_green_HR1{1,i}=Correct_green_HR;


Wrong_red=matrixData{i,3}(:,1:3);
Wrong_red(Wrong_red(:,3)~=0,:)=[];
Wrong_red=(Wrong_red./200).*100;
Wrong_red_1{1,i}=Wrong_red;

Wrong_red_HR=matrixData{i,3}(:,3:7);
Wrong_red_HR(Wrong_red_HR(:,1)~=0,:)=[];
Wrong_red_HR=Wrong_red_HR(:,2);
Wrong_red_HR1{1,i}=Wrong_red_HR;

        if subplotIndex > 64
            currentFigure = currentFigure + 1;
            subplotIndex = 1;
            figure(currentFigure);
        end
        subplot(8, 8, subplotIndex);
        subplotIndex = subplotIndex + 1;
 sz=300;
h(1)=scatter(Correct_green(:,1),Correct_green(:,2),sz,'+','g','LineWidth',2);
hold on
h(2)=scatter(Wrong_red(:,1),Wrong_red(:,2),sz,'s','r');
hold on
h(3)=plot(plot_x,plot_y,'-','LineWidth',4,'color','b');
hold on;
title(matrixData{i,2})
xlabel('Reward amount (%)','FontSize',10);
ylabel('Airpurff amount (%)','FontSize',10);
xlim([0 100]);
ylim([0 100]);



num_sessions=numel(Correct_green_1);
% Extend  grid boundaries to ensure coverage
x_edges = [2.5:2.5:100.5];
y_edges = [2.5:2.5:100.5];

% Initialize matrices to store counts for each session
num_sessions = length(Correct_green_1);
all_Correct_green_counts = zeros(length(y_edges)-1, length(x_edges)-1, num_sessions);
all_Wrong_red_counts = zeros(length(y_edges)-1, length(x_edges)-1, num_sessions);

% Loop for counting
for i = 1:num_sessions
    for x_idx = 1:length(x_edges)-1
        for y_idx = 1:length(y_edges)-1
            x_center = (x_edges(x_idx) + x_edges(x_idx+1))/2;
            y_center = (y_edges(y_idx) + y_edges(y_idx+1))/2;

            all_Correct_green_counts(y_idx, x_idx, i) = sum(sqrt((Correct_green_1{1,i}(:,1)-x_center).^2 + (Correct_green_1{1,i}(:,2)-y_center).^2) <= 2.5);
            all_Wrong_red_counts(y_idx, x_idx, i) = sum(sqrt((Wrong_red_1{1,i}(:,1)-x_center).^2 + (Wrong_red_1{1,i}(:,2)-y_center).^2) <= 2.5);
        end
    end
end

% Since we extended our grid, we'll remove the first and last rows/columns to get our original size back
% all_Correct_green_counts = all_Correct_green_counts(1:end-1, 1:end-1, :);
% all_Wrong_red_counts = all_Wrong_red_counts(1:end-1, 1:end-1, :);
end
% Average the count matrices across sessions
avg_Approach_All = mean(all_Correct_green_counts, 3);
avg_Avoidance_All = mean(all_Wrong_red_counts, 3);
DecisionDifference=avg_Approach_All-avg_Avoidance_All;


% Visualization
figure;
imagesc(x_edges, y_edges, DecisionDifference);  % Example: showing the difference between the two
colormap('jet');  % Or any other colormap of  choice
c = colorbar;
xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title(['Behavior'],'FontSize', 15, 'FontName', 'Arial'); 
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100,'FontSize', 18);
set(c, 'YTickLabel', {'Av', '', '', '', '', '', '', '', '', '', 'Ap'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11),'FontSize', 18);
box off;  % Removes the box around the figure
% Set the font size for the axes

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
%%% version b
% Visualization
figure;

% Apply a Gaussian smoothing filter
kernelSize = 5; % Size of the filter, can be adjusted
sigma = 1; % Standard deviation of the Gaussian kernel, adjust as needed
kernel = fspecial('gaussian', kernelSize, sigma); % Create a Gaussian kernel

% Smooth the data
smoothedData = imfilter(DecisionDifference, kernel, 'replicate');

% Display the smoothed heatmap
imagesc(x_edges, y_edges, smoothedData);  % Smoothed heatmap
colormap('jet');  % Or any other colormap of  choice
c = colorbar;

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Average across all sessions', 'FontSize', 15, 'FontName', 'Arial'); 
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100, 'FontSize', 18);
set(c, 'YTickLabel', {'Av', '', '', '', '', '', '', '', '', '', 'Ap'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11), 'FontSize', 18);
box off;  % Removes the box around the figure

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 5, 'Color', 'k');  % Average decision boundary
hold off;

%%% Approach and avoidance separately
% Visualization with two subplots
figure;

% Subplot for Approach Trials
subplot(1, 2, 1); % First subplot in a 1x2 grid
imagesc(x_edges, y_edges, avg_Approach_All);  % Displaying the Approach data
colormap('jet');  % Use a colormap of  choice
c=colorbar;
xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Approach Trials', 'FontSize', 17, 'FontName', 'Arial'); 
set(gca, 'YDir', 'normal');
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(c, 'YTickLabel', {'Low Frequency', '', '', '', '', '', '', '', '', '', 'High Frequency'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11));

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Subplot for Avoidance Trials
subplot(1, 2, 2); % Second subplot in a 1x2 grid
imagesc(x_edges, y_edges, avg_Avoidance_All);  % Displaying the Avoidance data
colormap('jet');
c=colorbar;
xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Avoidance Trials', 'FontSize', 17, 'FontName', 'Arial'); 
set(gca, 'YDir', 'normal');
xlim([0 100]);
ylim([0 100]);
set(c, 'YTickLabel', {'Low Frequency', '', '', '', '', '', '', '', '', '', 'High Frequency'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11));


% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Initialize matrices for heart rate sums and counts
heartRateSums = zeros(length(y_edges)-1, length(x_edges)-1);
heartRateCounts = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(matrixData, 1)
    % Extract heart rate data
    sessionData = matrixData{i, 3}; %  this is where  trial data is
    sessionHeartRates = matrixData{i, 3}(:,4); % Heart rate data

    % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage

    % Loop over each trial
    for j = 1:size(sessionData, 1)
        % Find the grid cell for each trial
        x_index = find(x_edges <= sessionData(j, 1), 1, 'last');
        y_index = find(y_edges <= sessionData(j, 2), 1, 'last');

        % Check if the trial falls within the grid boundaries and heart rate is within the specified range
        if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges) && sessionHeartRates(j) >= 100 && sessionHeartRates(j) <= 140
            % Add the heart rate to the sum for that cell and increment the count
            heartRateSums(y_index, x_index) = heartRateSums(y_index, x_index) + sessionHeartRates(j);
            heartRateCounts(y_index, x_index) = heartRateCounts(y_index, x_index) + 1;
        end
    end
end

% Calculate average heart rate for each grid cell
avgHeartRates = heartRateSums ./ heartRateCounts;
avgHeartRates(isnan(avgHeartRates)) = 0; % Replace NaN values with 0

% Visualization
figure;
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgHeartRates);
colormap('jet'); % Choose colormap as needed
colorbar;
caxis([100 130]); % Set the color scale to the heart rate range
xlabel('Reward amount (%)'); % Update with actual label
ylabel('Airpuff amount (%)'); % Update with actual label
title('Heart Rate Heatmap');
set(gca, 'YDir','normal'); % Set Y-axis to normal direction

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Additional plot formatting as needed
%%%%

% Initialize matrices for heart rate sums and counts for approach and avoidance decisions
approachHRSum = zeros(length(y_edges)-1, length(x_edges)-1);
approachHRCount = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRSum = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRCount = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(matrixData, 1)
    sessionData = matrixData{i, 3}; % Trial data
    sessionHeartRates = matrixData{i, 3}(:, 4); % Heart rate data
    decisions = sessionData(:, 3); % Decision data: 0 for avoidance, 1 for approach

    % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage

    % Loop over each trial
    for j = 1:size(sessionData, 1)
        % Only process heart rates within the specified range
        if sessionHeartRates(j) >= 100 && sessionHeartRates(j) <= 140
            x_index = find(x_edges <= sessionData(j, 1), 1, 'last');
            y_index = find(y_edges <= sessionData(j, 2), 1, 'last');

            if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges)
                if decisions(j) == 1 % Approach decision
                    approachHRSum(y_index, x_index) = approachHRSum(y_index, x_index) + sessionHeartRates(j);
                    approachHRCount(y_index, x_index) = approachHRCount(y_index, x_index) + 1;
                else % Avoidance decision
                    avoidanceHRSum(y_index, x_index) = avoidanceHRSum(y_index, x_index) + sessionHeartRates(j);
                    avoidanceHRCount(y_index, x_index) = avoidanceHRCount(y_index, x_index) + 1;
                end
            end
        end
    end
end

% Calculate average heart rates
avgApproachHR = approachHRSum ./ approachHRCount;
avgApproachHR(isnan(avgApproachHR)) = 0; % Replace NaN with 0
avgAvoidanceHR = avoidanceHRSum ./ avoidanceHRCount;
avgAvoidanceHR(isnan(avgAvoidanceHR)) = 0; % Replace NaN with 0

% Calculate the difference in heart rate
HRDifference = avgApproachHR - avgAvoidanceHR;

% Visualization for Difference in Heart Rate
figure;
imagesc(x_edges(1:end-1), y_edges(1:end-1), HRDifference);
colormap('jet');

% Create colorbar and set its properties
cb = colorbar;
caxis([-10 10]); % Set color scale range
set(cb, 'FontSize', 18); % Set font size for colorbar

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title(['Difference in Heart Rate (Approach - Avoidance)'], 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir', 'normal'); % Correct Y-axis orientation

% Set x and y ticks
xticks(0:10:100);
yticks(0:10:100);

% Set axis ticks' font size
set(gca, 'FontSize', 18);

% Decision boundary plot code
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Creating a figure with two subplots for Approach and Avoidance Trials
figure;

% Subplot for Approach Trials
subplot(1, 2, 1); % First subplot in a 1x2 grid
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgApproachHR);
colormap('jet');
colorbar;
caxis([100 130]); % Set color scale range
xlabel('Reward amount (%)'); % Update with actual label
ylabel('Airpuff amount (%)'); % Update with actual label
title('Average Heart Rate for Approach Trials');
set(gca, 'YDir','normal'); % Correct Y-axis orientation
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Subplot for Avoidance Trials
subplot(1, 2, 2); % Second subplot in a 1x2 grid
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgAvoidanceHR);
colormap('jet');
colorbar;
caxis([100 130]); % Set color scale range
xlabel('Reward amount (%)'); % Update with actual label
ylabel('Airpuff amount (%)'); % Update with actual label
title('Average Heart Rate for Avoidance Trials');
set(gca, 'YDir','normal'); % Correct Y-axis orientation
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Additional plot formatting as needed
% Variables to store the minimum values
minReward = Inf;
minPunishment = Inf;

% Loop through each session
for i = 1:size(matrixData, 1)
    sessionData = matrixData{i, 3}; % Trial data
        % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage
    % Check if sessionData is not empty
    if ~isempty(sessionData)
        % Update the minimum values
        
        currentMinReward = min(sessionData(:, 1));
        currentMinPunishment = min(sessionData(:, 2));

        if currentMinReward < minReward
            minReward = currentMinReward;
        end
        if currentMinPunishment < minPunishment
            minPunishment = currentMinPunishment;
        end
    end
end

% Display the minimum values
fprintf('Minimum Reward: %f\n', minReward);
fprintf('Minimum Punishment: %f\n', minPunishment);

%%%%%%%%%%%5
% Variables to store the sum of heart rates and the count of trials
totalHeartRate = 0;
countTrials = 0;

% Initialize total heart rate sum and trial count
totalHeartRate = 0;
countTrials = 0;

% Loop through each session
for i = 1:size(matrixData, 1)
    sessionData = matrixData{i, 3}; % Trial data
    sessionHeartRates = matrixData{i, 3}(:, 4); % Heart rate data
    decisions = sessionData(:, 3); % Decision data: 0 for avoidance, 1 for approach

    % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage

    % Loop over each trial
    for j = 1:size(sessionData, 1)
        rewardSize = sessionData(j, 1);
        airpuffSize = sessionData(j, 2);
        
        % Check if it's an approach trial, the sizes are within specified ranges, and heart rate is within the specified range
        if decisions(j) == 1 && rewardSize >= 2.5 && rewardSize <= 7.5 && airpuffSize >= 37.5 && airpuffSize <= 42.5 && sessionHeartRates(j) >= 100 && sessionHeartRates(j) <= 140
            totalHeartRate = totalHeartRate + sessionHeartRates(j);
            countTrials = countTrials + 1;
        end
    end
end

% Calculate the average heart rate
if countTrials > 0
    avgHeartRate = totalHeartRate / countTrials;
else
    avgHeartRate = NaN; % In case there are no such trials
end

% Display the average heart rate
fprintf('Average Heart Rate for approach trials with reward and airpuff size between 0-5: %f\n', avgHeartRate);

%%%%% for HRV total figure

% Initialize matrices for heart rate variability sums and counts for approach and avoidance decisions
approachHRVSum = zeros(length(y_edges)-1, length(x_edges)-1);
approachHRVCount = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRVSum = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRVCount = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(matrixData, 1)
    sessionData = matrixData{i, 3}; % Trial data
    sessionHRVs = matrixData{i, 3}(:, 6); % Heart rate variability data
    decisions = sessionData(:, 3); % Decision data: 0 for avoidance, 1 for approach

    % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage

    % Loop over each trial
    for j = 1:size(sessionData, 1)
        x_index = find(x_edges <= sessionData(j, 1), 1, 'last');
        y_index = find(y_edges <= sessionData(j, 2), 1, 'last');

        % Check if the trial falls within the grid boundaries and HRV is not greater than 80
        if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges) && sessionHRVs(j) <= 800
            if decisions(j) == 1 % Approach decision
                approachHRVSum(y_index, x_index) = approachHRVSum(y_index, x_index) + sessionHRVs(j);
                approachHRVCount(y_index, x_index) = approachHRVCount(y_index, x_index) + 1;
            else % Avoidance decision
                avoidanceHRVSum(y_index, x_index) = avoidanceHRVSum(y_index, x_index) + sessionHRVs(j);
                avoidanceHRVCount(y_index, x_index) = avoidanceHRVCount(y_index, x_index) + 1;
            end
        end
    end
end

% Calculate average heart rate variability
avgApproachHRV = approachHRVSum ./ approachHRVCount;
avgApproachHRV(isnan(avgApproachHRV)) = 0; % Replace NaN with 0
avgAvoidanceHRV = avoidanceHRVSum ./ avoidanceHRVCount;
avgAvoidanceHRV(isnan(avgAvoidanceHRV)) = 0; % Replace NaN with 0

% Calculate the difference in heart rate variability
HRVDifference = avgApproachHRV - avgAvoidanceHRV;

% Visualization for Difference in Heart Rate variability
figure;
imagesc(x_edges(1:end-1), y_edges(1:end-1), HRVDifference);
colormap('jet');
colorbar;
caxis([-40 40]); % Set color scale range
xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Difference in HRV (Approach - Avoidance)', 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir','normal'); % Correct Y-axis orientation
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Set x and y ticks
xticks(0:10:100);
yticks(0:10:100);
% Set axis ticks' font size
set(gca, 'FontSize', 18);

% Smoothed Visualization for Difference in HRV
figure;

% Apply a Gaussian smoothing filter
kernelSize = 5; % Size of the filter, can be adjusted
sigma = 1; % Standard deviation of the Gaussian kernel, adjust as needed
kernel = fspecial('gaussian', kernelSize, sigma); % Create a Gaussian kernel

% Smooth the data
smoothedHRVDifference = imfilter(HRVDifference, kernel, 'replicate');

% Display the smoothed heatmap
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedHRVDifference);  % Smoothed heatmap
colormap('jet');  % Or any other colormap of  choice
c = colorbar;
caxis([-40 40]); % Set color scale range

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Smoothed Difference in HRV (Approach - Avoidance)', 'FontSize', 15, 'FontName', 'Arial'); 
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100, 'FontSize', 18);
% set(c, 'YTickLabel', {'Av', '', '', '', '', '', '', '', '', '', 'Ap'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11), 'FontSize', 18);
box off;  % Removes the box around the figure

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 5, 'Color', 'k');  % Average decision boundary
hold off;

%%%%%%%%%%%%%%5
% Initialize matrices for HRV sums and counts for approach and avoidance decisions
approachHRVSum = zeros(length(y_edges)-1, length(x_edges)-1);
approachHRVCount = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRVSum = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRVCount = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(matrixData, 1)
    sessionData = matrixData{i, 3}; % Trial data

    % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage

    sessionHRV = matrixData{i, 3}(:, 6); % HRV data
    decisions = sessionData(:, 3); % Decision data: 0 for avoidance, 1 for approach

    % Loop over each trial
    for j = 1:size(sessionData, 1)
        rewardSize = sessionData(j, 1); % Reward size
        airpuffSize = sessionData(j, 2); % Airpuff size
        x_index = find(x_edges <= rewardSize, 1, 'last');
        y_index = find(y_edges <= airpuffSize, 1, 'last');

        % Check if the trial falls within the grid boundaries and HRV is not greater than 800
        if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges) && sessionHRV(j) <= 800
            if decisions(j) == 1 % Approach decision
                approachHRVSum(y_index, x_index) = approachHRVSum(y_index, x_index) + sessionHRV(j);
                approachHRVCount(y_index, x_index) = approachHRVCount(y_index, x_index) + 1;
            else % Avoidance decision
                avoidanceHRVSum(y_index, x_index) = avoidanceHRVSum(y_index, x_index) + sessionHRV(j);
                avoidanceHRVCount(y_index, x_index) = avoidanceHRVCount(y_index, x_index) + 1;
            end
        end
    end
end

% Calculate average HRV for Approach and Avoidance
avgApproachHRV = approachHRVSum ./ approachHRVCount;
avgApproachHRV(isnan(avgApproachHRV)) = 0; % Replace NaN with 0
avgAvoidanceHRV = avoidanceHRVSum ./ avoidanceHRVCount;
avgAvoidanceHRV(isnan(avgAvoidanceHRV)) = 0; % Replace NaN with 0
% Define a Gaussian smoothing kernel
kernelSize = 3; % Size of the kernel (3x3, 5x5, etc.)
sigma = 1; % Standard deviation of the Gaussian
[x, y] = meshgrid(-floor(kernelSize/2):floor(kernelSize/2), -floor(kernelSize/2):floor(kernelSize/2));
gaussianKernel = exp(-(x.^2 + y.^2) / (2*sigma^2));
gaussianKernel = gaussianKernel / sum(gaussianKernel(:)); % Normalize the kernel
% Smooth the data
smoothedApproachHRV = conv2(avgApproachHRV, gaussianKernel, 'same');
smoothedAvoidanceHRV = conv2(avgAvoidanceHRV, gaussianKernel, 'same');

% Creating a figure with two subplots for Approach and Avoidance HRV
figure;

% Subplot for Approach HRV
subplot(1, 2, 1); % First subplot in a 1x2 grid
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedApproachHRV);
colormap('jet');
caxis([500 700]); % Set color scale from 500 to 800
colorbar;
xlabel('Reward Percentage', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff Percentage', 'FontSize', 18, 'FontName', 'Arial'); 
title('Average HRV for Approach Trials', 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir','normal'); % Correct Y-axis orientation
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Subplot for Avoidance HRV
subplot(1, 2, 2); % Second subplot in a 1x2 grid
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedAvoidanceHRV);
colormap('jet');
caxis([500 700]); % Set color scale from 500 to 800
colorbar;
xlabel('Reward Percentage', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff Percentage', 'FontSize', 18, 'FontName', 'Arial'); 
title('Average HRV for Avoidance Trials', 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir','normal'); % Correct Y-axis orientation
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Creating a figure with two subplots for Approach and Avoidance HRV
figure;

% Subplot for Approach HRV
subplot(1, 2, 1); % First subplot in a 1x2 grid
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgApproachHRV);
colormap('jet');
caxis([500 700]); % Set color scale from 500 to 800
colorbar;
xlabel('Reward Percentage'); % Update with actual label
ylabel('Airpuff Percentage'); % Update with actual label
title('Average HRV for Approach Trials');
set(gca, 'YDir','normal'); % Correct Y-axis orientation
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Subplot for Avoidance HRV
subplot(1, 2, 2); % Second subplot in a 1x2 grid
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgAvoidanceHRV);
colormap('jet');
caxis([500 700]); % Set color scale from 500 to 800
colorbar;
xlabel('Reward Percentage'); % Update with actual label
ylabel('Airpuff Percentage'); % Update with actual label
title('Average HRV for Avoidance Trials');
set(gca, 'YDir','normal'); % Correct Y-axis orientation
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Variables to store the maximum values
maxReward = -Inf;
maxPunishment = -Inf;

% Loop through each session
for i = 1:size(matrixData, 1)
    sessionData = matrixData{i, 3}; % Trial data
    % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage

    % Check if sessionData is not empty
    if ~isempty(sessionData)
        % Update the maximum values
        currentMaxReward = max(sessionData(:, 1));
        currentMaxPunishment = max(sessionData(:, 2));

        if currentMaxReward > maxReward
            maxReward = currentMaxReward;
        end
        if currentMaxPunishment > maxPunishment
            maxPunishment = currentMaxPunishment;
        end
    end
end

% Convert to percentage
maxRewardPercentage = maxReward;
maxPunishmentPercentage = maxPunishment;

% Display the maximum values
fprintf('Maximum Reward Percentage: %f%%\n', maxRewardPercentage);
fprintf('Maximum Punishment Percentage: %f%%\n', maxPunishmentPercentage);

% Initialize empty cell arrays for Debbie and Prez with the same number of columns as matrixData
data_Debbie = cell(0, size(matrixData, 2));
data_Prez = cell(0, size(matrixData, 2));

% Loop through matrixData to separate data for Debbie and Prez
for i = 1:size(matrixData, 1)
    if strcmp(matrixData{i, 1}, 'Debbie')
        data_Debbie(end + 1, :) = matrixData(i, :);
    elseif strcmp(matrixData{i, 1}, 'Prez')
        data_Prez(end + 1, :) = matrixData(i, :);
    end
end

% Function to convert reward and punishment sizes to percentages
convertToPercentage = @(data) data / 2;

% Initialize matrices for heart rate sums and counts for approach and avoidance decisions
approachHRSum_Debbie = zeros(length(y_edges)-1, length(x_edges)-1);
approachHRCount_Debbie = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRSum_Debbie = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRCount_Debbie = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(data_Debbie, 1)
    sessionData_Debbie = data_Debbie{i, 3}; % Trial data
    sessionHeartRates_Debbie = data_Debbie{i, 3}(:, 4); % Heart rate data
    decisions = sessionData_Debbie(:, 3); % Decision data: 0 for avoidance, 1 for approach

    % Convert reward and punishment values to percentages
    sessionData_Debbie(:, 1) = sessionData_Debbie(:, 1) / 2; % Reward percentage
    sessionData_Debbie(:, 2) = sessionData_Debbie(:, 2) / 2; % Punishment percentage

    % Loop over each trial
    for j = 1:size(sessionData_Debbie, 1)
        x_index = find(x_edges <= sessionData_Debbie(j, 1), 1, 'last');
        y_index = find(y_edges <= sessionData_Debbie(j, 2), 1, 'last');

        % Check if the trial falls within the grid boundaries and heart rate is within the specified range
        if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges) && sessionHeartRates_Debbie(j) >= 100 && sessionHeartRates_Debbie(j) <= 140
            if decisions(j) == 1 % Approach decision
                approachHRSum_Debbie(y_index, x_index) = approachHRSum_Debbie(y_index, x_index) + sessionHeartRates_Debbie(j);
                approachHRCount_Debbie(y_index, x_index) = approachHRCount_Debbie(y_index, x_index) + 1;
            else % Avoidance decision
                avoidanceHRSum_Debbie(y_index, x_index) = avoidanceHRSum_Debbie(y_index, x_index) + sessionHeartRates_Debbie(j);
                avoidanceHRCount_Debbie(y_index, x_index) = avoidanceHRCount_Debbie(y_index, x_index) + 1;
            end
        end
    end
end

% Calculate average heart rates
avgApproachHR_Debbie = approachHRSum_Debbie ./ approachHRCount_Debbie;
avgApproachHR_Debbie(isnan(avgApproachHR_Debbie)) = 0; % Replace NaN with 0
avgAvoidanceHR_Debbie = avoidanceHRSum_Debbie ./ avoidanceHRCount_Debbie;
avgAvoidanceHR_Debbie(isnan(avgAvoidanceHR_Debbie)) = 0; % Replace NaN with 0

%% for Prez
% Initialize matrices for heart rate sums and counts for approach and avoidance decisions
approachHRSum_Prez = zeros(length(y_edges)-1, length(x_edges)-1);
approachHRCount_Prez = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRSum_Prez = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceHRCount_Prez = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(data_Prez, 1)
    sessionData_Prez = data_Prez{i, 3}; % Trial data
    sessionHeartRates_Prez = data_Prez{i, 3}(:, 4); % Heart rate data
    decisions = sessionData_Prez(:, 3); % Decision data: 0 for avoidance, 1 for approach

    % Convert reward and punishment values to percentages
    sessionData_Prez(:, 1) = sessionData_Prez(:, 1) / 2; % Reward percentage
    sessionData_Prez(:, 2) = sessionData_Prez(:, 2) / 2; % Punishment percentage

    % Loop over each trial
    for j = 1:size(sessionData_Prez, 1)
        x_index = find(x_edges <= sessionData_Prez(j, 1), 1, 'last');
        y_index = find(y_edges <= sessionData_Prez(j, 2), 1, 'last');

        % Check if the trial falls within the grid boundaries and heart rate is within the specified range
        if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges) && sessionHeartRates_Prez(j) >= 100 && sessionHeartRates_Prez(j) <= 140
            if decisions(j) == 1 % Approach decision
                approachHRSum_Prez(y_index, x_index) = approachHRSum_Prez(y_index, x_index) + sessionHeartRates_Prez(j);
                approachHRCount_Prez(y_index, x_index) = approachHRCount_Prez(y_index, x_index) + 1;
            else % Avoidance decision
                avoidanceHRSum_Prez(y_index, x_index) = avoidanceHRSum_Prez(y_index, x_index) + sessionHeartRates_Prez(j);
                avoidanceHRCount_Prez(y_index, x_index) = avoidanceHRCount_Prez(y_index, x_index) + 1;
            end
        end
    end
end

% Calculate average heart rates
avgApproachHR_Prez = approachHRSum_Prez ./ approachHRCount_Prez;
avgApproachHR_Prez(isnan(avgApproachHR_Prez)) = 0; % Replace NaN with 0
avgAvoidanceHR_Prez = avoidanceHRSum_Prez ./ avoidanceHRCount_Prez;
avgAvoidanceHR_Prez(isnan(avgAvoidanceHR_Prez)) = 0; % Replace NaN with 0

% Calculate the difference in heart rate
HRDifference = avgApproachHR_Prez - avgAvoidanceHR_Prez;

% Creating a figure with four subplots
figure;

% Subplot 1: Approach Trials for Debbie
subplot(2, 2, 1);
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgApproachHR_Debbie);
caxis([100 130]); % Set color scale range
colormap('jet');
colorbar;
title('Debbie: Average Heart Rate for Approach Trials');
xlabel('Reward Size (%)');
ylabel('Punishment Size (%)');
set(gca, 'YDir','normal');
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Subplot 2: Avoidance Trials for Debbie
subplot(2, 2, 2);
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgApproachHR_Debbie);
caxis([100 130]); % Set color scale range
colormap('jet');
colorbar;
title('Debbie: Average Heart Rate for Avoidance Trials');
xlabel('Reward Size (%)');
ylabel('Punishment Size (%)');
set(gca, 'YDir','normal');
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Subplot 3: Approach Trials for Prez
subplot(2, 2, 3);
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgApproachHR_Prez);
caxis([100 130]); % Set color scale range
colormap('jet');
colorbar;
title('Prez: Average Heart Rate for Approach Trials');
xlabel('Reward Size (%)');
ylabel('Punishment Size (%)');
set(gca, 'YDir','normal');
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Subplot 4: Avoidance Trials for Prez
subplot(2, 2, 4);
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgAvoidanceHR_Prez);
caxis([100 130]); % Set color scale range
colormap('jet');
colorbar;
title('Prez: Average Heart Rate for Avoidance Trials');
xlabel('Reward Size (%)');
ylabel('Punishment Size (%)');
set(gca, 'YDir','normal');
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

%%%%%%
%%% for Reaction Times

% Initialize matrices for Reaction Time sums and counts for approach and avoidance decisions
approachRTSum_Debbie = zeros(length(y_edges)-1, length(x_edges)-1);
approachRTCount_Debbie = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceRTSum_Debbie = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceRTCount_Debbie = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(data_Debbie, 1)
    sessionData_Debbie = data_Debbie{i, 3}; % Trial data
    sessionReactionTimes_Debbie = data_Debbie{i, 3}(:,8); % Reaction Time data
    decisions = sessionData_Debbie(:, 3); % Decision data: 0 for avoidance, 1 for approach

        % Convert reward and punishment values to percentages
    sessionData_Debbie(:, 1) = sessionData_Debbie(:, 1) / 2; % Reward percentage
    sessionData_Debbie(:, 2) = sessionData_Debbie(:, 2) / 2; % Punishment percentage
    sessionData_Debbie(:, 8) = sessionData_Debbie(:, 8) *1e3 ; % Reaction time in ms


    % Loop over each trial
    for j = 1:size(sessionData_Debbie, 1)
        x_index = find(x_edges <= sessionData_Debbie(j, 1), 1, 'last');
        y_index = find(y_edges <= sessionData_Debbie(j, 2), 1, 'last');

        if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges)
            if decisions(j) == 1 % Approach decision
                approachRTSum_Debbie(y_index, x_index) = approachRTSum_Debbie(y_index, x_index) + sessionReactionTimes_Debbie(j);
                approachRTCount_Debbie(y_index, x_index) = approachRTCount_Debbie(y_index, x_index) + 1;
            else % Avoidance decision
                avoidanceRTSum_Debbie(y_index, x_index) = avoidanceRTSum_Debbie(y_index, x_index) + sessionReactionTimes_Debbie(j);
                avoidanceRTCount_Debbie(y_index, x_index) = avoidanceRTCount_Debbie(y_index, x_index) + 1;
            end
        end
    end
end

% Calculate average Reaction Times
avgApproachRT_Debbie = approachRTSum_Debbie ./ approachRTCount_Debbie;
avgApproachRT_Debbie(isnan(avgApproachRT_Debbie)) = 0; % Replace NaN with 0
avgAvoidanceRT_Debbie = avoidanceRTSum_Debbie ./ avoidanceRTCount_Debbie;
avgAvoidanceRT_Debbie(isnan(avgAvoidanceRT_Debbie)) = 0; % Replace NaN with 0

RTDifference_Debbie = avgApproachRT_Debbie - avgAvoidanceRT_Debbie;

%  plot code
figure;
imagesc(x_edges(1:end-1), y_edges(1:end-1), RTDifference_Debbie);
colormap('jet');

% Create colorbar and set its properties
cb = colorbar;
caxis([-0.06 0.06]); % Set color scale range
set(cb, 'FontSize', 18); % Set font size for colorbar

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title(['Difference in Reaction Time (Approach - Avoidance) - Little Debbie'], 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir', 'normal'); % Correct Y-axis orientation

% Set x and y ticks
xticks(0:10:100);
yticks(0:10:100);

% Set axis ticks' font size
set(gca, 'FontSize', 18);

% Decision boundary plot code
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Smoothed Visualization for Difference in HRV
figure;

% Apply a Gaussian smoothing filter
kernelSize = 5; % Size of the filter, can be adjusted
sigma = 1; % Standard deviation of the Gaussian kernel, adjust as needed
kernel = fspecial('gaussian', kernelSize, sigma); % Create a Gaussian kernel

% Smooth the data
smoothedHRVDifference = imfilter(HRVDifference, kernel, 'replicate');

% Display the smoothed heatmap
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedHRVDifference);  % Smoothed heatmap
colormap('jet');  % Or any other colormap of  choice
c = colorbar;
caxis([-40 40]); % Set color scale range

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Smoothed Difference in HRV (Approach - Avoidance)', 'FontSize', 15, 'FontName', 'Arial'); 
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100, 'FontSize', 18);
set(c, 'YTickLabel', {'Av', '', '', '', '', '', '', '', '', '', 'Ap'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11), 'FontSize', 18);
box off;  % Removes the box around the figure

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 5, 'Color', 'k');  % Average decision boundary
hold off;

%% for Prez
% Initialize matrices for Reaction Time sums and counts for approach and avoidance decisions
approachRTSum_Prez = zeros(length(y_edges)-1, length(x_edges)-1);
approachRTCount_Prez = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceRTSum_Prez = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceRTCount_Prez = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(data_Prez, 1)
    sessionData_Prez = data_Prez{i, 3}; % Trial data
    sessionReactionTimes_Prez = data_Prez{i, 3}(:,8); % Reaction Time data
    decisions = sessionData_Prez(:, 3); % Decision data: 0 for avoidance, 1 for approach

        % Convert reward and punishment values to percentages
    sessionData_Prez(:, 1) = sessionData_Prez(:, 1) / 2; % Reward percentage
    sessionData_Prez(:, 2) = sessionData_Prez(:, 2) / 2; % Punishment percentage
    sessionData_Prez(:, 8) = sessionData_Prez(:, 8) *1e3 ; % Reaction time in ms


    % Loop over each trial
    for j = 1:size(sessionData_Prez, 1)
        x_index = find(x_edges <= sessionData_Prez(j, 1), 1, 'last');
        y_index = find(y_edges <= sessionData_Prez(j, 2), 1, 'last');

        if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges)
            if decisions(j) == 1 % Approach decision
                approachRTSum_Prez(y_index, x_index) = approachRTSum_Prez(y_index, x_index) + sessionReactionTimes_Prez(j);
                approachRTCount_Prez(y_index, x_index) = approachRTCount_Prez(y_index, x_index) + 1;
            else % Avoidance decision
                avoidanceRTSum_Prez(y_index, x_index) = avoidanceRTSum_Prez(y_index, x_index) + sessionReactionTimes_Prez(j);
                avoidanceRTCount_Prez(y_index, x_index) = avoidanceRTCount_Prez(y_index, x_index) + 1;
            end
        end
    end
end

% Calculate average Reaction Times
avgApproachRT_Prez = approachRTSum_Prez ./ approachRTCount_Prez;
avgApproachRT_Prez(isnan(avgApproachRT_Prez)) = 0; % Replace NaN with 0
avgAvoidanceRT_Prez = avoidanceRTSum_Prez ./ avoidanceRTCount_Prez;
avgAvoidanceRT_Prez(isnan(avgAvoidanceRT_Prez)) = 0; % Replace NaN with 0

% Calculate the difference in Reaction Time
RTDifference_Prez = avgApproachRT_Prez - avgAvoidanceRT_Prez;

%  plot code
figure;
imagesc(x_edges(1:end-1), y_edges(1:end-1), RTDifference_Prez);
colormap('jet');

% Create colorbar and set its properties
cb = colorbar;
caxis([-0.1 0.1]); % Set color scale range
set(cb, 'FontSize', 18); % Set font size for colorbar

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title(['Difference in Reaction Time (Approach - Avoidance) - Prez'], 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir', 'normal'); % Correct Y-axis orientation

% Set x and y ticks
xticks(0:10:100);
yticks(0:10:100);

% Set axis ticks' font size
set(gca, 'FontSize', 18);

% Decision boundary plot code
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Creating a figure with four subplots
figure;

% Subplot 1: Approach Trials for Debbie
subplot(2, 2, 1);
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgApproachRT_Debbie*1e3);
caxis([300 600]); % Set color scale range
colormap('jet');
colorbar;
title('Debbie: Average Reaction Time for Approach Trials');
xlabel('Reward Size (%)');
ylabel('Punishment Size (%)');
set(gca, 'YDir','normal');
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Subplot 2: Avoidance Trials for Debbie
subplot(2, 2, 2);
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgApproachRT_Debbie*1e3);
caxis([300 600]); % Set color scale range
colormap('jet');
colorbar;
title('Debbie: Average Reaction Time for Avoidance Trials');
xlabel('Reward Size (%)');
ylabel('Punishment Size (%)');
set(gca, 'YDir','normal');
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Subplot 3: Approach Trials for Prez
subplot(2, 2, 3);
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgApproachRT_Prez*1e3);
caxis([300 600]); % Set color scale range
colormap('jet');
colorbar;
title('Prez: Average Reaction Time for Approach Trials');
xlabel('Reward Size (%)');
ylabel('Punishment Size (%)');
set(gca, 'YDir','normal');
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Subplot 4: Avoidance Trials for Prez
subplot(2, 2, 4);
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgAvoidanceRT_Prez*1e3);
caxis([300 600]); % Set color scale range
colormap('jet');
colorbar;
title('Prez: Average Reaction Time  for Avoidance Trials');
xlabel('Reward Size (%)');
ylabel('Punishment Size (%)');
set(gca, 'YDir','normal');
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;
% Initialize matrices for ReactionTime sums and counts for approach and avoidance decisions
approachRTSum = zeros(length(y_edges)-1, length(x_edges)-1);
approachRTCount = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceRTSum = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceRTCount = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop tRTough each session
for i = 1:size(matrixData, 1)
    sessionData = matrixData{i, 3}; % Trial data
    sessionReactionTimes = matrixData{i, 3}(:,8); % ReactionTime data
    decisions = sessionData(:, 3); % Decision data: 0 for avoidance, 1 for approach

        % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage


    % Loop over each trial
    for j = 1:size(sessionData, 1)
        x_index = find(x_edges <= sessionData(j, 1), 1, 'last');
        y_index = find(y_edges <= sessionData(j, 2), 1, 'last');

        if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges)
            if decisions(j) == 1 % Approach decision
                approachRTSum(y_index, x_index) = approachRTSum(y_index, x_index) + sessionReactionTimes(j);
                approachRTCount(y_index, x_index) = approachRTCount(y_index, x_index) + 1;
            else % Avoidance decision
                avoidanceRTSum(y_index, x_index) = avoidanceRTSum(y_index, x_index) + sessionReactionTimes(j);
                avoidanceRTCount(y_index, x_index) = avoidanceRTCount(y_index, x_index) + 1;
            end
        end
    end
end

% Calculate average ReactionTimes
avgApproachRT = approachRTSum ./ approachRTCount;
avgApproachRT(isnan(avgApproachRT)) = 0; % Replace NaN with 0
avgAvoidanceRT = avoidanceRTSum ./ avoidanceRTCount;
avgAvoidanceRT(isnan(avgAvoidanceRT)) = 0; % Replace NaN with 0

% Calculate the difference in ReactionTime
RTDifference = (avgApproachRT - avgAvoidanceRT)*1e3;

%  plot code
figure;
imagesc(x_edges(1:end-1), y_edges(1:end-1), RTDifference);
colormap('jet');

% Create colorbar and set its properties
cb = colorbar;
caxis([-50 50]); % Set color scale range
set(cb, 'FontSize', 18); % Set font size for colorbar

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title(['Difference in Reaction Time (Approach - Avoidance)'], 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir', 'normal'); % Correct Y-axis orientation

% Set x and y ticks
xticks(0:10:100);
yticks(0:10:100);

% Set axis ticks' font size
set(gca, 'FontSize', 18);

% Decision boundary plot code
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 5, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Visualization of Reaction Time Differences with Decision Boundary
figure;

% Apply a Gaussian smoothing filter
kernelSize = 5; % Size of the filter, can be adjusted
sigma = 1; % Standard deviation of the Gaussian kernel, adjust as needed
kernel = fspecial('gaussian', kernelSize, sigma); % Create a Gaussian kernel

% Smooth the data
smoothedData = imfilter(RTDifference, kernel, 'replicate');  % Smooth the reaction time difference data

% Display the smoothed heatmap
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedData);  % Display the smoothed heatmap
colormap('jet');

% Create colorbar and set its properties
cb = colorbar;
caxis([-100 100]); % Set color scale range
set(cb, 'FontSize', 18); % Set font size for colorbarset(cb, 'FontSize', 18); % Set font size for colorbar
% Set the x and y axis properties according to the original specs
xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Reaction Time (Ap minus Av)', 'FontSize', 15, 'FontName', 'Arial'); 
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100, 'FontSize', 20);
% set(c, 'YTickLabel', {'Av', '', '', '', '', '', '', '', '', '', 'Ap'}, 'YTick', linspace(c.Limits(1), c.Limits(2), 11), 'FontSize', 20);
box off;  % Removes the box around the figure

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  betaValues is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 5, 'Color', [0 0 0]);  % Plot the average decision boundary
hold off;

% Creating a figure with two subplots for Approach and Avoidance Trials
figure;

% Subplot for Approach Trials
subplot(1, 2, 1); % First subplot in a 1x2 grid
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgApproachRT*1e3);
colormap('jet');
colorbar;
% caxis([100 130]); % Set color scale range
xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Average RT for Approach Trials', 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir','normal'); % Correct Y-axis orientation
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Subplot for Avoidance Trials
subplot(1, 2, 2); % Second subplot in a 1x2 grid
imagesc(x_edges(1:end-1), y_edges(1:end-1), avgAvoidanceRT*1e3);
colormap('jet');
colorbar;
% caxis([100 130]); % Set color scale range
xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Average RT for Avoidance Trials', 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir','normal'); % Correct Y-axis orientation
% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;


% for lick data
% Initialize matrices for Lick rate sums and counts for approach and avoidance decisions
approachLickSum = zeros(length(y_edges)-1, length(x_edges)-1);
approachLickCount = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceLickSum = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceLickCount = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(matrixData, 1)
    sessionData = matrixData{i, 3}; % Trial data
    sessionLickRates = matrixData{i, 3}(:, 9); % Lick rate data
    decisions = sessionData(:, 3); % Decision data: 0 for avoidance, 1 for approach

    % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage

    % Loop over each trial
    for j = 1:size(sessionData, 1)
        % Only process Lick rates within the specified range
%         if sessionLickRates(j) >= 100 && sessionLickRates(j) <= 140
            x_index = find(x_edges <= sessionData(j, 1), 1, 'last');
            y_index = find(y_edges <= sessionData(j, 2), 1, 'last');

            if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges)
                if decisions(j) == 1 % Approach decision
                    approachLickSum(y_index, x_index) = approachLickSum(y_index, x_index) + sessionLickRates(j);
                    approachLickCount(y_index, x_index) = approachLickCount(y_index, x_index) + 1;
                else % Avoidance decision
                    avoidanceLickSum(y_index, x_index) = avoidanceLickSum(y_index, x_index) + sessionLickRates(j);
                    avoidanceLickCount(y_index, x_index) = avoidanceLickCount(y_index, x_index) + 1;
                end
            end
%         end
    end
end

% Calculate average Lick rates
avgapproachLick = approachLickSum ./ approachLickCount;
avgapproachLick(isnan(avgapproachLick)) = 0; % Replace NaN with 0
avgAvoidanceLick = avoidanceLickSum ./ avoidanceLickCount;
avgAvoidanceLick(isnan(avgAvoidanceLick)) = 0; % Replace NaN with 0

% Calculate the difference in Lick rate
LickDifference = (avgapproachLick - avgAvoidanceLick);

% Visualization for Difference in Lick rate
figure;
imagesc(x_edges(1:end-1), y_edges(1:end-1), LickDifference);
colormap('jet');

% Create colorbar and set its properties
cb = colorbar;
caxis([-0.3e-04 0.3e-04]); % Set color scale range
set(cb, 'FontSize', 18); % Set font size for colorbar

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title(['Difference in Lick rate (Approach - Avoidance)'], 'FontSize', 15, 'FontName', 'Arial'); 
set(gca, 'YDir', 'normal'); % Correct Y-axis orientation

% Set x and y ticks
xticks(0:10:100);
yticks(0:10:100);

% Set axis ticks' font size
set(gca, 'FontSize', 18);

% Decision boundary plot code
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', [0 0 0]);  % Average decision boundary
hold off;

% Smoothed Visualization for Difference in Lick Rate
figure;

% Apply a Gaussian smoothing filter
kernelSize = 5; % Size of the filter, can be adjusted
sigma = 1; % Standard deviation of the Gaussian kernel, adjust as needed
kernel = fspecial('gaussian', kernelSize, sigma); % Create a Gaussian kernel

% Smooth the data
smoothedLickDifference = imfilter(LickDifference, kernel, 'replicate');

% Display the smoothed heatmap
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedLickDifference);  % Smoothed heatmap
colormap('jet');  % Or any other colormap of  choice

% Create colorbar and set its properties
cb = colorbar;
caxis([-0.3e-04 0.3e-04]); % Set color scale range
set(cb, 'FontSize', 18); % Set font size for colorbar

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Smoothed Difference in Lick Rate (Approach - Avoidance)', 'FontSize', 15, 'FontName', 'Arial'); 
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100, 'FontSize', 18);
box off;  % Removes the box around the figure

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 5, 'Color', 'k');  % Average decision boundary
hold off;

% For lick data
% Initialize matrices for Lick rate sums and counts for approach and avoidance decisions
approachLickSum = zeros(length(y_edges)-1, length(x_edges)-1);
approachLickCount = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceLickSum = zeros(length(y_edges)-1, length(x_edges)-1);
avoidanceLickCount = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:size(matrixData, 1)
    sessionData = matrixData{i, 3}; % Trial data
    sessionLickRates = matrixData{i, 3}(:, 9); % Lick rate data
    decisions = sessionData(:, 3); % Decision data: 0 for avoidance, 1 for approach

    % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage

    % Loop over each trial
    for j = 1:size(sessionData, 1)
        x_index = find(x_edges <= sessionData(j, 1), 1, 'last');
        y_index = find(y_edges <= sessionData(j, 2), 1, 'last');

        if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges)
            if decisions(j) == 1 % Approach decision
                approachLickSum(y_index, x_index) = approachLickSum(y_index, x_index) + sessionLickRates(j);
                approachLickCount(y_index, x_index) = approachLickCount(y_index, x_index) + 1;
            else % Avoidance decision
                avoidanceLickSum(y_index, x_index) = avoidanceLickSum(y_index, x_index) + sessionLickRates(j);
                avoidanceLickCount(y_index, x_index) = avoidanceLickCount(y_index, x_index) + 1;
            end
        end
    end
end

% Calculate average Lick rates
avgApproachLick = approachLickSum ./ approachLickCount;
avgApproachLick(isnan(avgApproachLick)) = 0; % Replace NaN with 0
avgAvoidanceLick = avoidanceLickSum ./ avoidanceLickCount;
avgAvoidanceLick(isnan(avgAvoidanceLick)) = 0; % Replace NaN with 0

% Apply a Gaussian smoothing filter
kernelSize = 5; % Size of the filter, can be adjusted
sigma = 1; % Standard deviation of the Gaussian kernel, adjust as needed
kernel = fspecial('gaussian', kernelSize, sigma); % Create a Gaussian kernel

% Smooth the average Lick rates
smoothedApproachLick = imfilter(avgApproachLick, kernel, 'replicate');
smoothedAvoidanceLick = imfilter(avgAvoidanceLick, kernel, 'replicate');

% Determine common color scale limits for consistency
minValue = min([smoothedApproachLick(:); smoothedAvoidanceLick(:)]);
maxValue = max([smoothedApproachLick(:); smoothedAvoidanceLick(:)]);

% Visualization: Create a figure with two subplots using tiledlayout
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.4]); % Adjust the figure size as needed

% Use tiledlayout for better control over subplot sizes
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

% First tile: Approach trials
ax1 = nexttile;
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedApproachLick);
colormap(ax1, 'jet');
cb1 = colorbar;
caxis([minValue maxValue]); % Use common color scale limits
set(cb1, 'FontSize', 18);

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial');
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial');
title('Smoothed Lick Rate - Approach Trials', 'FontSize', 15, 'FontName', 'Arial');
set(ax1, 'YDir', 'normal'); % Correct Y-axis orientation
set(ax1, 'FontSize', 18);
xticks(0:10:100);
yticks(0:10:100);
box off;

% Ensure the axes are square
axis equal;
ax1.PlotBoxAspectRatio = [1 1 1];

% Plot the decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2); %  betaValues is defined
avg_plot_y = (-avg_b(1) - avg_b(2) .* plot_x) ./ avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', 'k');
hold off;

% Second tile: Avoidance trials
ax2 = nexttile;
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedAvoidanceLick);
colormap(ax2, 'jet');
cb2 = colorbar;
caxis([minValue maxValue]); % Use common color scale limits
set(cb2, 'FontSize', 18);

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial');
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial');
title('Smoothed Lick Rate - Avoidance Trials', 'FontSize', 15, 'FontName', 'Arial');
set(ax2, 'YDir', 'normal'); % Correct Y-axis orientation
set(ax2, 'FontSize', 18);
xticks(0:10:100);
yticks(0:10:100);
box off;

% Ensure the axes are square
axis equal;
ax2.PlotBoxAspectRatio = [1 1 1];

% Plot the decision boundary
hold on;
plot_x = 0:100;
avg_plot_y = (-avg_b(1) - avg_b(2) .* plot_x) ./ avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 3, 'Color', 'k');
hold off;

% Flatten the matrices into vectors
HRDifference = HRDifference(:);
RTDifference = RTDifference(:);
DecisionDifference = DecisionDifference(:);
HRVDifference = HRVDifference(:); % Flattening HRVDifference
LickDifference = LickDifference(:); % Flattening LickDifference

% Calculate the correlation coefficients and p-values
[rHR_RT, pValueHR_RT] = corrcoef(HRDifference, RTDifference);
[rHR_Decision, pValueHR_Decision] = corrcoef(HRDifference, DecisionDifference);
[rRT_Decision, pValueRT_Decision] = corrcoef(RTDifference, DecisionDifference);
[rHRV_HR, pValueHRV_HR] = corrcoef(HRVDifference, HRDifference);
[rHRV_RT, pValueHRV_RT] = corrcoef(HRVDifference, RTDifference);
[rHRV_Decision, pValueHRV_Decision] = corrcoef(HRVDifference, DecisionDifference);
[rLick_HR, pValueLick_HR] = corrcoef(LickDifference, HRDifference);
[rLick_RT, pValueLick_RT] = corrcoef(LickDifference, RTDifference);
[rLick_Decision, pValueLick_Decision] = corrcoef(LickDifference, DecisionDifference);
[rLick_HRV, pValueLick_HRV] = corrcoef(LickDifference, HRVDifference);

% Extract the correlation coefficients and p-values
correlationHR_RT = rHR_RT(1,2);
correlationHR_Decision = rHR_Decision(1,2);
correlationRT_Decision = rRT_Decision(1,2);
correlationHRV_HR = rHRV_HR(1,2);
correlationHRV_RT = rHRV_RT(1,2);
correlationHRV_Decision = rHRV_Decision(1,2);
correlationLick_HR = rLick_HR(1,2);
correlationLick_RT = rLick_RT(1,2);
correlationLick_Decision = rLick_Decision(1,2);
correlationLick_HRV = rLick_HRV(1,2);

pValueHR_RT = pValueHR_RT(1,2);
pValueHR_Decision = pValueHR_Decision(1,2);
pValueRT_Decision = pValueRT_Decision(1,2);
pValueHRV_HR = pValueHRV_HR(1,2);
pValueHRV_RT = pValueHRV_RT(1,2);
pValueHRV_Decision = pValueHRV_Decision(1,2);
pValueLick_HR = pValueLick_HR(1,2);
pValueLick_RT = pValueLick_RT(1,2);
pValueLick_Decision = pValueLick_Decision(1,2);
pValueLick_HRV = pValueLick_HRV(1,2);


% with 8 decimals
% Display the results for HR, RT, Decision, HRV, and Lick correlations
fprintf('Correlation between HRDifference and RTDifference: %.5f (p-value: %.5f)\n', correlationHR_RT, pValueHR_RT);
fprintf('Correlation between HRDifference and DecisionDifference: %.5f (p-value: %.5f)\n', correlationHR_Decision, pValueHR_Decision);
fprintf('Correlation between RTDifference and DecisionDifference: %.5f (p-value: %.5f)\n', correlationRT_Decision, pValueRT_Decision);
fprintf('Correlation between HRVDifference and HRDifference: %.5f (p-value: %.5f)\n', correlationHRV_HR, pValueHRV_HR);
fprintf('Correlation between HRVDifference and RTDifference: %.5f (p-value: %.5f)\n', correlationHRV_RT, pValueHRV_RT);
fprintf('Correlation between HRVDifference and DecisionDifference: %.5f (p-value: %.5f)\n', correlationHRV_Decision, pValueHRV_Decision);
fprintf('Correlation between LickDifference and HRDifference: %.5f (p-value: %.5f)\n', correlationLick_HR, pValueLick_HR);
fprintf('Correlation between LickDifference and RTDifference: %.5f (p-value: %.5f)\n', correlationLick_RT, pValueLick_RT);
fprintf('Correlation between LickDifference and DecisionDifference: %.5f (p-value: %.5f)\n', correlationLick_Decision, pValueLick_Decision);
fprintf('Correlation between LickDifference and HRVDifference: %.5f (p-value: %.5f)\n', correlationLick_HRV, pValueLick_HRV);

% Define the trial sizes of interest
trialSizes = [30, 60, 90, 120, 150, 180];

% Initialize matrices to hold session-level data 
numSessions = size(matrixData, 1); %  matrixData is a cell array with session data
sessionMeansReward = NaN(numSessions, length(trialSizes));
sessionMeansAirpuff = NaN(numSessions, length(trialSizes));

% Loop through each session to calculate mean lick rates
for i = 1:numSessions
    pavlovianData = matrixData{i, 4}; % Extract session data
    for j = 1:length(trialSizes)
        trialSize = trialSizes(j);
        % Extract rewards and airpuff data for the current trial size
        rewardData = pavlovianData(pavlovianData(:,1) == trialSize & pavlovianData(:,2) == 1, 9);
        airpuffData = pavlovianData(pavlovianData(:,1) == trialSize & pavlovianData(:,2) == 0, 9);
        % Calculate means if data is present
        if ~isempty(rewardData)
            sessionMeansReward(i, j) = mean(rewardData, 'omitnan');
        end
        if ~isempty(airpuffData)
            sessionMeansAirpuff(i, j) = mean(airpuffData, 'omitnan');
        end
    end
end

% Calculate overall means across sessions
overallAvgReward = nanmean(sessionMeansReward);
overallAvgAirpuff = nanmean(sessionMeansAirpuff);

% Calculate SEM for each condition across trial sizes
semReward = nanstd(sessionMeansReward) ./ sqrt(sum(~isnan(sessionMeansReward)));
semAirpuff = nanstd(sessionMeansAirpuff) ./ sqrt(sum(~isnan(sessionMeansAirpuff)));

% Define colors for the plots
rewardColor = [0, 0, 1]; % Blue color for Reward
airpuffColor = [1, 0, 0]; % Red color for Airpuff

% Calculate the upper and lower bounds for the shaded areas
upperBoundReward = overallAvgReward + semReward;
lowerBoundReward = overallAvgReward - semReward;
upperBoundAirpuff = overallAvgAirpuff + semAirpuff;
lowerBoundAirpuff = overallAvgAirpuff - semAirpuff;

% Plotting with shaded error bars
figure; hold on;

rewardLine = plot(trialSizes, overallAvgReward,'-o', 'Color', rewardColor, 'LineWidth', 2, 'DisplayName', 'Reward');
airpuffLine = plot(trialSizes, overallAvgAirpuff,'-o', 'Color', airpuffColor, 'LineWidth', 2, 'DisplayName', 'Airpuff');

% Finalize the plot
xlabel('Outcome Size/Duration', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Lick activity (Hz)', 'FontSize', 18, 'FontName', 'Arial'); 
title(['Average lick activity during Pavlovian trials'],'FontSize', 14, 'FontName', 'Arial'); 

% Creating the legend and capturing the legend handle
lgd = legend([rewardLine, airpuffLine], {'Reward', 'Airpuff'}, 'Location', 'best');

% Setting the font size of the legend
lgd.FontSize = 15; % 
% Set axis ticks' font size
set(gca, 'FontSize', 18);
lgd.Box = 'off';
grid on;

% Define the trial sizes of interest
trialSizes = [30, 60, 90, 120, 150, 180];

% Initialize matrices to hold session-level data
numSessions = size(matrixData, 1); %  matrixData is a cell array with session data
sessionMeansHR_Reward = NaN(numSessions, length(trialSizes));
sessionMeansHR_Airpuff = NaN(numSessions, length(trialSizes));

% Loop through each session to calculate mean HR
for i = 1:numSessions
    pavlovianData = matrixData{i, 4}; % Extract session data
    for j = 1:length(trialSizes)
        trialSize = trialSizes(j);
        % Extract HR data for reward and airpuff for the current trial size
        rewardHRData = pavlovianData(pavlovianData(:,1) == trialSize & pavlovianData(:,2) == 1 & pavlovianData(:,6) >= 1 & pavlovianData(:,6) <= 800, 6);
        airpuffHRData = pavlovianData(pavlovianData(:,1) == trialSize & pavlovianData(:,2) == 0 & pavlovianData(:,6) >= 1 & pavlovianData(:,6) <= 800, 6);
        % Calculate means for this session and trial size, if data is present
        if ~isempty(rewardHRData)
            sessionMeansHR_Reward(i, j) = mean(rewardHRData, 'omitnan');
        end
        if ~isempty(airpuffHRData)
            sessionMeansHR_Airpuff(i, j) = mean(airpuffHRData, 'omitnan');
        end
    end
end

% Calculate overall means across sessions
overallAvgHR_Reward = nanmean(sessionMeansHR_Reward, 1);  % Mean across sessions (row-wise)
overallAvgHR_Airpuff = nanmean(sessionMeansHR_Airpuff, 1);  % Mean across sessions (row-wise)

% Calculate SEM for each condition across trial sizes
semHR_Reward = nanstd(sessionMeansHR_Reward, [], 1) ./ sqrt(sum(~isnan(sessionMeansHR_Reward), 1));
semHR_Airpuff = nanstd(sessionMeansHR_Airpuff, [], 1) ./ sqrt(sum(~isnan(sessionMeansHR_Airpuff), 1));

% Define colors for the plots
rewardColor = [0, 0, 1]; % Blue color for Reward
airpuffColor = [1, 0, 0]; % Red color for Airpuff

% Calculate the upper and lower bounds for the shaded areas
upperBoundHR_Reward = overallAvgHR_Reward + semHR_Reward;
lowerBoundHR_Reward = overallAvgHR_Reward - semHR_Reward;
upperBoundHR_Airpuff = overallAvgHR_Airpuff + semHR_Airpuff;
lowerBoundHR_Airpuff = overallAvgHR_Airpuff - semHR_Airpuff;

% Plotting with shaded error bars
figure; hold on;

% Reward - shaded error bar
fill([trialSizes, fliplr(trialSizes)], [lowerBoundHR_Reward, fliplr(upperBoundHR_Reward)], ...
    rewardColor, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
rewardLine = plot(trialSizes, overallAvgHR_Reward, 'Color', rewardColor, 'LineWidth', 2, 'DisplayName', 'Reward');

% Airpuff - shaded error bar
fill([trialSizes, fliplr(trialSizes)], [lowerBoundHR_Airpuff, fliplr(upperBoundHR_Airpuff)], ...
    airpuffColor, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
airpuffLine = plot(trialSizes, overallAvgHR_Airpuff, 'Color', airpuffColor, 'LineWidth', 2, 'DisplayName', 'Airpuff');

% Finalize the plot
xlabel('Outcome Size/Duration', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Average Heart Rate', 'FontSize', 18, 'FontName', 'Arial'); 
title(['Average Heart Rate for Different Trial Sizes'],'FontSize', 14, 'FontName', 'Arial'); 

% Creating the legend and capturing the legend handle
lgd = legend([rewardLine, airpuffLine], {'Reward', 'Airpuff'}, 'Location', 'best');

% Setting the font size of the legend
lgd.FontSize = 15;
lgd.Box = 'off';

% Set axis ticks' font size
set(gca, 'FontSize', 18);
grid on;

% for HRV
% Define the trial sizes of interest
trialSizes = [30, 60, 90, 120, 150, 180];

% Initialize matrices to hold session-level data
numSessions = size(matrixData, 1); %  matrixData is a cell array with session data
sessionMeansHRV_Reward = NaN(numSessions, length(trialSizes));
sessionMeansHRV_Airpuff = NaN(numSessions, length(trialSizes));

% Loop through each session to calculate mean HRV for each trial size
for i = 1:numSessions
    pavlovianData = matrixData{i, 4}; % Extract Pavlovian data for session i
    
    for j = 1:length(trialSizes)
        trialSize = trialSizes(j);
        
        % Filter data for the current trial size for reward and airpuff
        rewardData = pavlovianData(pavlovianData(:,1) == trialSize & pavlovianData(:,2) == 1  & pavlovianData(:,6) <= 800, 6);
        airpuffData = pavlovianData(pavlovianData(:,1) == trialSize & pavlovianData(:,2) == 0  & pavlovianData(:,6) <= 800, 6);
        
        % Calculate the mean HRV for this session and trial size, if data is present
        if ~isempty(rewardData)
            sessionMeansHRV_Reward(i, j) = mean(rewardData, 'omitnan');
        end
        if ~isempty(airpuffData)
            sessionMeansHRV_Airpuff(i, j) = mean(airpuffData, 'omitnan');
        end
    end
end

% Calculate the overall mean across sessions for each trial size
meanHRV_Reward_TrialSizes = nanmean(sessionMeansHRV_Reward, 1);
meanHRV_Airpuff_TrialSizes = nanmean(sessionMeansHRV_Airpuff, 1);

% Now calculate the overall mean across all levels (trial sizes)
overallMeanHRV_Reward = mean(meanHRV_Reward_TrialSizes, 'omitnan');
overallMeanHRV_Airpuff = mean(meanHRV_Airpuff_TrialSizes, 'omitnan');

% Plotting
figure;
plot(trialSizes, meanHRV_Reward_TrialSizes, '-o', 'Color', rewardColor, 'LineWidth', 2, 'DisplayName', 'Reward');
hold on;
plot(trialSizes, meanHRV_Airpuff_TrialSizes, '-o', 'Color', airpuffColor, 'LineWidth', 2, 'DisplayName', 'Airpuff');
xlabel('Outcome Size/Duration');
ylabel('Mean HRV (ms)');
title('Average Heart Rate Variability during Pavlovian trials', 'FontSize', 18, 'FontName', 'Arial');
legend('show');
grid on;

% Display overall means across all levels (for reference)
disp(['Overall Mean HRV for Reward: ', num2str(overallMeanHRV_Reward)]);
disp(['Overall Mean HRV for Airpuff: ', num2str(overallMeanHRV_Airpuff)]);

% Set axis ticks' font size
set(gca, 'FontSize', 18);
grid on;

% Setting the font size of the legend
lgd.FontSize = 15;
lgd.Box = 'off';


% Initialize the index for the loop
i = 1;

% Loop through the cell array
while i <= size(trialdata, 1)
    % Check if the size of the element in the 4th column is 1
    if size(trialdata{i,4}, 1) < 3
        % Remove the entire row
        trialdata(i,:) = [];
    else
        % Move to the next row if the row is not removed
        i = i + 1;
    end
end

trialData=trialdata;
% Define the trial sizes of interest
trialSizes = [30, 60, 90, 120, 150, 180];

% Initialize matrices to hold session-level data
numSessions = size(trialData, 1); %  trialData is a cell array with session data
sessionMeansPupil_Reward = NaN(numSessions, length(trialSizes));
sessionMeansPupil_Airpuff = NaN(numSessions, length(trialSizes));

% Loop through each session to calculate mean pupil size for each trial size
for i = 1:numSessions
    sessionData = trialData{i, 4}; % Extract the data for this session
    
    for j = 1:length(trialSizes)
        trialSize = trialSizes(j);
        
        % Filter data for the current trial size for reward and airpuff
        rewardPupilData = sessionData(sessionData(:,1) == trialSize & sessionData(:,2) == 1, 3);
        airpuffPupilData = sessionData(sessionData(:,1) == trialSize & sessionData(:,2) == 0, 3);
        
        % Calculate the mean pupil size for this session and trial size, if data is present
        if ~isempty(rewardPupilData)
            sessionMeansPupil_Reward(i, j) = mean(rewardPupilData, 'omitnan');
        end
        if ~isempty(airpuffPupilData)
            sessionMeansPupil_Airpuff(i, j) = mean(airpuffPupilData, 'omitnan');
        end
    end
end

% Calculate the overall mean across sessions for each trial size
meanPupil_Reward_TrialSizes = nanmean(sessionMeansPupil_Reward, 1);
meanPupil_Airpuff_TrialSizes = nanmean(sessionMeansPupil_Airpuff, 1);

% Now calculate the overall mean across all levels (trial sizes)
overallMeanPupil_Reward = mean(meanPupil_Reward_TrialSizes, 'omitnan');
overallMeanPupil_Airpuff = mean(meanPupil_Airpuff_TrialSizes, 'omitnan');

% Plotting
figure;

% Plot data for Reward and Airpuff
plot(trialSizes, meanPupil_Reward_TrialSizes, '-o', 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Reward');
hold on;
plot(trialSizes, meanPupil_Airpuff_TrialSizes, '-o', 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Airpuff');

% Labels and title
xlabel('Outcome Size/Duration', 'FontSize', 15, 'FontName', 'Arial');
ylabel('Average Pupil Diameter (z-score)', 'FontSize', 15, 'FontName', 'Arial');
title('Average Pupil Diameter for Different Trial Sizes', 'FontSize', 15, 'FontName', 'Arial');

% Set x-axis range and ticks
xlim([0 200]);
xticks([0 50 100 150 200]);

% Adjust the y-axis format to remove leading zeros (e.g., .01 instead of 0.01)
ax = gca; % Get current axes
ax.YAxis.Exponent = 0; % Disable scientific notation
ax.YTickLabel = strrep(ax.YTickLabel, '0.', '.'); % Modify y-tick labels to remove leading zeros

% Set x and y axis tick label font size to 4
ax.FontSize = 14; % Applies to both x and y axis values
ax.FontName = 'Arial';

% Customize the legend
legend('show', 'FontSize', 14, 'FontName', 'Arial', 'Box', 'off'); % Remove box from the legend

% Enable grid
grid on;
% Display overall means across all levels (for reference)
disp(['Overall Mean Pupil Diameter for Reward: ', num2str(overallMeanPupil_Reward)]);
disp(['Overall Mean Pupil Diameter for Airpuff: ', num2str(overallMeanPupil_Airpuff)]);

% Initialize matrices for pupil size sums and counts for approach and avoidance decisions
approachPupilSum = zeros(length(y_edges)-1, length(x_edges)-1);
approachPupilCount = zeros(length(y_edges)-1, length(x_edges)-1);
avoidancePupilSum = zeros(length(y_edges)-1, length(x_edges)-1);
avoidancePupilCount = zeros(length(y_edges)-1, length(x_edges)-1);

% Loop through each session
for i = 1:length(trialData)
    sessionData = trialData{i, 3}; % Trial data
    sessionPupilSizes = trialData{i, 3}(:, 4); % Pupil size data
    decisions = sessionData(:, 3); % Decision data: 0 for avoidance, 1 for approach

    % Convert reward and punishment values to percentages
    sessionData(:, 1) = sessionData(:, 1) / 2; % Reward percentage
    sessionData(:, 2) = sessionData(:, 2) / 2; % Punishment percentage

    % Loop over each trial
    for j = 1:length(sessionData)
        % Only process pupil sizes within a specified range
        % if sessionPupilSizes(j) >= some_lower_limit && sessionPupilSizes(j) <= some_upper_limit  % Define limits
            x_index = find(x_edges <= sessionData(j, 1), 1, 'last');
            y_index = find(y_edges <= sessionData(j, 2), 1, 'last');

            if ~isempty(x_index) && ~isempty(y_index) && x_index < length(x_edges) && y_index < length(y_edges)
                if decisions(j) == 1 % Approach decision
                    approachPupilSum(y_index, x_index) = approachPupilSum(y_index, x_index) + sessionPupilSizes(j);
                    approachPupilCount(y_index, x_index) = approachPupilCount(y_index, x_index) + 1;
                else % Avoidance decision
                    avoidancePupilSum(y_index, x_index) = avoidancePupilSum(y_index, x_index) + sessionPupilSizes(j);
                    avoidancePupilCount(y_index, x_index) = avoidancePupilCount(y_index, x_index) + 1;
                end
            end
        % end
    end
end

% Calculate average pupil sizes
avgApproachPupil = approachPupilSum ./ approachPupilCount;
avgApproachPupil(isnan(avgApproachPupil)) = 0; % Replace NaN with 0
avgAvoidancePupil = avoidancePupilSum ./ avoidancePupilCount;
avgAvoidancePupil(isnan(avgAvoidancePupil)) = 0; % Replace NaN with 0

% Calculate the difference in pupil size
PupilDifference = avgApproachPupil - avgAvoidancePupil;

% Smoothed Visualization for Difference in Pupil Size
figure;

% Apply a Gaussian smoothing filter
kernelSize = 5; % Size of the filter, can be adjusted
sigma = 1; % Standard deviation of the Gaussian kernel, adjust as needed
kernel = fspecial('gaussian', kernelSize, sigma); % Create a Gaussian kernel

% Smooth the data
smoothedPupilDifference = imfilter(PupilDifference, kernel, 'replicate');

% Display the smoothed heatmap
imagesc(x_edges(1:end-1), y_edges(1:end-1), smoothedPupilDifference);  % Smoothed heatmap
colormap('jet');  % Or any other colormap of  choice

% Create colorbar and set its properties
cb = colorbar;
caxis([-0.02 0.02]); % Set color scale range
set(cb, 'FontSize', 18); % Set font size for colorbar

xlabel('Reward amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
ylabel('Airpuff amount (%)', 'FontSize', 18, 'FontName', 'Arial'); 
title('Smoothed Difference in Pupil Size (Approach - Avoidance)', 'FontSize', 15, 'FontName', 'Arial'); 
xlim([0 100]);  % Setting x-axis limits
ylim([0 100]);  % Setting y-axis limits
set(gca, 'YDir', 'normal');  % Adjust y-axis direction
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', 'XGrid', 'off', 'XTick', 0:10:100, 'YTick', 0:10:100, 'FontSize', 18);
% set(cb, 'YTickLabel', {'Av', '', '', '', '', '', '', '', '', '', 'Ap'}, 'YTick', linspace(cb.Limits(1), cb.Limits(2), 11), 'FontSize', 18);
box off;  % Removes the box around the figure

% Decision boundary
hold on;
plot_x = 0:100;
avg_b = mean(cell2mat(betaValues), 2);  %  b_1 is a cell array containing 'b' vectors for all sessions
avg_plot_y = (-avg_b(1) - avg_b(2).*plot_x)./avg_b(3);
plot(plot_x, avg_plot_y, '-', 'LineWidth', 5, 'Color', 'k');  % Average decision boundary
hold off;

% PupilDifference, RTDifference, DecisionDifference, HRVDifference, LickDifference

% Calculate the correlation coefficients and p-values for Pupil data
[rPupil_RT, pValuePupil_RT] = corrcoef(PupilDifference, RTDifference);
[rPupil_Decision, pValuePupil_Decision] = corrcoef(PupilDifference, DecisionDifference);
[rRT_Decision, pValueRT_Decision] = corrcoef(RTDifference, DecisionDifference);
[rHRV_Pupil, pValueHRV_Pupil] = corrcoef(HRVDifference, PupilDifference);
[rHRV_RT, pValueHRV_RT] = corrcoef(HRVDifference, RTDifference);
[rHRV_Decision, pValueHRV_Decision] = corrcoef(HRVDifference, DecisionDifference);
[rLick_Pupil, pValueLick_Pupil] = corrcoef(LickDifference, PupilDifference);
[rLick_RT, pValueLick_RT] = corrcoef(LickDifference, RTDifference);
[rLick_Decision, pValueLick_Decision] = corrcoef(LickDifference, DecisionDifference);
[rLick_HRV, pValueLick_HRV] = corrcoef(LickDifference, HRVDifference);

% Extract the correlation coefficients and p-values
correlationPupil_RT = rPupil_RT(1,2);
correlationPupil_Decision = rPupil_Decision(1,2);
correlationRT_Decision = rRT_Decision(1,2);
correlationHRV_Pupil = rHRV_Pupil(1,2);
correlationHRV_RT = rHRV_RT(1,2);
correlationHRV_Decision = rHRV_Decision(1,2);
correlationLick_Pupil = rLick_Pupil(1,2);
correlationLick_RT = rLick_RT(1,2);
correlationLick_Decision = rLick_Decision(1,2);
correlationLick_HRV = rLick_HRV(1,2);

pValuePupil_RT = pValuePupil_RT(1,2);
pValuePupil_Decision = pValuePupil_Decision(1,2);
pValueRT_Decision = pValueRT_Decision(1,2);
pValueHRV_Pupil = pValueHRV_Pupil(1,2);
pValueHRV_RT = pValueHRV_RT(1,2);
pValueHRV_Decision = pValueHRV_Decision(1,2);
pValueLick_Pupil = pValueLick_Pupil(1,2);
pValueLick_RT = pValueLick_RT(1,2);
pValueLick_Decision = pValueLick_Decision(1,2);
pValueLick_HRV = pValueLick_HRV(1,2);

% Display the results for Pupil, RT, Decision, HRV, and Lick correlations
fprintf('Correlation between PupilDifference and RTDifference: %.5f(p-value: %.5f)\n', correlationPupil_RT, pValuePupil_RT);
fprintf('Correlation between PupilDifference and DecisionDifference: %.5f(p-value: %.5f)\n', correlationPupil_Decision, pValuePupil_Decision);
fprintf('Correlation between RTDifference and DecisionDifference: %.5f(p-value: %.5f)\n', correlationRT_Decision, pValueRT_Decision);
fprintf('Correlation between HRVDifference and PupilDifference: %.5f(p-value: %.5f)\n', correlationHRV_Pupil, pValueHRV_Pupil);
fprintf('Correlation between HRVDifference and RTDifference: %.5f(p-value: %.5f)\n', correlationHRV_RT, pValueHRV_RT);
fprintf('Correlation between HRVDifference and DecisionDifference: %.5f(p-value: %.5f)\n', correlationHRV_Decision, pValueHRV_Decision);
fprintf('Correlation between LickDifference and PupilDifference: %.5f (p-value: %.5f)\n', correlationLick_Pupil, pValueLick_Pupil);
fprintf('Correlation between LickDifference and RTDifference: %.5f (p-value: %.5f)\n', correlationLick_RT, pValueLick_RT);
fprintf('Correlation between LickDifference and DecisionDifference: %.5f (p-value: %.5f)\n', correlationLick_Decision, pValueLick_Decision);
fprintf('Correlation between LickDifference and HR Difference: %.5f (p-value: %.5f)\n', correlationLick_HRV, pValueLick_HRV);

% Correlation between PupilDifference and RTDifference
[correlationPupil_RT, pValuePupil_RT] = computeCorrelation(PupilDifference, RTDifference);

% Correlation between PupilDifference and DecisionDifference
[correlationPupil_Decision, pValuePupil_Decision] = computeCorrelation(PupilDifference, DecisionDifference);

% Correlation between RTDifference and DecisionDifference
[correlationRT_Decision, pValueRT_Decision] = computeCorrelation(RTDifference, DecisionDifference);

% Correlation between HRVDifference and PupilDifference
[correlationHRV_Pupil, pValueHRV_Pupil] = computeCorrelation(HRVDifference, PupilDifference);

% Correlation between HRVDifference and RTDifference
[correlationHRV_RT, pValueHRV_RT] = computeCorrelation(HRVDifference, RTDifference);

% Correlation between HRVDifference and DecisionDifference
[correlationHRV_Decision, pValueHRV_Decision] = computeCorrelation(HRVDifference, DecisionDifference);

% Correlation between LickDifference and PupilDifference
[correlationLick_Pupil, pValueLick_Pupil] = computeCorrelation(LickDifference, PupilDifference);

% Correlation between LickDifference and RTDifference
[correlationLick_RT, pValueLick_RT] = computeCorrelation(LickDifference, RTDifference);

% Correlation between LickDifference and DecisionDifference
[correlationLick_Decision, pValueLick_Decision] = computeCorrelation(LickDifference, DecisionDifference);

% Correlation between LickDifference and HRVDifference
[correlationLick_HRV, pValueLick_HRV] = computeCorrelation(LickDifference, HRVDifference);

% Display the results
fprintf('Correlation between PupilDifference and RTDifference: %.5f (p-value: %.5f)\n', correlationPupil_RT, pValuePupil_RT);
fprintf('Correlation between PupilDifference and DecisionDifference: %.5f (p-value: %.5f)\n', correlationPupil_Decision, pValuePupil_Decision);
fprintf('Correlation between RTDifference and DecisionDifference: %.5f (p-value: %.5f)\n', correlationRT_Decision, pValueRT_Decision);
fprintf('Correlation between HRVDifference and PupilDifference: %.5f (p-value: %.5f)\n', correlationHRV_Pupil, pValueHRV_Pupil);
fprintf('Correlation between HRVDifference and RTDifference: %.5f (p-value: %.5f)\n', correlationHRV_RT, pValueHRV_RT);
fprintf('Correlation between HRVDifference and DecisionDifference: %.5f (p-value: %.5f)\n', correlationHRV_Decision, pValueHRV_Decision);
fprintf('Correlation between LickDifference and PupilDifference: %.5f (p-value: %.5f)\n', correlationLick_Pupil, pValueLick_Pupil);
fprintf('Correlation between LickDifference and RTDifference: %.5f (p-value: %.5f)\n', correlationLick_RT, pValueLick_RT);
fprintf('Correlation between LickDifference and DecisionDifference: %.5f (p-value: %.5f)\n', correlationLick_Decision, pValueLick_Decision);
fprintf('Correlation between LickDifference and HRVDifference: %.5f (p-value: %.5f)\n', correlationLick_HRV, pValueLick_HRV);

% Function to compute correlation with data cleaning
function [correlationCoefficient, pValue] = computeCorrelation(var1, var2)
    % Flatten the variables into vectors (if they are matrices)
    var1Vector = var1(:);
    var2Vector = var2(:);

    % Remove NaN and zero values from both variables
    validIndices = ~isnan(var1Vector) & ~isnan(var2Vector) & var1Vector ~= 0 & var2Vector ~= 0;

    % Extract the valid data
    var1Valid = var1Vector(validIndices);
    var2Valid = var2Vector(validIndices);

    % Check if there are enough data points
    if length(var1Valid) >= 2
        % Compute the correlation coefficient and p-value
        [rMatrix, pValueMatrix] = corrcoef(var1Valid, var2Valid);
        correlationCoefficient = rMatrix(1, 2);
        pValue = pValueMatrix(1, 2);
    else
        % Not enough data to compute correlation
        correlationCoefficient = NaN;
        pValue = NaN;
        warning('Not enough valid data points to compute correlation between variables.');
    end
end


