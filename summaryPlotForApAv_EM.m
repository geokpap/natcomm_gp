% Updated data for cOFC and pACC from the provided table
cOFC = [34, 41, 90, 41, 35, 39, 29, 26, 64, 24, 68, 16, 26];
pACC = [29, 37, 88, 56, 16, 33, 29, 15, 55, 22, 66, 12, 27];

% Updated event labels
events = {'Fix point ON BOTH', 'Choice cue ON noEM', 'Choice cue ON EM', 'Ap Air ON BOTH', ...
          'Av rew ON BOTH', 'Ap highlight BOTH', 'Av highlight BOTH', 'Red bar ON noEM', ...
          'Red bar ON EM', 'Yellow bar ON noEM', 'Yellow bar ON EM', 'Red bar Rew ON BOTH', 'Yellow bar Air ON BOTH'};

% Set up the figure size for a landscape style
figure('Units', 'inches', 'Position', [0.5, 0.5, 14, 8]); % 14x8 inches figure

% Calculate positions and width for the bars to accommodate two groups
numEvents = length(events);
positions = 1:numEvents;
width = 0.3; % Width for two groups

% Define subtle colors for cOFC and pACC
colors = [0.7, 0.7, 0.9; 0.8, 0.8, 0.7]; % Subtle colors for two data sets

hold on;

% Create the bar plots for cOFC and pACC
bar(positions - width/2, cOFC, width, 'FaceColor', colors(1,:));
bar(positions + width/2, pACC, width, 'FaceColor', colors(2,:));

hold off;

% Customize the plot
set(gca, 'XTick', positions, 'XTickLabel', events, 'XTickLabelRotation', 45);
ylabel('Neurons (%)', 'FontSize', 20, 'FontWeight', 'bold');
title('Neural Activity by Event and Brain Region during the Ap-Av EM task', 'FontSize', 22, 'FontWeight', 'bold');
legend({'cOFC', 'pACC'}, 'Location', 'bestoutside', 'Box', 'off', 'FontSize', 18);
set(gca, 'FontSize', 18, 'FontName', 'Arial');
set(gcf, 'Color', 'w'); % Set background color to white
box off;

% Ensure the plot is displayed appealingly
axis tight;
grid on;
