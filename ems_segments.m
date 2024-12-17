% Combined Script
% Initialize
clc;
close all;
clear;

% Define the prefixes and dates for both scripts
names = {'deb', 'deb1', 'prez', 'prez1'};

% Dates for the first script
deb_dates1 = {'20191104', '20200428', '20200508', ...
              '20200522', '20201215', '20210819', '20210820', ...
              '20211231', '20220112', '20220126', '20220128', '20220227', ...
              '20220303', '20220520'};
prez_dates1 = {'20190627', '20190721', '20190726', '20190806', '20190808', ...
               '20190821', '20190823', '20190826', '20190906', '20190922', ...
               '20191003', '20191009', '20191029'};

% Dates for the second script
deb_dates2 = {'20191116', '20191130', '20200401', '20200404', '20200407', ...
              '20200408', '20200411', '20200415', '20200416', '20200417', ...
              '20200506', '20200507', '20200514', '20200515', '20200520', ...
              '20200521', '20200527', '20200528', '20200529', '20200603', ...
              '20200604', '20200610', '20200611', '20200612', '20200617', ...
              '20200619', '20200624', '20201216', '20201224', '20201229', ...
              '20210106', '20210114', '20210727', '20210728', '20210729', ...
              '20210730', '20210802', '20210803', '20210804', '20210805', ...
              '20210813', '20210817', '20210818', '20210823', '20210825', ...
              '20211229', '20211230', '20220111', '20220122', '20220125', ...
              '20220127', '20220202', '20220302', '20220326', '20220519', ...
              '20220520'};
prez_dates2 = {'20190531', '20190601', '20190812', '20190814', '20190819', ...
               '20190829', '20190901', '20190903', '20190911', '20190915', ...
               '20190919', '20190922', '20191001', '20191007', '20191011', ...
               '20191103'};

% Process both sets of data
[relative_change_avoidance_all1, relative_change_approach_all1] = process_dates(deb_dates1, prez_dates1);
[relative_change_avoidance_all2, relative_change_approach_all2] = process_dates(deb_dates2, prez_dates2);

% Plot the results
plot_combined_results(relative_change_avoidance_all1, relative_change_approach_all1, relative_change_avoidance_all2, relative_change_approach_all2);

%% Function Definitions

% Main function to process a list of dates
function [relative_change_avoidance_all, relative_change_approach_all] = process_dates(deb_dates, prez_dates)
    relative_change_avoidance_all = [];
    relative_change_approach_all = [];
    
    % Combine dates into a single list and process each date
    combined_dates = [deb_dates, prez_dates];
    for date_str = combined_dates
        date_str = date_str{1};  % Convert cell to string
        try
            % Load data and process blocks
            [Included_trials_1, Included_trials_2] = load_and_process_blocks(date_str, deb_dates, prez_dates);
            if isempty(Included_trials_1) || isempty(Included_trials_2)
                continue;
            end
            
            % Calculate and accumulate relative changes
            [relative_change_avoidance, relative_change_approach] = calculate_relative_changes(Included_trials_1, Included_trials_2, 5);
            relative_change_avoidance_all = [relative_change_avoidance_all; relative_change_avoidance];
            relative_change_approach_all = [relative_change_approach_all; relative_change_approach];
        catch ME
            fprintf('Error processing date %s: %s\n', date_str, ME.message);
        end
    end
end

% Function to load and process data blocks
function [Included_trials_1, Included_trials_2] = load_and_process_blocks(date_str, deb_dates, prez_dates)
    prefix_list = get_prefix_list(date_str, deb_dates, prez_dates);
    name = get_existing_prefix(prefix_list, date_str);
    
    % Load data
    A = load_data(name, date_str, 'BarSize');
    B = load_data(name, date_str, 'TrialType');
    
    % Determine block sizes and split indices
    blockSizes = determine_block_sizes(length(A));
    blockIndices = mat2cell(1:length(A), 1, blockSizes);
    
    % Process blocks
    Included_trials_1 = process_block(A, B, blockIndices{1});
    Included_trials_2 = process_block(A, B, blockIndices{2});
end

% Function to plot the results from both data sets
function plot_combined_results(relative_change_avoidance_all1, relative_change_approach_all1, relative_change_avoidance_all2, relative_change_approach_all2)
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 0.5]);
    
    % First subplot
    subplot(1, 2, 1);
    plot_stacked_relative_changes(relative_change_avoidance_all1, relative_change_approach_all1, 'EMS (p < 0.05)', relative_change_avoidance_all1, relative_change_avoidance_all2, relative_change_approach_all1, relative_change_approach_all2, true);
    axis square;

    % Second subplot
    subplot(1, 2, 2);
    plot_stacked_relative_changes(relative_change_avoidance_all2, relative_change_approach_all2, 'EMS (p > 0.05)', relative_change_avoidance_all1, relative_change_avoidance_all2, relative_change_approach_all1, relative_change_approach_all2, false);
    axis square;
end

% Function to get the prefix list based on the date
function prefix_list = get_prefix_list(date_str, deb_dates, prez_dates)
    if ismember(date_str, deb_dates)
        prefix_list = {'deb', 'deb1'};
    elseif ismember(date_str, prez_dates)
        prefix_list = {'prez', 'prez1'};
    else
        error('Date %s not found in any prefix list.', date_str);
    end
end

% Function to get the existing prefix by checking file existence
function name = get_existing_prefix(prefix_list, date_str)
    for prefix = prefix_list
        if exist([prefix{1}, '_BarSize_', date_str, '.txt'], 'file')
            name = prefix{1};
            return;
        end
    end
    error('File for date %s not found with any prefix.', date_str);
end

% Function to load data based on prefix, date, and data type
function data = load_data(prefix, date_str, data_type)
    filename = [prefix, '_', data_type, '_', date_str, '.txt'];
    data = load(filename);
end

% Function to determine block sizes based on the length of data
function blockSizes = determine_block_sizes(lenA)
    switch lenA
        case 2700, blockSizes = [900, 900, 900];
        case 2400, blockSizes = [800, 800, 800];
        case 2100, blockSizes = [700, 700, 700];
        case 4000, blockSizes = [1500, 1500, 1000];
        case 3298, blockSizes = [1500, 1500, 298];
        case 3148, blockSizes = [1500, 1500, 148];
        case 2602, blockSizes = [1300, 1300, 2];
        case 2560, blockSizes = [800, 800, 960];
        case 2416, blockSizes = [800, 800, 816];
        otherwise
            blockSizes = [floor(lenA / 3), floor(lenA / 3), lenA - 2 * floor(lenA / 3)];
    end
end

% Function to process blocks and filter included trials
function Included_trials = process_block(A, B, block_indices)
    Red_bar = A(block_indices(1:2:end), :);
    Yellow_bar = A(block_indices(2:2:end), :);
    Trial_type = B(block_indices(1:2:end), :);
    Included_trials = [Red_bar, Yellow_bar, Trial_type];
    Included_trials(any(Included_trials(:, 1:2) == 0, 2) | ismember(Included_trials(:, 3), 1:5), :) = [];
end

% Function to calculate relative changes between two sets of included trials
function [relative_change_avoidance, relative_change_approach] = calculate_relative_changes(Included_trials_1, Included_trials_2, num_parts)
    [avoidance_percentages_1, approach_percentages_1] = calculate_percentages(Included_trials_1, num_parts);
    [avoidance_percentages_2, approach_percentages_2] = calculate_percentages(Included_trials_2, num_parts);

    valid_avoidance_indices = ~isnan(avoidance_percentages_1) & ~isnan(avoidance_percentages_2);
    valid_approach_indices = ~isnan(approach_percentages_1) & ~isnan(approach_percentages_2);

    relative_change_avoidance = calculate_relative_change(avoidance_percentages_1(valid_avoidance_indices), avoidance_percentages_2(valid_avoidance_indices));
    relative_change_approach = calculate_relative_change(approach_percentages_1(valid_approach_indices), approach_percentages_2(valid_approach_indices));
end

% Function to calculate percentages of avoidance and approach trials
function [avoidance_percentages, approach_percentages] = calculate_percentages(Included_trials, num_parts)
    part_size = floor(length(Included_trials) / num_parts);
    avoidance_percentages = NaN(1, num_parts);
    approach_percentages = NaN(1, num_parts);

    for i = 1:num_parts
        window_start = (i - 1) * part_size + 1;
        window_end = min(i * part_size, length(Included_trials));
        window_trials = Included_trials(window_start:window_end, :);
        if isempty(window_trials)
            avoidance_percentages(i) = NaN;
            approach_percentages(i) = NaN;
        else
            avoidance_percentages(i) = (sum(window_trials(:, 3) == 6) / length(window_trials)) * 100;
            approach_percentages(i) = (sum(window_trials(:, 3) == 0) / length(window_trials)) * 100;
        end
    end
end

% Function to calculate relative changes
function relative_change = calculate_relative_change(percentages_1, percentages_2)
    relative_change = ((percentages_2 - percentages_1) ./ percentages_1) * 100;
end

% Function to plot relative changes as stacked bars with significance markers
function plot_stacked_relative_changes(relative_change_avoidance_all, relative_change_approach_all, name, relative_change_avoidance_all1, relative_change_avoidance_all2, relative_change_approach_all1, relative_change_approach_all2, plot_significance)
    % Calculate average relative changes
    if ~isempty(relative_change_avoidance_all)
        avg_relative_change_avoidance = mean(relative_change_avoidance_all, 1, 'omitnan');
    else
        avg_relative_change_avoidance = NaN(1, 5);
    end

    if ~isempty(relative_change_approach_all)
        avg_relative_change_approach = mean(relative_change_approach_all, 1, 'omitnan');
    else
        avg_relative_change_approach = NaN(1, 5);
    end

    % Combine the average relative changes into a single matrix
    combined_avg_relative_changes = [avg_relative_change_avoidance; avg_relative_change_approach]';

    % Plot the average relative changes as stacked bar plots
    b = bar(combined_avg_relative_changes, 'stacked');
    
    % Adjust y-axis limits
    ylim([-50 100]);

    % Add labels and title
    title(name);
    xlabel('Segments');
    ylabel('Average Relative Change (%)');
    legend('Avoidance Trials', 'Approach Trials');
    grid on;
    set(gca, 'FontSize', 18);

    % Set the background color to white
    set(gca, 'Color', 'white');
    box off

    % Add significance markers if plot_significance is true
    if plot_significance
        for i = 1:size(relative_change_avoidance_all1, 2)
            [~, p_avoidance] = ttest2(relative_change_avoidance_all1(:, i), relative_change_avoidance_all2(:, i));
            [~, p_approach] = ttest2(relative_change_approach_all1(:, i), relative_change_approach_all2(:, i));
            
            % Add significance markers for avoidance trials (blue bars)
            if p_avoidance < 0.01
                if combined_avg_relative_changes(i, 1) > 0
                    text(i - 0.15, combined_avg_relative_changes(i, 1) + 10, '**', 'FontSize', 18, 'Color', 'black', 'HorizontalAlignment', 'center');
                else
                    text(i - 0.15, combined_avg_relative_changes(i, 1) - 10, '**', 'FontSize', 18, 'Color', 'black', 'HorizontalAlignment', 'center');
                end
            elseif p_avoidance < 0.05
                if combined_avg_relative_changes(i, 1) > 0
                    text(i - 0.15, combined_avg_relative_changes(i, 1) + 5, '*', 'FontSize', 18, 'Color', 'black', 'HorizontalAlignment', 'center');
                else
                    text(i - 0.15, combined_avg_relative_changes(i, 1) - 5, '*', 'FontSize', 18, 'Color', 'black', 'HorizontalAlignment', 'center');
                end
            end
            
            % Add significance markers for approach trials (red bars)
            if p_approach < 0.01
                if combined_avg_relative_changes(i, 2) > 0
                    text(i + 0.15, combined_avg_relative_changes(i, 2) + 10, '**', 'FontSize', 18, 'Color', 'black', 'HorizontalAlignment', 'center');
                else
                    text(i + 0.15, combined_avg_relative_changes(i, 2) - 10, '**', 'FontSize', 18, 'Color', 'black', 'HorizontalAlignment', 'center');
                end
            elseif p_approach < 0.05
                if combined_avg_relative_changes(i, 2) > 0
                    text(i + 0.15, combined_avg_relative_changes(i, 2) + 5, '*', 'FontSize', 18, 'Color', 'black', 'HorizontalAlignment', 'center');
                else
                    text(i + 0.15, combined_avg_relative_changes(i, 2) - 5, '*', 'FontSize', 18, 'Color', 'black', 'HorizontalAlignment', 'center');
                end
            end
        end
    end
end
