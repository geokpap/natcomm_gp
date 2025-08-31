function lfp_perf_t2(trials, perfflag, perfblocksize)
% calculates the moving average versus chunks of average

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

%%%%%%%%%%%NEED TO FIX THIS WHOLE THING SO IT NO LONGER USES 'CTASK' BUT
%%%%%%%%%%%USES CTASK NAME INSTEAD

lfp_declareGlobals;

trial = 1; %solely to keep lfp_getBDParams_theresa happy.
lfp_getBDParams_theresa;

% use this just to plot vertical lines to plot trials when there were task
% parameter changes
plottaskchanges = true;
% taskchanges = [160 209 306 307 329 510 518 1393 1668]; %G092404
% taskchanges = [140 147 283]; %G092604
% taskchanges = [1171 1211 1473]; %G092904
% taskchanges = [287 525 1628]; %G100404
% taskchanges = [256 459 1213]; %G100704
% taskchanges = [606 833 989 1239 1455 1581]; %G100804
% taskchanges = [549 576 606 685]; %G101004
% taskchanges = [51 103 107 114 286 377 540 727 1115 1390]; %G101104
taskchanges = [];

% Scan task changes
ScanTargsIndex_old = ScanTargsIndex; %for trial 1
ScanTargsChanges = [];
ScanDelayTime_edit_old = ScanDelayTime_edit;
ScanDelayTime_editChanges = [];

% initialize all the variables needed

% number of columns throughout these is 10 because ctask = 0-9 currently
%vector to hold the counts of how many trials have been counted for each
%task type:
trialindex = zeros(1,10); 
%matrix to hold the performance on individual trials per task
% trialperf = zeros((2 * perfblocksize), 10);
trialperf = [];
%matrix to hold the performance in blocks - not sure how big it's going to
%be, but perhaps can figure out from total number of trials?
taskperf = [];
%matrix to hold the total number of correct and trials in a chunk of trials
%for an accurate percentage average for that chunk
correct_over_total = zeros(2,10);

% trying this for new section to keep track of calibration changes
ploteyechanges = true;

XCalibEye_old = XCalibEye; %starting value for calibration - should be trial 1
XCalibEyeChanges = [];
YCalibEye_old = YCalibEye;
YCalibEyeChanges = [];

XGainEyeL_old = XGainEye1L/XGainEye2L;
XGainEyeLChanges = [];
XGainEyeR_old = XGainEye1R/XGainEye2R;
XGainEyeRChanges = [];
YGainEyeU_old = YGainEye1U/YGainEye2U;
YGainEyeUChanges = [];
YGainEyeD_old = YGainEye1D/YGainEye2D;
YGainEyeDChanges = [];

% these are the coloros that each respective ctask will be displayed in.
% (Separated them to make them easier to read)
taskcolor(0 +1,1:3) = [0 0 1];
taskcolor(1 +1,1:3) = [0 .5 0];
taskcolor(2 +1,1:3) = [1 0 0];
taskcolor(3 +1,1:3) = [0 .75 .75];
taskcolor(4 +1,1:3) = [.75 0 .75];
taskcolor(5 +1,1:3) = [.75 .75 0];
taskcolor(6 +1,1:3) = [.5 .5 .5];
taskcolor(7 +1,1:3) = [.8 .4 0];
taskcolor(8 +1,1:3) = [0 .8 .2];
taskcolor(9 +1,1:3) = [0 0 0];
% initialize the ctask_old to the first trial
ctask_old = ctask;

% calculate the variables needed for the moving average
a = 1;
% b = [1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10];
for i = 1:perfblocksize
    b(i) = 1/perfblocksize;
end

perfblock = 1;
perfblockref = 1;

totalindex = 0;
totalindexref = 1;

tempavg = 0;
tempsum = 0;

figure;
hold on;

for trial = trials
    if lfp_SelectedTrials(trial) == 1 % do this loop only for selected trials
        lfp_getBDParams_theresa;        
        
        trialindex(ctask+1) = trialindex(ctask+1) + 1;
        switch trialperformance
            case 0
                thistrialperf = 0;
            case 1
                thistrialperf = 1;
            case 2
                thistrialperf = 0;
        end
        trialperf(trialindex(ctask+1),(ctask+1)) = thistrialperf;
        
        % calculate the percentage and put it in the array along with
        % shuffling things around if it's full
%         if (ctask == ctask_old) && (trialindex(ctask+1) == (2 * perfblocksize))
%             
%             %sum the total correct trials
%             tempsum = sum(trialperf(1:perfblocksize,(ctask+1)));
%             
%             %put the totals in the matrix for future averageing use
%             correct_over_total(1,(ctask+1)) = correct_over_total(1,(ctask+1)) + tempsum;
%             correct_over_total(2,(ctask+1)) = correct_over_total(2,(ctask+1)) + perfblocksize;
%             
%             % put the performance for that group in the task performance
%             % matrix in the column for the appropriate ctask
%             taskperf(perfblock,(ctask+1)) = tempsum / perfblocksize;
%             perfblock = perfblock + 1;
%             
%             % move the matrix around so that the ones just calculated are
%             % removed
%             trialperf(1:perfblocksize, (ctask+1)) =...
%                 trialperf((perfblocksize+1):(2 * perfblocksize), (ctask+1));
%             trialperf((perfblocksize+1):(2 * perfblocksize), (ctask+1)) = zeros(perfblocksize, 1);
%             
%             % reset the index
%             trialindex(ctask+1) = perfblocksize;
%             
%         elseif ctask ~= ctask_old


        if ctask ~= ctask_old
             
            % add up the ones that haven't been tabulated yet
            tempsum = sum(trialperf(1:trialindex(ctask_old+1),(ctask_old+1)));
            
            %put the totals in the matrix for future averageing use
            correct_over_total(1,(ctask_old+1)) =...
                correct_over_total(1,(ctask_old+1)) + tempsum;
            correct_over_total(2,(ctask_old+1)) =...
                correct_over_total(2,(ctask_old+1)) + trialindex(ctask_old+1);            
            
%             % put the performance for the remainder of this chunk in task
%             % performance
%             taskperf(perfblock,(ctask_old+1)) = tempsum / trialindex(ctask_old+1);
%             
%             % plot this current chunk of data
%             plot(perfblockref:perfblock, taskperf(perfblockref:perfblock,(ctask_old+1)), ...
%                 'Color', taskcolor((ctask_old+1), 1:3), 'Marker', '*');
            
            % try the moving average stuff           
            x = trialperf(1:trialindex(ctask_old+1), (ctask_old+1));
            y = filter(b, a, x);
            
            % only want to plot the points that have a complete average
            % behind them
            totalindex = totalindex + length(x) - perfblocksize;
            plot(totalindexref:totalindex, y(perfblocksize+1:length(y)), ...
                'Color', taskcolor((ctask_old+1), 1:3), 'Marker', '*');


            % calculate the average for the chunk just plotted, print
            % text above the plot, and zero the matrix
            tempavg = (correct_over_total(1,(ctask_old+1)) / correct_over_total(2,(ctask_old+1))) * 100;
%             text(perfblockref, taskperf(perfblockref,(ctask_old+1)) + .025, ...
%                 [num2str(tempavg, 3) '%']);
            % puts the label above the first max value for this chunk
            text(totalindexref + min(find(y == max(y)) - perfblocksize), max(y) + .025, [num2str(tempavg, 3) '%']);
            correct_over_total(:,(ctask_old+1)) = 0;
            
%             perfblock = perfblock + 1;
%             % do this after increment it so it starts on the next one
%             perfblockref = perfblock;
            
            %reset the ref so it works the next time through
            totalindexref = totalindex + 1;
            
            % no longer want to hold in "memory" because they have been
            % tabulated so just zero the variables
            trialperf(1:trialindex(ctask_old+1), (ctask_old+1)) = zeros(trialindex(ctask_old+1), 1);
            trialindex(ctask_old+1) = 0;          

        end
        
        ctask_old = ctask;

    end
    
    % add a section on here to keep track of calibration changes
    if XCalibEye ~= XCalibEye_old
        XCalibEyeChanges = [XCalibEyeChanges; trial, XCalibEye];
    end
    XCalibEye_old = XCalibEye;
    if YCalibEye ~= YCalibEye_old
        YCalibEyeChanges = [YCalibEyeChanges; trial, YCalibEye];
    end
    YCalibEye_old = YCalibEye;
    
    if (XGainEye1L/XGainEye2L) ~= XGainEyeL_old
        XGainEyeLChanges = [XGainEyeLChanges; trial, (XGainEye1L/XGainEye2L)];
    end
    XGainEyeL_old = (XGainEye1L/XGainEye2L);
    if (XGainEye1R/XGainEye2R) ~= XGainEyeR_old
        XGainEyeRChanges = [XGainEyeRChanges; trial, (XGainEye1R/XGainEye2R)];
    end
    XGainEyeR_old = (XGainEye1R/XGainEye2R);
    if (YGainEye1U/YGainEye2U) ~= YGainEyeU_old
        YGainEyeUChanges = [YGainEyeUChanges; trial, (YGainEye1U/YGainEye2U)];
    end
    YGainEyeU_old = (YGainEye1U/YGainEye2U);
    if (YGainEye1D/YGainEye2D) ~= YGainEyeD_old
        YGainEyeDChanges = [YGainEyeDChanges; trial, (YGainEye1D/YGainEye2D)];
    end
    YGainEyeD_old = (YGainEye1D/YGainEye2D);
    
    % Section to keep track of changes in Scan Task
    switch CurrentCtaskName
        case {'Sc'}
            if ScanTargsIndex_old ~= ScanTargsIndex
                ScanTargsChanges = [ScanTargsChanges; trial, ScanTargsIndex];
            end
            ScanTargsIndex_old = ScanTargsIndex;
            if ScanDelayTime_edit_old ~= ScanDelayTime_edit
                ScanDelayTime_editChanges = [ScanDelayTime_editChanges; trial, ScanDelayTime_edit];
            end
            ScanDelayTime_edit_old = ScanDelayTime_edit;
    end
end

% print out when the calibration changed
XCalibEyeChanges
YCalibEyeChanges
XGainEyeLChanges
XGainEyeRChanges
YGainEyeUChanges
YGainEyeDChanges
ScanTargsChanges
ScanDelayTime_editChanges

% need something here that will calculate one last taskperf for the last
% block of ctask (trialindex > 11)
for i = 1:10
    if (trialindex(i) ~= 0) && (trialindex(i) ~= perfblocksize)
        
        % add up the ones that haven't been tabulated yet
        tempsum = sum(trialperf(1:trialindex(i),(i)));
        
        %put the totals in the matrix for future averageing use
        correct_over_total(1,i) =...
            correct_over_total(1,i) + tempsum;
        correct_over_total(2,i) =...
            correct_over_total(2,i) + trialindex(i);
        
        %         % put the performance for the remainder of this chunk in task
        %         % performance        
        %         taskperf(perfblock,i) = tempsum / trialindex(i);
        %         
        %         % plot this current chunk of data
        %         plot(perfblockref:perfblock, taskperf(perfblockref:perfblock, i), ...
        %             'Color', taskcolor(i, 1:3), 'Marker', '*');
        
        %moving average stuff
        x = trialperf(1:trialindex(i), i);
        y = filter(b, a, x);
        
        % only want to plot the points that have a complete average
        % behind them
        totalindex = totalindex + length(x) - perfblocksize;
        plot(totalindexref:totalindex, y(perfblocksize+1:length(y)), ...
            'Color', taskcolor(i, 1:3), 'Marker', '*');
        
        % calculate the average for the chunk just plotted, print
        % text above the plot, and zero the matrix
        tempavg = (correct_over_total(1,i) / correct_over_total(2,i)) * 100;
        %         text(perfblockref, taskperf(perfblockref,i) + .025, ...
        %             [num2str(tempavg, 3) '%']);
        % puts the label above the first max value for this chunk
        text(totalindexref + min(find(y == max(y)) - perfblocksize), max(y) + .025, [num2str(tempavg, 3) '%']);
        correct_over_total(:,i) = 0;
        
        %         perfblock = perfblock + 1;
        
        %reset the ref so it works the next time through
        totalindexref = totalindex + 1;
        
        % no longer want to hold in "memory" because they have been
        % tabulated so just zero the variables
        trialperf(1:trialindex(i), (i)) = zeros(trialindex(i), 1);
        trialindex(i) = 0;
    end
end

% figure;
% hold on;

% plot(taskperf,'LineStyle', '-', 'Marker', '*')

% -1 position places the legend outside of the axes
% legend('ctask0','ctask1','ctask2','ctask3','ctask4','ctask5','ctask6','ctask7','ctask8','ctask9', -1);

for i = 1:10
    line([totalindex-40, totalindex-15],[1-(.025*i), 1-(.025*i)], ...
        'LineWidth', 2, 'Color', taskcolor(i, 1:3));
    text(totalindex-10, 1-(.025*i), int2str(i-1));
end

% want a way to draw lines where significant task parameters changed -
% maybe come up with a way to eventually pass this as a variable?
if plottaskchanges == true
    for i = 1:length(taskchanges)
        line([taskchanges(i)-perfblocksize, taskchanges(i)-perfblocksize],[0, 1], 'LineStyle', '--', 'color', 'b');
    end
    for i = 1:size(ScanDelayTime_editChanges,1) %number of rows in XCalibEyeChanges
        line([ScanDelayTime_editChanges(i,1)-perfblocksize, ScanDelayTime_editChanges(i,1)-perfblocksize],[0, 1], 'LineStyle', '--', 'color', 'b');
    end    
    for i = 1:size(ScanTargsChanges,1) %number of rows in XCalibEyeChanges
        line([ScanTargsChanges(i,1)-perfblocksize, ScanTargsChanges(i,1)-perfblocksize],[0, 1], 'LineStyle', '--', 'color', [0 .5 0]);
    end    
end

% draw lines also where the eye calibration changed - differentiating
% between offset changes and gain changes
if ploteyechanges == true
    for i = 1:size(XCalibEyeChanges,1) %number of rows in XCalibEyeChanges
        line([XCalibEyeChanges(i,1)-perfblocksize, XCalibEyeChanges(i,1)-perfblocksize],[0, 1], 'LineStyle', '--', 'color', 'y');
    end
    for i = 1:size(YCalibEyeChanges,1) %number of rows in YCalibEyeChanges
        line([YCalibEyeChanges(i,1)-perfblocksize, YCalibEyeChanges(i,1)-perfblocksize],[0, 1], 'LineStyle', '--', 'color', 'y');
    end
    for i = 1:size(XGainEyeLChanges,1)
        line([XGainEyeLChanges(i,1)-perfblocksize, XGainEyeLChanges(i,1)-perfblocksize],[0, 1], 'LineStyle', '--', 'color', 'r');
    end
    for i = 1:size(XGainEyeRChanges,1)
        line([XGainEyeRChanges(i,1)-perfblocksize, XGainEyeRChanges(i,1)-perfblocksize],[0, 1], 'LineStyle', '--', 'color', 'r');
    end
    for i = 1:size(YGainEyeUChanges,1)
        line([YGainEyeUChanges(i,1)-perfblocksize, YGainEyeUChanges(i,1)-perfblocksize],[0, 1], 'LineStyle', '--', 'color', 'r');
    end
    for i = 1:size(YGainEyeDChanges,1)
        line([YGainEyeDChanges(i,1)-perfblocksize, YGainEyeDChanges(i,1)-perfblocksize],[0, 1], 'LineStyle', '--', 'color', 'r');
    end
end

hA = gca;
hT = get(hA, 'Title');
figtitle = [lfp_DataDir ' Performance (groups of ' int2str(perfblocksize) ' trials)'];
set(hT, 'Interpreter', 'none', 'String', figtitle);

axis([0 totalindex 0 1]);

xlabel(['Trial Number (minus the first ' int2str(perfblocksize) ' trials of each trace)']);
ylabel('% Correct (moving average)');

% plot(fraction_samples_inside);
% trialavg_fraction_samples_inside = mean(fraction_samples_inside)
% line([0 length(find(lfp_SelectedTrials))], [trialavg_fraction_samples_inside trialavg_fraction_samples_inside], 'color', 'r', 'linewidth', 3);
% trialstd_fraction_samples_inside = std(fraction_samples_inside)
% line([0 length(find(lfp_SelectedTrials))], [trialavg_fraction_samples_inside...
%         - trialstd_fraction_samples_inside trialavg_fraction_samples_inside...
%         - trialstd_fraction_samples_inside], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
% line([0 length(find(lfp_SelectedTrials))], [trialavg_fraction_samples_inside...
%         + trialstd_fraction_samples_inside trialavg_fraction_samples_inside...
%         + trialstd_fraction_samples_inside], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
