function histogram2D = lfp_hist_j3(trials, filenums, avgflag, histflag, figflag, calflag)

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

    alltrialsamples = [];
    penis = 1;
    for trial = trials
        if lfp_SelectedTrials(trial) == 1 % do this loop only for selected trials
            eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
            
            % find 'reftime', the lfp_AlignmentRef absolute timestamp:
            trialevents = lfp_Events(eventrange,:);
            reftime = trialevents( ...
                find(ismember(trialevents(:,2), lfp_AlignmentRef)), ...
                1 );
            if length(reftime) == 0
                reftime = 0;
            else
                reftime = reftime(1);
            end
            refpoint = lfp_time2index(reftime);
            
            % The time range of data we want to process is the union of the nominal
            % trial interval with the lfp_XLimAll interval, limited to the bounds
            % of the trial's recorded time segment.
            if isempty(lfp_XLimAll)
                xlimpoints = [0 0];
            else
                xlimpoints = round(lfp_XLimAll/lfp_SamplePeriod);
            end
            if reftime ~= 0
                xlimabspoints = xlimpoints + refpoint;
            else
                % There was no reference event in this trial. The following becomes
                % a harmless redundant value when we compute startsample and
                % endsample:
                xlimabspoints = [lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4)];
            end 
            startsample = max(...
                min( lfp_TrialIndex(trial,3), xlimabspoints(1)), ...
                lfp_TrialRec(trial,1) );
            endsample = min(...
                max( lfp_TrialIndex(trial,4), xlimabspoints(2)), ...
                lfp_TrialRec(trial,2) );
            trialsamplerange = startsample : endsample;
            
            if histflag 
                % Collect trial info to be used below for computing averaged
                % waveform. First find the sample index that is closest in time to
                % the lfp_AlignmentRef event:
                if (reftime - lfp_index2time(refpoint-1)) ...
                        < (lfp_index2time(refpoint) - reftime)
                    refpoint = refpoint - 1;
                end
                % Add to the trialinfo table, which contains the sample indices of
                % the start of trial, end of trial, and reference event:
                trialinfo = [ trialinfo
                    [ lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4) refpoint ] ];
                
                % need to change these indeces for BDFormatNumber = 8
                trialstarttime = trialevents(find(trialevents(:,2) == 2),1);
                
                % extract target position from BDOutput
                lfp_getTaskParams;
                
                % extract a range of sample indeces between two events of
                % interest, in order to plot the data in color (coding for time)
                % add option of user-input for these events
                coloroneventtime = trialevents(find(trialevents(:,2) == 16),1);
                coloroffeventtime = trialevents(find(trialevents(:,2) == 17),1);
                
                coloroneventindex = lfp_time2index(coloroneventtime);
                coloroffeventindex = lfp_time2index(coloroffeventtime);
                
                colortrialsamplerange = coloroneventindex : coloroffeventindex;
                greytrialsamplerange = startsample : (coloroneventindex - 1);
                blacktrialsamplerange = (coloroffeventindex + 1) : endsample;
                
                c = lfp_getEyeCalib(trial);
                calibrated_samples1c = ...
                    lfp_calibEye(lfp_Samples{1}(colortrialsamplerange), 'x', c);
                calibrated_samples2c = ...
                    lfp_calibEye(lfp_Samples{2}(colortrialsamplerange), 'y', c);
                
                % Minus sign because imagesc shows Y increasing down by
                % default:
                calibrated_samplesc = [calibrated_samples1c; -calibrated_samples2c]';
                
                % find maximum time range present in all trials
                pointsbefore = min(trialinfo(:,3) - trialinfo(:,1));
                pointsafter = min(trialinfo(:,2) - trialinfo(:,3));
                if ~isempty(lfp_XLimAll)
                    % convert lfp_XLimAll into points, truncate to maximum range
                    % present
                    pointsbefore = min(pointsbefore, ...
                        -round(lfp_XLimAll(1)/lfp_SamplePeriod));
                    pointsafter = min(pointsafter, ...
                        round(lfp_XLimAll(2)/lfp_SamplePeriod));
                end
                % compute averaged waveform for each filenum
                % avgwave = zeros(pointsbefore + pointsafter + 1, length(filenums));
                colnum = 0;
                
                allchannelsamples = [];
                
            end
            
            %distance from center of center target to center of peripheral
            %target
            ShapeSpacing = 200; % should read this value from BDoutput
            
            % extract target position from BDOutput
            lfp_getTaskParams2;
            
            Distance = [];
            for targetindex = 1:4
                
                switch (Target(targetindex))
                    case {0,1,3}
                        y_offset = 0;
                    case {5,6}
                        y_offset = -(2^0.5)*ShapeSpacing/2;
                    case {7,8}
                        y_offset = (2^0.5)*ShapeSpacing/2; %lower shapes should be negative
                    case {2}
                        y_offset = -ShapeSpacing;
                    case {4}
                        y_offset = 1*ShapeSpacing;
                end
                switch (Target(targetindex))
                    case {0,2,4}
                        x_offset = 0;
                    case {6,7}
                        x_offset = (2^0.5)*ShapeSpacing/2;
                    case {5,8}
                        x_offset = -(2^0.5)*ShapeSpacing/2;
                    case {3}
                        x_offset = ShapeSpacing;
                    case {1}
                        x_offset = -1*ShapeSpacing;
                end
                
                
                %calculate distance of samples from each of the four cues
                %(center, red, green, blue)
                
                Distance(targetindex,:) = (((calibrated_samplesc(:,1))' - x_offset).^2 + ((calibrated_samplesc(:,2))' - y_offset).^2).^0.5;
            end
             
            InTargets = Distance < 60;
            
            %fraction of samples in the interval of interest in this trial
            %that are inside cue(s)
            fraction_samples_inside(penis) = sum(sum(InTargets)) / size(calibrated_samplesc,1); 
            
            alltrialsamples = [alltrialsamples; calibrated_samplesc];
            penis = penis + 1;
        end    
    end
    figure;
    hold on;
    plot(fraction_samples_inside);
    trialavg_fraction_samples_inside = mean(fraction_samples_inside)
    line([0 length(find(lfp_SelectedTrials))], [trialavg_fraction_samples_inside trialavg_fraction_samples_inside], 'color', 'r', 'linewidth', 3);
    trialstd_fraction_samples_inside = std(fraction_samples_inside)
    line([0 length(find(lfp_SelectedTrials))], [trialavg_fraction_samples_inside...
            - trialstd_fraction_samples_inside trialavg_fraction_samples_inside...
            - trialstd_fraction_samples_inside], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
    line([0 length(find(lfp_SelectedTrials))], [trialavg_fraction_samples_inside...
            + trialstd_fraction_samples_inside trialavg_fraction_samples_inside...
            + trialstd_fraction_samples_inside], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
    title([lfp_DataDir 'Eye Fraction in Center(2-10)']);
    xlabel('Correct Trials');
    ylabel('% Time eyes within center cue');
    
    %calculate and plot 2D histogram of eye position over the interval of
    %interest
    figure;
    hist = hist2d([alltrialsamples(:,2) alltrialsamples(:,1)], linspace(-385,385,385+1),linspace(-513,513,513+1));
    pcolor(hist) % should add option of showing only "significant" bins, i.e., above some threshold
    % need to flip y axis because imagesc flips it for some reason
    imagesc(hist,[min(min(hist)) max(max(hist))])
    cmap = colormap;
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    title([lfp_DataDir 'EyeHist(2-10)']);
    %COLORBAR('horiz')
    
    hold on;
    
    targetarraycolor = [.5 .5 .5];
    % plot the target array
    lw = 2;
    R = 30; %target radius
    x = [-R+256:1:R+256];
    y_pos = (R^2 - (x-256).^2).^0.5 + 192;        
    y_neg = -(R^2 - (x-256).^2).^0.5 + 192;
    
    plot(x , y_pos, 'Color', [1, 1, 0], 'linewidth', lw);
    plot(x , y_neg, 'Color', [1, 1, 0], 'linewidth', lw);

            
    