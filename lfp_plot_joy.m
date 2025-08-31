function lfp_plot_joy(trialevents, trial, startsample, endsample)
lfp_declareGlobals;
%draw joy pos trial, assume fnum 3 = JoyX
%need to read gain and calibration (i.e., center coords) vals from
%BDdata (must read lo and hi parts and combine them to get calib
%vals)

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

% Joystick gain is hard-coded and equals 3 for both x and y;

lfp_getJoyCalibs;

% extract target position from BDOutput
lfp_getTaskParams;

figure('color',[0 0 0]);

axis equal
%         axis([-400 400 -400 400]);
hold on;
% plot the target array --must figure out the color and
% shape of each target on this trial.

% CHANGE THIS TO BE A SEPARATE ROUTINE THAT WE CAN CALL

R = 60; %target radius
x = [-R:1:R];
y = (R^2 - x.*x).^0.5;        


R2 = 30; %inner target radius to draw white circle
x2 = [-R2:1:R2];
y2 = (R2^2 - x2.*x2).^0.5;

ShapeSpacing = 200;

for i = 0:8
    switch (i)
        case {0,1,3}
            y_offset = 0;
        case {5,6}
            y_offset = (2^0.5)*ShapeSpacing/2;
        case {7,8}
            y_offset = -(2^0.5)*ShapeSpacing/2; %lower shapes should be negative
        case {2}
            y_offset = ShapeSpacing;
        case {4}
            y_offset = -1*ShapeSpacing;
    end
    switch (i)
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
    
    switch (i)
        case {1,2,3,4,5,6,7,8}
            % plot the target circles
            % need two plots to do top and bottom half of
            % target circle in red
            plot(x + x_offset, y + y_offset, 'Color', [1, 1, 1]);
            plot(x + x_offset, -y + y_offset, 'Color', [1, 1, 1]);
            
        case {0}
            % plot the center target in yellow    
            plot(x + x_offset, y + y_offset, 'Color', [1, 1, 0]);
            plot(x + x_offset, -y + y_offset, 'Color', [1, 1, 0]);
            if CenterHold1 == 1
                plot(x2 + x_offset, y2 + y_offset, 'Color', [1, 1, 0]);
                plot(x2 + x_offset, -y2 + y_offset, 'Color', [1, 1, 0]);
            end                    
    end
    if i == Target1
        plot(x + x_offset, y + y_offset, 'Color', [1, 0, 0]);
        plot(x + x_offset, -y + y_offset, 'Color', [1, 0, 0]);
        if CenterHold2 == 1
            plot(x2 + x_offset, y2 + y_offset, 'Color', [1, 0, 0]);
            plot(x2 + x_offset, -y2 + y_offset, 'Color', [1, 0, 0]);
        end
    end
    if i == Target2
        plot(x + x_offset, y + y_offset, 'Color', [0, 1, 0]);
        plot(x + x_offset, -y + y_offset, 'Color', [0, 1, 0]);
        if CenterHold3 == 1
            plot(x2 + x_offset, y2 + y_offset, 'Color', [0, 1, 0]);
            plot(x2 + x_offset, -y2 + y_offset, 'Color', [0, 1, 0]);
        end                
    end 
    if i == Target3
        plot(x + x_offset, y + y_offset, 'Color', [0, 0, 1]);
        plot(x + x_offset, -y + y_offset, 'Color', [0, 0, 1]);
        if CenterHold4 == 1
            plot(x2 + x_offset, y2 + y_offset, 'Color', [0, 0, 1]);
            plot(x2 + x_offset, -y2 + y_offset, 'Color', [0, 0, 1]);
        end                
    end
    % number the targets
    %             text(x_offset, y_offset, num2str(i), 'HorizontalAlignment', 'center');
end

% extract a range of sample indeces between two events of
% interest, in order to plot the data in color (coding for time)
coloroneventtime = trialevents(find(trialevents(:,2) == 16),1);
coloroffeventtime = trialevents(find(trialevents(:,2) == 34),1);

coloroneventindex = lfp_time2index(coloroneventtime);
coloroffeventindex = lfp_time2index(coloroffeventtime);

colortrialsamplerange = coloroneventindex : coloroffeventindex;
greytrialsamplerange = startsample : (coloroneventindex - 1);
blacktrialsamplerange = (coloroffeventindex + 1) : endsample;

oC = 0;       
aD = 4096 / 10; % Delphi (behavioral control computer) A/D board gain
aC = (4096 / 4)*(1/5); % Cheetah Analog-In (BNC) gain (can also be read from (ADBITVOLT*AmpGain) from the CSC file-header)
% with 5:1 attenuator
oD = 2048;

% Reconstruct the output of the A/D board in the behavioral control
% computer, based on the recorded Cheetah data
Delphi_samples4c = ((lfp_Samples{4}(colortrialsamplerange)) - oC) * aD/aC + oD;
Delphi_samples3c = ((lfp_Samples{3}(colortrialsamplerange)) - oC) * aD/aC + oD;

Delphi_samples4b = ((lfp_Samples{4}(blacktrialsamplerange)) - oC) * aD/aC + oD;
Delphi_samples3b = ((lfp_Samples{3}(blacktrialsamplerange)) - oC) * aD/aC + oD;

% collect the data preceding the color to plot in grey, and the
% data following the color to be in black
% Now perform same calibration as is done in Delphi, but add the
% offset from the upper left corner of the monkeysee (or mirror)
% form, so that the data appear at the center of the Matlab figures
% (-1930,-1755);

calibrated_samples4c = -3*(Delphi_samples4c - XCalib) %- 1930;
calibrated_samples3c = -(Delphi_samples3c - YCalib) * 3% - 1755;

% flipdim here will cause the plot to portray elapsed time from RED to BLUE.
%         proportional_trialtime = flipdim([0:1:length(calibrated_samples4c) - 1]', 1);

proportional_trialtime = [1:1:length(calibrated_samples4c)];

% flip the colormap
%colormap(flipdim(colormap,1));

scatter(calibrated_samples4c', calibrated_samples3c', 4, proportional_trialtime, 'filled');
colorbar;

% JoyY (I think JoyX and Y were flipped on 1/11/04)

% btw, should consider flipping the order of these data s.t. when
% plotted in color, using scatter, they will run from red thru
% green to blue, following the actual progession of the colors of
% the target-capture sequence.

title([lfp_DataDir 'Trial' lfp_getTrialID(trial) 'Joy']);
