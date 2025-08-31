function [allvalues, taskrelated, nonrelated] = dg_PH2MIDI(dirname)
% Yasuo PH*.xls files to MIDI track files

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

myevents = [10 11 13 39 54 56 78];  % order determines output order
% myevents = [10  78];  % order determines output order
acqstages = 1:15;
OTstages = 17:22;
reacqstages = 24:29;
% acqstages = 1:2;
% OTstages = []
% reacqstages = 29;

allstages = [acqstages OTstages reacqstages];
taskrelatedfmt = 'PH%dData6-1.xls';
nonrelatedfmt = 'PH%dData26-1.xls';

taskrelated = cell(length(myevents), length(allstages));
nonrelated = cell(length(myevents), length(allstages));

% Read spreadsheets
for evtix = 1:length(myevents)
    evt = myevents(evtix);
    for sheetix = 1:length(allstages)
        sheet = allstages(sheetix);
        filename = sprintf(taskrelatedfmt, evt);
        taskrelated{evtix, sheetix} = ...
            xlsread(fullfile(dirname, filename), sprintf('%d', sheet));
        filename = sprintf(nonrelatedfmt, evt);
        nonrelated{evtix, sheetix} = ...
            xlsread(fullfile(dirname, filename), sprintf('%d', sheet));
    end
end

if nargout == 0
    save('dg_PH2MIDI_values', 'taskrelated', 'nonrelated');
    return
end

allvalues = dg_writeMIDI(dirname, myevents, taskrelated, nonrelated, allstages);
allvalues(isnan(allvalues)) = [];

end


