function rawvalues = dg_writeMIDI(dirname, myevents, taskrelated, nonrelated, allstages)
% Create output files, one unit per file, names like <stage><t|n><FileID>
% where <t|n> denotes t for task-related, n for non-task-related.

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

rawvalues = NaN(14*30*50,1);
rawvaluesoffset = 0;
outdir = fullfile(dirname, 'output');
mymax = -Inf;
mymin = Inf;
for stageix = 1:length(allstages)
    stage = allstages(stageix);
    multivalues = [];
    multifileIDs = cell(0);
    for tn = 'tn'
        if tn == 't'
            blankrow = find(all(isnan(taskrelated{1, stageix}), 2));
            if numel(blankrow) ~= 1
                warning('not exactly one blank t row');
                blankrow = blankrow(1);
            end
        else
            blankrow = find(all(isnan(nonrelated{1, stageix}), 2));
            if numel(blankrow) ~= 1
                warning('not exactly one blank n row');
                blankrow = blankrow(1);
            end
        end
        for unitrow = 2:blankrow-1
            if tn == 't'
                fileID = taskrelated{1, stageix}(unitrow, 1);
            else
                fileID = nonrelated{1, stageix}(unitrow, 1);
            end
            values = NaN(length(myevents)*40, 1);
            for evtix = 1:length(myevents)
                if tn == 't'
                    values((1:40) + (evtix-1) * 40) = ...
                        taskrelated{evtix, stageix}( ...
                        unitrow, 135:174);
                else
                    values((1:40) + (evtix-1) * 40) = ...
                        nonrelated{evtix, stageix}( ...
                        unitrow, 135:174);
                end
            end
            mymax = max(mymax, max(values));
            mymin = min(mymin, min(values));
            rawvalues(rawvaluesoffset + (1:numel(values))) = values;
            rawvaluesoffset = rawvaluesoffset + numel(values);
            values = dg_MIDIscalefunc(values);
            if any(isnan(values))
                warning('No file output due to missing values for unit %d event %d stage %d', ...
                    fileID, myevents(evtix), stage);
            else
                if any(values<1)
                    warning('clipping on minus side');
                    values(values<1) = 1;
                end
                if any(values>127)
                    warning('clipping on plus side');
                    values(values>127) = 127;
                end
                outfilename = sprintf('%02d%c%d.txt', stage, tn, fileID);
                fid = fopen(fullfile(outdir, outfilename), 'w');
                fprintf(fid, '%d\n', values);
                fprintf(fid, '0\n');
                fclose(fid);
            end
            multivalues = [ multivalues values ];
            multifileIDs = [ multifileIDs {[tn sprintf('%d', fileID)]} ];
        end
    end
    % Write multicolumn file and track names file
    colheaders = multifileIDs{1};
    if length(multifileIDs) > 1
        for colnum = 2:length(multifileIDs)
            colheaders = sprintf('%s\t%s', ...
                colheaders, multifileIDs{colnum} );
        end
    end
    outfilename = sprintf('%02dmultihdr.txt', stage);
    fid = fopen(fullfile(outdir, outfilename), 'w');
    fprintf(fid, '%s\n', colheaders);
    fclose(fid);
    outfilename = sprintf('%02dmulti.txt', stage);
    fid = fopen(fullfile(outdir, outfilename), 'w');
    for row = 1:size(multivalues,1)
        rowstr = sprintf('%d\t', multivalues(row,:));
        rowstr(end) = sprintf('\n');
        fwrite(fid, rowstr);
    end
    fclose(fid);
end
disp(sprintf('mymax=%d, mymin=%d', mymax, mymin));

