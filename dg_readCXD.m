function [TS, frame, framenum] = dg_readCXD(filepath)
%INPUTS
% filepath: absolute or relative pathname to a Hamamatsu CXD format video
%   file.
%OUTPUTS
% TS: time since start of recording in seconds, column vector.
% frame: image data in cell column vector.
% framenum: row vector of integer frame numbers derived from <TS>.
%NOTES
% This function requires the Bio-Formats toolbox for Matlab, available at:
% http://downloads.openmicroscopy.org/bio-formats/5.5.3/
% Further info at:
% https://www-legacy.openmicroscopy.org/site/support/bio-formats5.1/developers/matlab-dev.html

%$Rev: 262 $
%$Date: 2017-10-16 16:32:55 -0400 (Mon, 16 Oct 2017) $
%$Author: dgibson $ 

data = bfopen(filepath);
metadata = data{1, 2};
metadataKeys = metadata.keySet().iterator();
key = cell(metadata.size(), 1);
for k=1:metadata.size();
    key{k} = metadataKeys.nextElement();
end
Time_From_Start_idx = find( ~cellfun(@isempty, ...
    strfind(key, 'Time_From_Start')) );
TFSvalue = zeros(length(Time_From_Start_idx), 1);
for k = 1:length(Time_From_Start_idx)
    TFSvalue(k) = metadata.get(key{Time_From_Start_idx(k)});
end
toks = regexp(key(Time_From_Start_idx)', ...
    'Global Field (\d+) Time_From_Start', 'tokens');
TFSnumeric = reshape(str2double(cellfun(@(x) x{1}, toks)), [], 1);
[~, sortTFSnumericidx] = sort(TFSnumeric);
TS = TFSvalue(sortTFSnumericidx);
dTS = diff(TS);
framedur = median(diff(TS)); 
% Check for missing frames:
toolongidx = find(dTS > 1.5 * framedur);
framenum = [];
nextFrameNum = 1;
msg = 'Dropped frame(s) between timestamps:';
if isempty(toolongidx)
    framenum = 1:length(TS);
else
    for idx2 = 1:length(toolongidx)
        numdropped = round(dTS(toolongidx(idx2))/framedur) - 1;
        if idx2 == 1
            numgood = toolongidx(idx2);
        else
            numgood = toolongidx(idx2) - toolongidx(idx2-1);
        end
        msg = sprintf( '%s\n%.3f \t%.3f', msg, ...
            TS(toolongidx(idx2)), TS(toolongidx(idx2) + 1) );
        framenum = [framenum nextFrameNum + (1:numgood) - 1]; %#ok<AGROW>
        nextFrameNum = framenum(end) + 1 + numdropped;
    end
    numgood = length(TS) - toolongidx(end);
    framenum = [framenum nextFrameNum + (1:numgood) - 1];
end
if ~isempty(toolongidx)
    warning('dg_readCXD:drop', '%s', msg);
end
if ~isequal(framenum, TFSnumeric(sortTFSnumericidx))
    warning('dg_readCXD:framenum', ...
        '<framenum> disagrees with recorded frame numbers.');
end

frame = data(:, 1);


