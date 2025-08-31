function lfp_readBursts(filenames, fragname, varargin)
%lfp_readBursts(filenames)
% Create a logical CSC channel containing bursts read from each file listed
% in <filenames>.
%INPUTS
% filenames:  cell array of strings containing absolute or relative
%   pathnames to the files to be read.  Each file must contain the
%   following variables:
%       burstbounds:  a cell array with one element for each fragment in
%           the file.  Each element is a two column array of timestamps,
%           with burst start timestamps in column 1 and burst end
%           timestamps in column 2.
%       fragment:  a cell string array containing the names of the
%           fragments in the file.
% fragname:  a string or cell string array that matches one or more of the
%   elements of <fragments> as read from <filename>, which controls which
%   fragment(s) are read from the file.  
%OPTIONS
% 'dur' - do not create logical CSC channels; instead, create CSC channels
%   containing the burst duration, propagated to fill the non-bursting time
%   from the end of the burst to the onset of the next one.  The time
%   before the first burst in each rec segment is NaN, and if the segment
%   ends within a burst, then so is the time from the last burst onset to
%   the end.
% 'rate' - do not create logical CSC channels; instead, create CSC channels
%   containing the "instantaneous" burst rate, defined as the reciprocal of
%   the time between burst onsets between each pair of burst onsets, and
%   NaN before the first and after the last burst onset in each rec
%   segment.

%$Rev: 287 $
%$Date: 2012-10-26 20:14:45 -0400 (Fri, 26 Oct 2012) $
%$Author: dgibson $

global lfp_ActiveFilenums lfp_Samples lfp_RecSegments

durflag = false;
rateflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'dur'
            durflag = true;
        case 'rate'
            rateflag = true;
        otherwise
            error('lfp_readBursts:badoption', ...
                ['The option "' ...
                dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

for fileidx = 1:length(filenames)
    load(filenames{fileidx}, '-mat');
    if ~exist('burstbounds', 'var') || ~exist('fragment', 'var')
        error('lfp_readBursts:file', ...
            '%s is not a valid burst bounds file', filenames{fileidx});
    end
    fragidx = find(ismember(fragment, fragname)); %#ok<*USENS>
    if isempty(fragidx)
        fragstr = sprintf('\n%s', fragment{1});
        for k = 2:length(fragment)
            fragstr = sprintf('%s\n%s', fragstr, fragment{k});
        end
        error('lfp_readBursts:fragname', ...
            'For file %s, <fragname> must match one of:%s', ...
            filenames{fileidx}, fragstr);
    end
    [p, n] = fileparts(filenames{fileidx}); %#ok<ASGLU>
    for fragnum = reshape(fragidx, 1, [])
        if ~durflag && ~rateflag
            name = sprintf('%s:%s', n, fragment{fragnum});
            lfp_createWave(@lfp_mkBurstMask, burstbounds{fragnum}, ...
                'name', name);
            lfp_Samples{lfp_ActiveFilenums(end)} = ...
                logical(lfp_Samples{lfp_ActiveFilenums(end)});
        else
            if rateflag
                name = sprintf('%s:%s_rate', n, fragment{fragnum});
                data = NaN(size(lfp_Samples{lfp_ActiveFilenums(1)}));
                IBIs = diff(burstbounds{fragnum}, 1);
                for segidx = 1:size(lfp_RecSegments,1)
                    % Fill in the data within this segment
                    segstartTS = lfp_index2time(lfp_RecSegments(segidx,1));
                    segendTS = lfp_index2time(lfp_RecSegments(segidx,2));
                    insegidx = find(burstbounds{fragnum}(:,1) >= segstartTS & ...
                        burstbounds{fragnum}(:,1) <= segendTS);
                    for burstidx = insegidx(1) : insegidx(end-1)
                        data(lfp_time2index( ...
                            burstbounds{fragnum}(burstidx,1) ) ...
                            : lfp_time2index( ...
                            burstbounds{fragnum}(burstidx+1,1) ) ) = ...
                            1/IBIs(burstidx);
                    end
                end
                lfp_createWave(@lfp_waverecord, lfp_ActiveFilenums(1), ...
                    data, 'name', name, 'units', 'Hz');
            end
            if durflag
                name = sprintf('%s:%s_dur', n, fragment{fragnum});
                data = NaN(size(lfp_Samples{lfp_ActiveFilenums(1)}));
                burstdurs = diff(burstbounds{fragnum}, 1, 2);
                for segidx = 1:size(lfp_RecSegments,1)
                    % Fill in the data within this segment, handling the
                    % propagation of the last value at the end of the
                    % segment appropriately
                    segstartTS = lfp_index2time(lfp_RecSegments(segidx,1));
                    segendTS = lfp_index2time(lfp_RecSegments(segidx,2));
                    insegidx = find(burstbounds{fragnum}(:,1) >= segstartTS & ...
                        burstbounds{fragnum}(:,1) <= segendTS & ...
                        burstbounds{fragnum}(:,2) >= segstartTS & ...
                        burstbounds{fragnum}(:,2) <= segendTS);
                    for burstidx = insegidx(1) : insegidx(end-1)
                        data(lfp_time2index( ...
                            burstbounds{fragnum}(burstidx,1) ) ...
                            : lfp_time2index( ...
                            burstbounds{fragnum}(burstidx+1,1) ) ) = ...
                            burstdurs(burstidx);
                    end
                    if insegidx(end) >= size(burstbounds{fragnum}, 1) ...
                            || burstbounds{fragnum}(insegidx(end)+1, 1) ...
                            > segendTS
                        % The next burst is in a different rec segment, so
                        % we propagate the value to the end of this segment
                        data(lfp_time2index( ...
                            burstbounds{fragnum}(insegidx(end),1) ) ...
                            : lfp_RecSegments(segidx,2) ) = ...
                            burstdurs(insegidx(end));
                    else
                        % The next burst is in this segment, so we use it
                        % to terminate the value propagation
                        data(lfp_time2index( ...
                            burstbounds{fragnum}(insegidx(end),1) ) ...
                            : burstbounds{fragnum}(insegidx(end)+1, 1) ) ...
                            = burstdurs(insegidx(end));
                    end
                end
                lfp_createWave(@lfp_waverecord, lfp_ActiveFilenums(1), ...
                    data, 'name', name, 'units', 's');
            end
        end
    end
end
