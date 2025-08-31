function [sesreport, sesreportary] = dg_sessionClustTest(sessiondir, ...
    retrig, exceptions)
% Look recursively for .mat cluster files in <sessiondir> and submit each
% one to 'dg_runClustTests'.  All .mat files are assumed to be cluster
% files except for filenames listed in <exceptions>.
%INPUTS
% sessiondir: root of tree to search for cluster files to score.
% retrig: value for 'retrig' option to 'dg_shortISI', i.e. the "dead time"
%   after a spike trigger, during which there cannot be another trigger.
% exceptions: cell string array containing literal filenames of .mat files
%   that should be skipped.
%OUTPUT
% sesreport: same format as returned by dg_runClustTests, but with a new
%   column inserted before the first column returned by dg_runClustTests
%   containing the basename of the .mat file analyzed.  Note that this
%   means that all the column numbers specified in 'dg_runClustTests' have
%   to be incremented by one when addressing <sesreport>, <sesreportary>.
% sesreportary: same content as <sesreport> but formatted as a cell array
%   with one row per cluster.

%$Rev: 307 $
%$Date: 2023-10-06 13:06:15 -0400 (Fri, 06 Oct 2023) $
%$Author: dgibson $

sesreport = '';
sesreportary = {};
files = dir(sessiondir);
s = warning('query', 'dg_shortISI:nonstat');
warning('off', 'dg_shortISI:nonstat');
for fidx = 1:length(files)
    if files(fidx).isdir
        if ismember(files(fidx).name, {'.', '..'})
            continue
        end
        [newreportstr, newreportary] = dg_sessionClustTest( ...
            fullfile(sessiondir, files(fidx).name), retrig, exceptions );
    elseif ismember(files(fidx).name, exceptions)
        continue
    else
        [~, basename, ext] = fileparts(files(fidx).name);
        if isequal(lower(ext), '.mat') && isequal( ...
                upper(basename(1:3)), 'CSC' )
            [reportstr, reportary] = dg_runClustTests( ...
                fullfile(sessiondir, files(fidx).name), retrig );
            reportlines = strsplit(reportstr, '\n');
            % The last one is always empty because <reportstr> always ends
            % with '\n':
            reportlines(end) = [];
            % Insert new columns:
            newreportstr = '';
            for linenum = 1:length(reportlines)
                newreportstr = sprintf( '%s%s\t%s\n', ...
                    newreportstr, basename, ...
                    reportlines{linenum} );
            end
            newreportary = [ repmat({basename}, size(reportary, 1), 1) ...
                reportary ];
        else
            continue
        end
    end
    sesreport = sprintf('%s%s', sesreport, newreportstr);
    sesreportary = [ sesreportary
        newreportary ]; %#ok<AGROW>
end
warning(s.state, 'dg_shortISI:nonstat');


