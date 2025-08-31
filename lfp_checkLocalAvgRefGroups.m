function [filenames, isbadfile, grpfiles, grpresults, reportstr] = ...
    lfp_checkLocalAvgRefGroups(localavgrefs, sessiondir, alignment, varargin)
% Performs various tests on the CSC files in each local average reference
% group specified by <localavgrefs>, and determines for each file whether
% its contents look sufficiently plausibly LFP-like to include in the local
% average reference.  Specifically, lfp_markBadCSC is run (or the contents
% of 'lfp_markBadCSC.mat' is read if it already exists), and then each
% group of "good" channels is tested as a group to determine whether the
% files look sufficiently similar to each other to be plausibly derived
% from a single local neighborhood.  A human-readable report of the results
% is produced on standard output and is also returned in <reportstr>.
%INPUTS
% localavgrefs: as for dg_makeLocalAvgRefs.
% sessiondir: as for dg_makeLocalAvgRefs.
% alignment: a value to use for lfp_AlignmentRef when running the
%   similarity tests.
%OUTPUTS
% filenames: cell column vector of strings as returned by lfp_markBadCSc.
% isbadfile: boolean column vector as returned by lfp_markBadCSc.
% grpfiles: cell array of the same shape as <localavgrefs>, where each cell
%   contains a column vector of the literal filenames belonging to the
%   group.
% grpresults: cell array of the same shape as <localavgrefs>, where each
%   cell contains a logical column vector the same size as the
%   corresponding cell in <grpfiles>, with the value set to <true> if the 
%   corresponding file was marked bad by lfp_markBadCSC OR by the
%   within-group tests.
% reportstr: the human-readable report.
%OPTIONS
% 'destdir', destdir - see lfp_markBadCSC.
% 'matchExpr', matchExpr - see lfp_markBadCSC.
% 'selection', selectionrule - specifies a <selectionrule> for submission
%   to lfp_selectByRule prior to applying lfp_selectByDuration (see NOTES).
%   Default is 'HasEvent(lfp_AlignmentRef)'.
% 'verbose' - somewhat.
%NOTES
%   There are two similarity tests:
%   1. Compute the zero-lag result of lfp_xcov on whole trials between
% every pair of files within each local avg ref group; then for each
% electrode, average those values across all pairs of which the electrode
% was a member.  The electrode is marked "bad" if its average xcov is below
% 0.25.  Note that this result depends on the value of <alignment> and the
% 'selection' option (if given).  lfp_selectByDuration is used to eliminate
% unusually short trials with beginnings or endings in the [-10 10] range.
% Each channel is highpass filtered at 1 Hz before calling lfp_xcov.
%   2. For each channel, compute the 95th percentile of the absolute
% difference between sample values and their mean value.  The same set of
% samples is used as is used for the lfp_xcov test. Any channel whose 95th
% prctile differs from the group median by more than a factor of two is
% marked "bad".
%   Note that <grpfiles> and <grpresults> may contain a different total
% number of file references from <filenames>/<isbadfile> because the same
% file could potentially be repeated in different groups, and the files
% that are included in all groups put together may be only a subset of
% <filenames>.

%$Rev: 389 $
%$Date: 2017-06-12 15:07:03 -0400 (Mon, 12 Jun 2017) $
%$Author: dgibson $

global lfp_AlignmentRef lfp_XLimAll lfp_SelectedTrials lfp_EventNames ...
    lfp_Xxcov lfp_SamplePeriod lfp_ActiveFilenums

markBadCSC_opts = {};
maxfactor = 3;
minxcov = 0.25;
selectionrule = 'HasEvent(lfp_AlignmentRef)';
verboseflag = false;

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'destdir'
            markBadCSC_opts{end+1} = 'destdir';
            argnum = argnum + 1;
            markBadCSC_opts{end+1} = varargin{argnum};
        case 'matchExpr'
            markBadCSC_opts{end+1} = 'matchExpr';
            argnum = argnum + 1;
            markBadCSC_opts{end+1} = varargin{argnum};
        case 'selection'
            argnum = argnum + 1;
            selectionrule = varargin{argnum};
        case 'verbose'
            verboseflag = true;
        otherwise
            error('lfp_checkLocalAvgRefGroups:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

got_markBadCSC = exist(fullfile(sessiondir, 'lfp_markBadCSC.mat'), 'file');

if verboseflag
    if got_markBadCSC
        fprintf('Loading lfp_markBadCSC.mat\n');
    else
        fprintf('Running lfp_markBadCSC\n');
    end
end

reportstr = sprintf('Session directory: %s\n', sessiondir);
if got_markBadCSC
    badCSC = load(fullfile(sessiondir, 'lfp_markBadCSC.mat'));
    filenames = badCSC.filenames;
    isbadfile = badCSC.isbadfile;
    clear('badCSC');
    reportstr = [reportstr, sprintf('Loaded lfp_markBadCSC.mat.\n')];
else
    [filenames, isbadfile] = lfp_markBadCSC(sessiondir, markBadCSC_opts{:});
    reportstr = [reportstr, sprintf('Ran lfp_markBadCSC.\n')];
end
for grpidx = 1:length(localavgrefs)
    if isempty(localavgrefs{grpidx})
        warning('lfp_checkLocalAvgRefGroups:skipping', ...
            'Skipping empty group #%d', grpidx);
        continue
    end
    missingfn = ~ismember(localavgrefs{grpidx}, filenames);
    if any(missingfn)
        error('lfp_checkLocalAvgRefGroups:badfn', ...
            'File(s) %s is not in list of filenames from lfp_markBadCSC.', ...
            dg_thing2str(localavgrefs{grpidx}(missingfn)) );
    end
end
for badidx = reshape(find(isbadfile), 1, [])
    reportstr = [reportstr, sprintf( ...
        'Bad File: %s flagged by lfp_markBadCSC.\n', ...
        filenames{badidx} )]; %#ok<*AGROW>
end
badfns = filenames(isbadfile);

% Find the events file name <evtfname>:
evtfname = dg_findEventsFile(sessiondir);
reportstr = [reportstr, sprintf('Event File: %s\n', evtfname)];
grpfiles = cell(size(localavgrefs));
grpresults = cell(size(localavgrefs));
% Iterate over local average refs:
oldalign = lfp_AlignmentRef;
oldselect = lfp_SelectedTrials;
oldxlim = lfp_XLimAll;
win = [];

for refidx = 1:length(localavgrefs)
    grpresults{refidx} = false(0,0);
    CSCfilespec = localavgrefs{refidx};
    if isempty(CSCfilespec)
        continue
    end
    % Do the within-group tests on one local avg ref group.
    % <grpfnames> is the list of all files in the group.
    if ischar(CSCfilespec)
        filelist = dir(fullfile(sessiondir, CSCfilespec));
        grpfnames = reshape({filelist.name}, [], 1);
    elseif iscell(CSCfilespec)
        grpfnames = reshape(CSCfilespec, [], 1);
    else
        error('lfp_checkLocalAvgRefGroups:CSCfilespec', ...
            'Each element of <localavgrefs> must be a cell array or a string.');
    end
    grpfiles{refidx} = grpfnames;
    grpresults{refidx} = false(size(grpfnames));
    for grpfnidx = 1:length(grpfnames)
        grpresults{refidx}(grpfnidx) = isbadfile( ...
            ismember(filenames, grpfnames{grpfnidx}) );
    end
    % Remove bad files from <grpfnames> before analyzing.  This is where
    % <grpfnames> ceases to be the same as grpfiles{refidx}.
    grpfnames = setdiff(grpfnames, badfns);
    if length(grpfnames) < 2
        reportstr = [ reportstr, sprintf( ...
            'No group tests possible, group#%d=%s\n', ...
            refidx, dg_thing2str(grpfnames) ) ];
        continue
    end
    % iterate over pairs of files from <grpfnames>:
    matrix = NaN(length(grpfnames));
    prct = NaN(length(grpfnames), 1);
    for ix1 = 1:(length(grpfnames) - 1)
        for ix2 = ix1+1 : length(grpfnames)
            lfp_read2('preset', sessiondir, [{evtfname} ...
                grpfnames(ix1) grpfnames(ix2)]);
            lfp_AlignmentRef = alignment;
            lfp_selectByRule(selectionrule);
            lfp_XLimAll = [0 10]; %#ok<*NASGU>
            lfp_selectByDuration('and', 'noexpand');
            lfp_XLimAll = [-10 0];
            lfp_selectByDuration('and', 'noexpand');
            lfp_Xxcov = [-1 1] * lfp_SamplePeriod;
            lfp_XLimAll = [];
            lfp_createWave(@lfp_bandpass2, lfp_ActiveFilenums(1), 1, ...
                'high', 'replace', lfp_ActiveFilenums(1));
            lfp_createWave(@lfp_bandpass2, lfp_ActiveFilenums(2), 1, ...
                'high', 'replace', lfp_ActiveFilenums(2));
            [~, result] = lfp_xcov([], [], [], 'norm', 'noplot');
            rowzero = find(result(:,1) == 0);
            matrix(ix1, ix2) = result(rowzero, 2);
            matrix(ix2, ix1) = result(rowzero, 2);
            if verboseflag
                fprintf('Finished %s vs %s\n', ...
                    grpfnames{ix1}, grpfnames{ix2});
            end
        end
        [sampledata, timepts] = lfp_getSamples( [], ...
            lfp_ActiveFilenums(1), [] );
        prct(ix1) = prctile(abs(reshape( ...
            sampledata(:, :, 1) - nanmean(reshape( ...
            sampledata(:, :, 1), [], 1 )), ...
            [], 1 )), 95);
        if isempty(win)
            win = [timepts(1) timepts(end)];
            reportstr = [reportstr, sprintf( ...
                'lfp_AlignmentRef: %d (%s)\n', ...
                lfp_AlignmentRef, ...
                lfp_EventNames{lfp_AlignmentRef} )];
            reportstr = [reportstr, sprintf( ...
                'Trial analysis window: [%d %d]\n', ...
                timepts(1), timepts(end) )];
            reportstr = [reportstr, sprintf( ...
                'Number of selected trials: %d\n', ...
                length(lfp_enabledTrials) )];
        end
    end
    % We have still not computed the prct magnitude for the very last file:
    [sampledata, timepts] = lfp_getSamples( [], ...
        lfp_ActiveFilenums(2), [] );
    prct(ix1+1) = prctile(abs(reshape( ...
        sampledata(:, :, 1) - nanmean(reshape( ...
        sampledata(:, :, 1), [], 1 )), ...
        [], 1 )), 95);
    absmean = nanmean(abs(matrix));
    grpisbadfile = reshape(absmean < minxcov, [], 1) ...
        | prct < nanmedian(prct) / maxfactor | prct > maxfactor * nanmedian(prct);
    for grpbadidx = reshape(find(grpisbadfile), 1, [])
        if absmean(grpbadidx) < minxcov
            reportstr = [reportstr, sprintf( ...
                'Bad File: %s mean xcov = %.1d\n', ...
                grpfnames{grpbadidx}, absmean(grpbadidx) )]; %#ok<*AGROW>
        end
        if prct(grpbadidx) < median(prct) / maxfactor ...
                || prct(grpbadidx) > maxfactor * median(prct)
            reportstr = [reportstr, sprintf( ...
                'Bad File: %s 95th prct mag = %.1d (median = %.1d)\n', ...
                grpfnames{grpbadidx}, prct(grpbadidx), median(prct) )]; %#ok<*AGROW>
        end
        % Since <grpfnames> can be different from grpfiles{refidx} at this
        % point, we need to look up the correct entry in grpfiles{refidx}
        % by name: 
        grpfnidx = ismember(grpfiles{refidx}, grpfnames{grpbadidx});
        grpresults{refidx}(grpfnidx) = true;
    end
end
badfilevector = cell2mat(reshape(grpresults, [], 1));
fracbadfiles = sum(badfilevector)/length(badfilevector);
if fracbadfiles > 0.1
    reportstr = [ reportstr, sprintf( ...
        'Large fraction of files are "bad": %.3f\n', fracbadfiles ) ];
else
    reportstr = [ reportstr, sprintf( ...
        'Small fraction of files are "bad": %.3f\n', fracbadfiles ) ];
end
fprintf('%s', reportstr);

% Restore original global values.
lfp_AlignmentRef = oldalign;
lfp_SelectedTrials = oldselect;
lfp_XLimAll = oldxlim;


