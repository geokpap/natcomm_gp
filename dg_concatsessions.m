function dg_concatsessions(dir1, dir2, dir3, varargin)
% Concatenates files from session directories <dir1> and <dir2> into <dir3>
% (which gets created if it does not already exist).  Only works on .ncs,
% .nev, .ntt, and .nse files (ignores all others).  The sessions must each
% contain at least one CSC file with the same name. Before writing any
% files to <dir3>, checks first timestamp of the lexically first CSC file
% in <dir2> to verify that it is at least one frame later than the last
% timestamp in the matching CSC file in <dir1>, and raises an error if not.
% Raises a warning if there are files that exist in only one of <dir1> and
% <dir2>, and proceeds to ignore them.
%OPTIONS
% 'ext', extfilter - <extfilter> is a string that must match the file
%   extension exactly, e.g. '.nse'.
%NOTES
% This function currently only works under Windows.

%$Rev: 271 $
%$Date: 2020-07-24 15:52:50 -0400 (Fri, 24 Jul 2020) $
%$Author: dgibson $

argnum = 0;
extfilter = '';
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'ext'
            argnum = argnum + 1;
            extfilter = varargin{argnum};
        otherwise
            error('dg_concatsessions:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

files1 = dir(dir1);
files2 = dir(dir2);
for k = 1:length(files1)
    [~, ~, ext1{k, 1}] = fileparts(files1(k).name); %#ok<AGROW>
end
for k = 1:length(files2)
    [~, ~, ext2{k, 1}] = fileparts(files2(k).name); %#ok<AGROW>
end

% Find files to concatenate:
exts = {'.ncs' '.nev' '.ntt' '.nse'};  % '.ncs' must be first
for e = 1:length(exts)
    ext = exts{e};
    gotext1 = ismember(ext1, ext);
    fn1 = cell(sum(gotext1), 1);
    [fn1{:}] = deal(files1(gotext1).name);
    gotext2 = ismember(ext2, ext);
    fn2 = cell(sum(gotext2), 1);
    [fn2{:}] = deal(files2(gotext2).name);
    commonnames{e} = intersect(fn1, fn2); %#ok<AGROW>
    if length(fn1) > length(commonnames{e})
        warning('dg_concatsessions:extra1', ...
            '%s contains %s files that do not exist in %s', ...
            dir1, ext, dir2);
    end
    if length(fn2) > length(commonnames{e})
        warning('dg_concatsessions:extra2', ...
            '%s contains %s files that do not exist in %s', ...
            dir2, ext, dir1);
    end
end

if isempty(commonnames{1})
    error('dg_concatsessions:noCSC', ...
        'There are no CSC files to concatenate.');
end

% Do the timestamp test:
ts1 = Nlx2MatCSC_411(fullfile(dir1, commonnames{1}{1}), [1 0 0 0 0], 0, 1);
ts2 = Nlx2MatCSC_411(fullfile(dir2, commonnames{1}{1}), [1 0 0 0 0], 0, 1);
if ts2(1) < ts1(end) + median(diff(ts1))
    error('dg_concatsessions:ts', ...
            '%s does not start after end of %s', dir2, dir1);
end
clear ts1 ts2

if ~exist(dir3, 'dir')
    mkdir(dir3);
end

% Do the concatenations:
batchsize = 2^16;
for e = 1:length(exts)
    for f = 1:length(commonnames{e})
        values = cell(4,1);
        ext = exts{e};
        if ~isempty(extfilter) && ~isequal(ext, extfilter)
            continue
        end
        fn = commonnames{e}{f};
        disp(['Concatenating ' fn '...']);
        copyfile(fullfile(dir1, fn), fullfile(dir3, fn));
        startrec = 0;
        endrec = -1;
        switch ext
            case '.ncs'
                ts = Nlx2MatCSC_411(fullfile(dir2, commonnames{e}{f}), ...
                    [1 0 0 0 0], 0, 1);
            case '.nev'
                % bizarrely, Mat2NlxEV_411 does not write event strings in
                % append mode unless you also "write" the header (which of
                % course *never* writes in append mode).
                [ts, hdr] = Nlx2MatEV_411(fullfile(dir2, commonnames{e}{f}), ...
                    [1 0 0 0 1], 0, 1);
            case {'.ntt' '.nse'}
                ts = Nlx2MatSpike_411(fullfile(dir2, commonnames{e}{f}), ...
                    [1 0 0 0 0], 0, 1);
        end
        numrecs2 = length(ts);
        while endrec < numrecs2 - 1
            endrec = min(startrec + batchsize - 1, numrecs2 - 1);
            switch ext
                case '.ncs'
                    [values{:}] = ...
                        Nlx2MatCSC_411(fullfile(dir2, commonnames{e}{f}), ...
                        [0 1 1 1 1], 0, 2, [startrec endrec]);
                    if ~ispc
                        error('dg_concatsessions:notPC', ...
                            'Neuralynx format files can only be created on Windows machines.');
                    end
                    Mat2NlxCSC_411(fullfile(dir3, fn), 1, 1, 1, ...
                        endrec - startrec + 1, [1 1 1 1 1 0], ...
                        ts((startrec:endrec)+1), ...
                        values{:});
                case '.nev'
                    [values{:}] = ...
                        Nlx2MatEV_411(fullfile(dir2, commonnames{e}{f}), ...
                        [0 1 1 1 1], 0, 2, [startrec endrec]);
                % bizarrely, Mat2NlxEV_411 does not write event strings in
                % append mode unless you also "write" the header (which of
                % course *never* writes in append mode).
                    if ~ispc
                        error('dg_concatsessions:notPC', ...
                            'Neuralynx format files can only be created on Windows machines.');
                    end
                    Mat2NlxEV_411(fullfile(dir3, fn), 1, 1, 1, ...
                        endrec - startrec + 1, [1 1 1 1 1 1], ...
                        ts((startrec:endrec)+1), ...
                        values{:}, hdr);
                case {'.ntt' '.nse'}
                    [values{:}] = ...
                        Nlx2MatSpike_411(fullfile(dir2, commonnames{e}{f}), ...
                        [0 1 1 1 1], 0, 2, [startrec endrec]);
                    if ~ispc
                        error('dg_concatsessions:notPC', ...
                            'Neuralynx format files can only be created on Windows machines.');
                    end
                    if isequal(ext, '.ntt')
                        Mat2NlxTT_411(fullfile(dir3, fn), 1, 1, 1, ...
                            endrec - startrec + 1, [1 1 1 1 1 0], ...
                            ts((startrec:endrec)+1), ...
                            values{:});
                    else % '.nse'
                        Mat2NlxSE_411(fullfile(dir3, fn), 1, 1, 1, ...
                            endrec - startrec + 1, [1 1 1 1 1 0], ...
                            ts((startrec:endrec)+1), ...
                            values{:});
                    end
            end
            startrec = endrec + 1;
        end
    end
end

disp('Done');


            
