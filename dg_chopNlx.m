function dg_chopNlx(src, dest, varargin)
%dg_chopNlx(src, dest)
% Finds disordered timestamps in <src> file, and copies the portion after
% the last disordered timestamp to <dest>.

%$Rev: 207 $
%$Date: 2014-10-16 19:07:56 -0400 (Thu, 16 Oct 2014) $
%$Author: dgibson $

maxrecs = 10000;
argnum = 1;
ncsflag = false;
nevflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'ncs'
            ncsflag = true;
            tmpext = '.ncs';
        case 'nev'
            nevflag = true;
            tmpext = '.nev';
        otherwise
            error('funcname:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
if ncsflag || nevflag
    tempname = dg_mktempname('dg_chopNlx', tmpext);
    copyfile(src, tempname);
    src = tempname;
end
[path, name, ext] = fileparts(src);
switch lower(ext)
    case '.ncs'
        [Timestamps, Header] = Nlx2MatCSC_411( src, [1 0 0 0 0], 1, 1, 1 );
    case '.nev'
        [Timestamps, Header] = Nlx2MatEV_411( src, [1 0 0 0 0], 1, 1, 1 );
    otherwise
        error('dg_chopNlx:filetype', ...
            'Cannot handle files of type %s', ext);
end
if ~isempty(regexp(Header{end}, '^\s*$'))
    Header(end) = [];
end
numrecs = length(Timestamps);
backwardsTSidx = find(diff(Timestamps) < 0);
disp(sprintf('There are %d reverse-ordered timestamps', ...
    length(backwardsTSidx) ));
if length(backwardsTSidx) > 0
    disp(sprintf('Records %s', mat2str(backwardsTSidx)));
    lastbackwards = backwardsTSidx(end);
    firstrec = lastbackwards + 1;
    lastrec = lastbackwards;
else
    firstrec = 0;
    lastrec = 0;
end
appendflag = 0;
hWaitBar = waitbar(0, 'Progress');
while lastrec < numrecs
    lastrec = min(lastrec + maxrecs, numrecs);
    disp(sprintf('Reading %d of %d records', lastrec, numrecs));
    switch lower(ext)
        case '.ncs'
            [Timestamps, ChanNum, SampleFrequency, ...
                NumberValidSamples, Samples] = Nlx2MatCSC_411( ...
                src, [1 1 1 1 1], 0, 2, [firstrec lastrec] );
            if ~ispc
                error('dg_chopNlx:notPC', ...
                    'Neuralynx format files can only be created on Windows machines.');
            end
            Mat2NlxCSC_411( dest, appendflag, 1, 1, length(Timestamps), ...
                [1 1 1 1 1 1], Timestamps, ChanNum, SampleFrequency, ...
                NumberValidSamples, Samples, Header );
        case '.nev'
            [Timestamps, EventIDs, Nttls, Extras, EventStrings] = ...
                Nlx2MatEV_411( ...
                src, [1 1 1 1 1], 0, 2, [firstrec lastrec] );
            if ~ispc
                error('dg_chopNlx:notPC', ...
                    'Neuralynx format files can only be created on Windows machines.');
            end
            Mat2NlxEV_411( dest, appendflag, 1, 1, length(Timestamps), ...
                [1 1 1 1 1 1], Timestamps, EventIDs, Nttls, Extras, ...
                EventStrings, Header );
    end
    if ~appendflag
        appendflag = 1;
    end
    waitbar( lastrec/numrecs, hWaitBar, ...
        sprintf('Read %d of %d records', lastrec, numrecs), hWaitBar );
    firstrec = lastrec + 1;
end
close(hWaitBar);
if ncsflag || nevflag
    delete(tempname);
end
