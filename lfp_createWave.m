function lfp_createWave(fhandle, filenums, varargin)
%lfp_createWave(fhandle, filenums, varargin)
%lfp_createWave(..., 'name', name)
%lfp_createWave(..., 'replace', filenum)

%  Creates a new waveform (i.e. a new filenum) in the same format as CSC
%  files.  <fhandle> is a function handle that points to a function that
%  accepts a vector of filenums and possibly additional arguments, and
%  produces a vector of length (length(lfp_TimeStamps) *
%  lfp_SamplesPerFrame) that is used as the value of the new waveform.
%  <fhandle> may optionally return a second result which
%  is a string representing the units of the new waveform; this string
%  is assigned to the corresponding element of lfp_SamplesUnits, which
%  would otherwise be 'arbs'.
%  Options recognized here together with any values that may follow will be
%  excised from varargin before passing varargin on to <fhandle>.
%  <filenums> can represent anything, it is simply the first arg passed to
%  <fhandle>.
% Options:
%   'name' - must be followed by another string to use as the
%       lfp_FileName of the new channel.  
%   'replace' - must be followed by a filenum whose contents will be
%       DELETED AND REPLACED with the new waveform.
%   'units', string - ignores the second return value from <fhandle> (if
%       any), and sets the units of the output signal to be <string>.

%$Rev: 379 $
%$Date: 2016-03-30 15:21:20 -0400 (Wed, 30 Mar 2016) $
%$Author: dgibson $

lfp_declareGlobals;
nameflag = false;
optunits = '';
replaceflag = false;
unitsflag = false;
argnum = 1;
args2delete = [];
while argnum <= length(varargin)
    if strcmp(class(varargin{argnum}), 'char')
        switch varargin{argnum}
            case 'name'
                nameflag = true;
                args2delete = [args2delete argnum:argnum+1];
                argnum = argnum + 1;
                name = varargin{argnum};
            case 'replace'
                replaceflag = true;
                args2delete = [args2delete argnum:argnum+1];
                argnum = argnum + 1;
                newfilenum = varargin{argnum};
            case 'units'
                args2delete = [args2delete argnum:argnum+1];
                argnum = argnum + 1;
                optunits = varargin{argnum};
        end
    end
    argnum = argnum + 1;
end
if ~isempty(args2delete)
    varargin(args2delete) = [];
end

if ~nameflag
    name = [func2str(fhandle) '(' dg_thing2str(filenums) ')'];
end
if replaceflag
    % Save undo info in case of error
    old_lfp_SelectedFile = lfp_SelectedFiles(newfilenum);
    old_lfp_FileName = lfp_FileNames{newfilenum};
else
    newfilenum = max(lfp_ActiveFilenums) + 1;
    old_ActiveFilenums = lfp_ActiveFilenums;
    lfp_ActiveFilenums(end+1) = newfilenum;
end
lfp_SelectedFiles(newfilenum) = true;
lfp_FileNames{newfilenum} = name;
newunits = '';

try
    try
        [lfp_Samples{newfilenum} newunits] = ...
            feval(fhandle, filenums, varargin{:});
    catch
        s = lasterror;
        if ismember(s.identifier, ...
                {'MATLAB:maxlhs' 'MATLAB:TooManyOutputs'} )
            lfp_Samples{newfilenum} = feval(fhandle, filenums, varargin{:});
        else
            rethrow(s);
        end
    end
catch e
    % undo the changes and rethrow the error
    if replaceflag
        lfp_SelectedFiles(newfilenum) = old_lfp_SelectedFile;
        lfp_FileNames{newfilenum} = old_lfp_FileName;
    else
        lfp_ActiveFilenums = old_ActiveFilenums;
        if newfilenum == length(lfp_Samples)
            lfp_Samples(newfilenum) = [];
        end
    end
    logmsg = sprintf('Error: while processing %s', char(fhandle));
    logmsg = sprintf('%s\n%s\n%s', ...
        logmsg, e.identifier, e.message);
    for stackframe = 1:length(e.stack)
        logmsg = sprintf('%s\n%s\nline %d', ...
            logmsg, e.stack(stackframe).file, ...
            e.stack(stackframe).line);
    end
    disp(logmsg);
    rethrow(e);
end
if ~isempty(optunits)
    newunits = optunits;
end
if isempty(newunits)
    lfp_SamplesUnits{newfilenum} = 'arbs';
else
    lfp_SamplesUnits{newfilenum} = newunits;
end
if numel(lfp_Samples{newfilenum}) ~= numel(lfp_Samples{lfp_ActiveFilenums(1)})
    warning('lfp_createWave:badlength', ...
        'The function returned length %d, which differs from others %d', ...
        numel(lfp_Samples{newfilenum}), numel(lfp_Samples{lfp_ActiveFilenums(1)}) );
end
arginstring = dg_thing2str(filenums);
for argnum = 1:length(varargin)
    arginstring = [ arginstring ', ' dg_thing2str(varargin{argnum}) ];
end
lfp_log(sprintf('Created new filenum %d: %s(%s)', newfilenum, ...
    func2str(fhandle), arginstring ));
