function trigTS = dg_convertStimSession(ctrlfilepath, fnames, varargin)
%trigTS = dg_convertStimSession(ctrlfilepath, fnames)
% Invokes dg_rmStimArtifact on a set of CSC files and also creates a
% stimulation-based events file suitable for use by lfp_lib.  Each "trial"
% starts at <trigTS> as returned by dg_rmStimArtifact, and ends halfway to
% the next <trigTS>, except for the last which is made equal to the median
% trial duration.  The events file is compatible with the lfp_lib setup
% named 'stim', whose setup and getEvtIDs files can be found in the
% 'lfp_setups' folder.  Its name is 'events.mat' and it will overwrite any
% pre-existing 'events.mat'.
%   The default value of <offsets> given to dg_rmStimArtifact is [0 .0015].
%INPUTS
% ctrlfilepath: path to a CSC file that contains data used to control the
%   timing of the interpolations.
% fnames: if empty, only a new .mat format events file is created; if a
%   string, then any file in <sessiondir> is processed if its name matches
%   the string according to the operating system's matching conventions; if
%   a cell array, it is passed verbatim to dg_rmStimArtifact.
%OPTIONS
% All options other than the ones listed below are passed through to
% dg_rmStimArtifact. 
%   'offsets', offsets - overrides the default value.
%   'dest', dest - the directory to which 'events.mat' gets written, which
%       by default is the same directory as <ctrlfilepath>.  This option is
%       also passed through to dg_rmStimArtifact.
%EXAMPLES
% dg_convertStimSession( ...
%   '/smbshare/monkey/v9/sakura/2011-08-08_23-11-21/csc124.ncs', ...
%   {'csc15.ncs'}, 'offsets', [.0005 .001]);
%       Creates a new 'events.mat' file in
%       /smbshare/monkey/v9/sakura/2011-08-08_23-11-21/csc124.ncs, with one
%       trial per stimulation pulse, and removes the corresponding artifact
%       from 'csc15.ncs' and saves the resulting waveform in
%       'csc15rmstim.mat'.  The artifact is presumed to extend from .0005 s
%       after the trigger signal to .001 s after, and that time period is
%       replaced by a linear interpolation between the end points.
%NOTES
% someday should handle .nev and .nse files as trig sources too

%$Rev:  $
%$Date:  $
%$Author: dgibson $

ctrldir = fileparts(ctrlfilepath);
offsets = [0 .0015];
dest = ctrldir;

opts2delete = [];
argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'dest'
                argnum = argnum + 1;
                dest = varargin{argnum};
            case 'offsets'
                opts2delete(end+1) = argnum; %#ok<*AGROW>
                argnum = argnum + 1;
                offsets = varargin{argnum};
                opts2delete(end+1) = argnum; %#ok<*AGROW>
        end
    else
        error('funcname:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end
varargin(opts2delete) = [];

if isempty(fnames)
    filenames = {};
else
    switch class(fnames)
        case 'char'
            files = dir(fullfile(ctrldir, fnames));
            filenames = {files.name};
        case 'cell'
            filenames = fnames;
        otherwise
            error('dg_convertStimSession:fnames', ...
                '<fnames> must be empty, a string, or a cell array.');
    end
end

trigTS = dg_rmStimArtifact(ctrlfilepath, filenames, ...
    offsets, varargin{:});
endTS = (trigTS(1:end-1) + trigTS(2:end)) / 2;
endTS(end+1) = trigTS(end) + median(endTS);

dg_Nlx2Mat_Timestamps(2:2:2*length(endTS)) = endTS;
dg_Nlx2Mat_Timestamps(1:2:end) = trigTS; %#ok<NASGU>
dg_Nlx2Mat_TTL(2:2:2*length(endTS)) = 2;
dg_Nlx2Mat_TTL(1:2:end) = 1; %#ok<NASGU>
dg_Nlx2Mat_EventStrings = repmat({'dg_convertStimSession'}, ...
    2*length(endTS), 1); %#ok<NASGU>
matfilepath = fullfile(dest, 'events.mat');
save(matfilepath, 'dg_Nlx2Mat_Timestamps', 'dg_Nlx2Mat_TTL', ...
    'dg_Nlx2Mat_EventStrings');


