% lfp_getEvtIDs
% Script to import named event ID constants into current workspace by
% running lfp_getEvtIDs_<setupname>; also sets value of lfp_SetupName.
%
% 11-Nov/2007 DG adopted convention that lfp_getEvtIDs_<setupname> is
% responsible for setting ALL user-specific symbolic constants.  Aside from
% event IDs, symbolic constants include variables that are specifically
% used by individual programs, e.g. lfp_EyeTabulation uses BlinkStart,
% BlinkEnd, etc., but these are not used anywhere else (as of 11-Nov-2007).
% These constants are (for the most part) noted in the comments of
% the programs that use them.  Symbols that are widely used in many
% programs and not necessarily commented as such include:
%   CSCFileRegexp: see lfp_loadSetup
%   CSCFileExt: see lfp_loadSetup
%   EVFilename1: see lfp_loadSetup
%   EVFilename2: see lfp_loadSetup
%   UnitlessFileRegexp: matched against absolute pathname of CSC file to
%       disable conversion from A/D units
%   useFileSelect: if true, the lfp_fileSelect GUI is used to select files
%       to load; otherwise, files are loaded automatically based on
%       CSCFileRegexp, etc. (see lfp_loadSetup).
%   usestrobe: if not specified, its default value is <true>.  Controls
%       whether or not the highest order bit ina 16-bit TTL event code is
%       interpreted as a strobe bit or a data bit when processing the event
%       data.

%$Rev: 351 $
%$Date: 2015-05-28 12:32:51 -0400 (Thu, 28 May 2015) $
%$Author: dgibson $

global lfp_SetupName

if isempty(lfp_SetupName)
    if exist('lfp_SetupName.mat', 'file')
        load('lfp_SetupName.mat');
    else
        error('lfp_getEvtIDs:noSetup', ...
            'You must use a non-default setup to use this function');
    end
end
myEvtIDs = ['lfp_getEvtIDs_' lfp_SetupName];
if ~exist(myEvtIDs, 'file')
    warning('lfp_getEvtIDs:noEvtIDs', ...
        'Could not find EvtIDs file %s.', myEvtIDs );
else
    try
        eval([myEvtIDs ';']);
    catch e
        logmsg = sprintf('%s\n%s', ...
            e.identifier, e.message);
        for stackframe = 1:length(e.stack)
            logmsg = sprintf('%s\n%s\nline %d', ...
                logmsg, e.stack(stackframe).file, ...
                e.stack(stackframe).line);
        end
        error('lfp_getEvtIDs:badEvtIDs', [ ...
            'Could not run the EvtIDs file %s.\n' ...
            '%s' ], ...
            myEvtIDs, logmsg );
    end
end
