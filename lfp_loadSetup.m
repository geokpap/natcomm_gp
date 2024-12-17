% lfp_loadSetup

%About lfp_setup
%
% Alternative versions of lfp_setup can be invoked by name; see
% lfp_changeSetup. 
%
% At a minimum, the setup script lfp_setup.m must assign values to the
% following variables:
%   CSCFileRegexp - a regular expression that matches any CSC file name
%       (just the filename and the dot that separates the filename and
%       extension; no directory) (not required if useFileSelect is set)
%   CSCFileExt - CSC file extension, without the '.'; may not contain any
%       characters that are special to regular expressions
%   EVFilename1 - the name of the preferred event file, relative to
%       lfp_DataDir (not required if useFileSelect is set)
%   EVFilename2 - an alternate event file name (not required if
%       useFileSelect is set)
%   read_mode:
%       'nlx' for Neuralynx files
%       'mat' for Matlab files
%       'preset' to skip GUI (see lfp_read2)
%   lfp_NominalTrialStart - (see lfp_read2)
%   lfp_NominalTrialEnd -  (see lfp_read2)
%   lfp_SetupType - (see lfp_read2)
%
% 11-Nov-2007:  DG adopted convention that lfp_setup_<name> should
% explicitly assign values only to globals, and any other assignments
% should reside in lfp_getEvtIDs_<name>, which should be implicitly called
% by calling lfp_getEvtIDs early in the setup file; any variable whose name
% is not a member of lfp_GlobalNames should thus be moved from
% lfp_setup_<name> to lfp_getEvtIDs_<name>.  This makes it possible to
% obtain values from lfp_getEvtIDs without altering the values of any
% globals.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

if isempty(lfp_SetupName)
    if exist('lfp_SetupName.mat') == 2
        load('lfp_SetupName.mat');
    else
        error('lfp_getEvtIDs:noSetup', ...
            'You must use a non-default setup to use this function');
    end
end

% Read lfp_setup script.  Note that any values set in that script for
% persistent global variables will overwrite the saved values.
if isempty(lfp_SetupName)
    disp('Please enter an lfp_SetupName name surround by single-quote marks,');
    lfp_SetupName = input('or just hit Enter to use default: ');
end
if isempty(lfp_SetupName)
    mysetup = 'lfp_setup';
else
    mysetup = ['lfp_setup_' lfp_SetupName];
end
try
    eval([mysetup ';']);
catch
    error('lfp_loadSetup:badsetup', [ ...
            'Could not run the setup file %s.\nThis could mean that the' ...
            ' file does not exist or contains an error.' ], mysetup );
end

