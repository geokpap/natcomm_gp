function dg_Nlx2Mat(pathname, varargin)
%dg_Nlx2Mat(pathname)
% A wrapper for converting monkey data from Neuralynx format to .MAT format
% so they can be read on a Mac.

% <pathname> is either an absolute pathname to the file that needs
% converting, or it is relative to the current Matlab working directory.
% This function reads the file into as many Matlab variables as are
% required to hold all the data that are used by lfp_lib, and then saves
% them to a .MAT file of the same name in the same directory as the source
% file.  The variable names are as specified in the Neuralynx .M files but
% with 'dg_Nlx2Mat_' prepended.  Attempts to convert CSC files from AD
% units to Volts using the value of ADBitVolts from the CSC file header.
% Raises a warning and leaves in AD units if there is no ADBitVolts value
% in the header.
%
% OPTIONS
% 'dest', destdir - writes output to <destdir> instead of to the same
%   directory as the source file.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

destdir = '';
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'dest'
            argnum = argnum + 1;
            destdir = varargin{argnum};
        otherwise
            error('dg_Nlx2Mat:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

[pathstr,name,ext] = fileparts(pathname);
if isempty(destdir)
    destdir = pathstr;
end
ext = upper(ext);
if isequal(ext, '.DAT')
    if ~isempty(regexpi(name, 'events'))
        filetype = '.NEV';
    elseif ~isempty(regexpi(name, 'CSC|LFP'))
        filetype = '.NCS';
    end
else
    filetype = ext;
end
matfilename = fullfile(destdir, [ name '.mat' ]);
switch filetype
    case '.NCS'
        [dg_Nlx2Mat_Timestamps dg_Nlx2Mat_Samples header] ...
            = dg_readCSC(pathname);
        dg_Nlx2Mat_SamplesUnits = 'AD';
        for k = 1:length(header)
            if regexp(header{k}, '^\s*-ADBitVolts\s+')
                ADBitVoltstr = regexprep(header{k}, ...
                    '^\s*-ADBitVolts\s+', '');
                ADBitVolts = str2double(ADBitVoltstr);
                if isempty(ADBitVolts)
                    warning('dg_Nlx2Mat:badADBitVolts', ...
                        'Could not convert number from:\n%s', ...
                        header{k} );
                else
                    dg_Nlx2Mat_Samples = ADBitVolts ...
                        * dg_Nlx2Mat_Samples;
                    dg_Nlx2Mat_SamplesUnits = 'V';
                end
            end
        end
        save(matfilename, 'dg_Nlx2Mat_Timestamps', ...
            'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits', '-v7.3');
    case '.NEV'
        [dg_Nlx2Mat_Timestamps, dg_Nlx2Mat_TTL, dg_Nlx2Mat_EventStrings] ...
            = dg_readEvents(pathname);
        save(matfilename, 'dg_Nlx2Mat_Timestamps', ...
            'dg_Nlx2Mat_TTL', 'dg_Nlx2Mat_EventStrings');
end