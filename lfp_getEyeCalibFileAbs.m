function eyecalibfilename = lfp_getEyeCalibFileAbs
%eyecalibfilename = lfp_getEyeCalibFilename
%   Returns the absolute pathname of the eye calibration file for the
%   current session.

%   The eye calib file resides in the session root directory.  We assume
%   that lfp_DataDir is pointing at the current session or session fragment
%   directory.  Starting at lfp_DataDir, we apply the following iterative
%   formula:  if the directory contains an lfp_FragmentFiles.mat, then do
%   another iteration looking at the parent directory; otherwise the
%   current directory is assumed to be the root.

%$Rev: 175 $
%$Date: 2010-09-30 19:29:27 -0400 (Thu, 30 Sep 2010) $
%$Author: dgibson $

global lfp_DataDir;

mydirectory = lfp_DataDir;
sessionroot = [];
while isempty(sessionroot)
    fraginfofilepath = fullfile(mydirectory, 'lfp_FragmentFiles.mat');
    [parent,name,ext] = fileparts(mydirectory);
    if exist(fraginfofilepath) == 2
        mydirectory = parent;
    else
        if isempty(name)
            % Oops, hit the file system root
            sessionroot = 0;
        else
            sessionroot = mydirectory;
        end
    end
end
if sessionroot
    eyecalibfilename = fullfile(sessionroot, 'lfp_ManualEyeCalib.Mat');
else
    error('lfp_getEyeCalibFileAbs:noroot', 'Could not find session root');
end
