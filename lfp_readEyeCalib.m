function eyecalibs = lfp_readEyeCalib
%eyecalibs = lfp_readEyeCalib
%   If there is an eye calibration values file in the root directory for
%   the current session, return its contents.  Otherwise, return [].

%$Rev: 287 $
%$Date: 2012-10-26 20:14:45 -0400 (Fri, 26 Oct 2012) $
%$Author: dgibson $

lfp_declareGlobals;

eyecalibs = [];

eyecalibfile = lfp_getEyeCalibFileAbs;
if exist(eyecalibfile) == 2
    try load(eyecalibfile, '-mat');
    catch
        warning('lfp_readEyeCalib:badfile', ...
            'Failed to load eye calib file %s', ...
            eyecalibfile );
        return
    end
    lfp_log(sprintf('Loaded eye calibrations from %s', ...
        eyecalibfile));
end