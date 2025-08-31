function [result, units] = lfp_linearizeRodentTracker(fnums, varargin)
%lfp_linearizeRodentTracker is intended for use with lfp_createWave.
% result = lfp_linearizeRodentTracker(fnums)
%   Convert T-maze (x,y) to linear distance in cm from starting position.
%   Starting position is defined as the front edge of the starting box.
%   Conversion is done simple-mindedly by using x position relative to
%   start until the animal reaches the left edge of the T; thereafter
%   position is calculated as the maximal x distance plus the distance in
%   the y direction from the center of the stem of the T.  <fnums> should
%   contain the x channel first and the y channel second.  For each trial,
%   the data are assumed to start at 2 s before lfp_NominalTrialStart and
%   to end at lfp_NominalTrialEnd.
%OPTIONS
% 'mazespec', <mazespec>
%   <mazespec> is a structure that overrides the default maze layout
%   parameters.  Its field names are 'X1', 'X2', 'Y0', 'Y1', 'cmperpixel',
%   and they must ALL be provided.  Default values were informally
%   estimated from s23acq07. See code for details.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~isequal(size(fnums), [1 2])
    error('lfp_velo:badfilenums', ...
        '<fnums> must have exactly 2 elements on one row' );
end
argnum = 1;
mazespec = [];
while argnum <= length(varargin)
    if strcmp(class(varargin{argnum}), 'char')
        switch varargin{argnum}
            case 'mazespec'
                argnum = argnum + 1;
                mazespec = varargin{argnum};
        end
    end
    argnum = argnum + 1;
end


X1 = 110;   % R (front) end of start box
X2 = 535;   % L edge of arms
Y0 = 220;   % top edge of stem
Y1 = 200;   % bottom edge of stem
% "C:\Documents and Settings\dgibsonGraybiel\My Documents\Graybiel Lab\LFP
% Analysis\Bill Paper 2005\Supp Fig 3 Re.doc" states that the rat runs
% about 110 cm along the stem and about 65 cm along the arm, which led me
% to settle on 0.27 cm/pixel.
cmperpixel = 0.27;

if ~isempty(mazespec)
    try
        X1 = mazespec.X1;   % R (front) end of start box
        X2 = mazespec.X2;   % L edge of arms
        Y0 = mazespec.Y0;   % top edge of stem
        Y1 = mazespec.Y1;   % bottom edge of stem
        cmperpixel = mazespec.cmperpixel;
    catch
        error('lfp_linearizeRodentTracker:mazespec', ...
            'Problem using the mazespec structure.' );
    end
end


result = zeros(size(lfp_Samples{fnums(1)}));
for trial = 1:size(lfp_TrialIndex,1)
    samples = (lfp_TrialIndex(trial,3) - round(2/lfp_SamplePeriod)) ...
        : lfp_TrialIndex(trial,4);
    xdata = lfp_Samples{fnums(1)}(samples);
    ydata = lfp_Samples{fnums(2)}(samples);
    inarms = find(xdata > X2);
    if isempty(inarms)
        result(samples) = [xdata - X1];
    else
        startarmidx = inarms(1);
        endstemidx = inarms(1) - 1;
        result(samples) = [xdata(1:endstemidx) - X1,  ...
            abs(ydata(startarmidx:end) - (Y0+Y1)/2) + X2 - X1 ];
    end
end
result = cmperpixel * result;
units = 'cm';

