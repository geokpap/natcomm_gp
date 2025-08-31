function [mazespec, numborderline, numtotal] = dg_measureTMaze(filename)
% Measures T-maze parameters without reference to event data.  This is done
% by fitting a horizontal line to the stem of the T and a vertical line to
% the crossbar, and computing the length of the crossbar as the distance
% between the 2nd percentile and 98th percentile of all points that fall in
% the crossbar portion of the video mask.  cmperpixel is calculated from
% the crossbar length.  All of the other parameters are then calculated
% from the fitted centerlines and cmperpixel.
%INPUTS
% <filename> is a Neuralynx Video Tracker file name.
%OUTPUTS
% <mazespec> is a lfp_linearizeRodentTracker-compatible T-maze layout
%   specification.  Fields: 'X1', 'X2', 'Y0', 'Y1', 'cmperpixel', 'StimOn'
%   (x coordinate), 'RGoal' (y coordinate), 'LGoal', 'trackwidth' (width of
%   of central (1-2*tailwidth) of actual trace).
% <numborderline>: number of video tracker points that are within 10 pixels
%   of the video mask
% <numtotal>: total number of video tracker points

mazespec = [];
[x, y] = Nlx2MatVT_v3(filename, [0 1 1 0 0 0], 0, 1);
x = reshape(x, [], 1);
y = reshape(y, [], 1);

% Values from lfp_add Rev: 69, validated by overlaying s18\acq15,
% s36\acq10, d19\acq05, a98\acq28, h01\h01acq10, with X0 increased to mask
% out the start box "blob" if present:
maskParameters = struct('maskLabel','c02-based default mask',...
    'X0',200,'X2',480,'X3',625,'Y0',140,'Y1',280,'Y4',30,'Y3',390);

borderwidth = 10;
borderline = ( ...
    ( x>maskParameters.X0 - borderwidth & ...
    x<=maskParameters.X2 & ...
    y>maskParameters.Y0 - borderwidth & ...
    y<maskParameters.Y1 + borderwidth ) | ...
    ( x>maskParameters.X2 & ...
    x<maskParameters.X3 + borderwidth & ...
    y>maskParameters.Y4 - borderwidth & ...
    y<maskParameters.Y3 + borderwidth ) ...
    ) & ( ...
    ( x<maskParameters.X0 + borderwidth & ...
    x<=maskParameters.X2 & ...
    y<maskParameters.Y0 + borderwidth & ...
    y<maskParameters.Y1 - borderwidth ) | ...
    ( x>maskParameters.X2 & ...
    x>maskParameters.X3 - borderwidth & ...
    y<maskParameters.Y4 + borderwidth & ...
    y>maskParameters.Y3 - borderwidth ));
numborderline = sum(borderline);
numtotal = numel(x);

onstem = (x>maskParameters.X0 & ...
    x<=maskParameters.X2 & ...
    y>maskParameters.Y0 & ...
    y<maskParameters.Y1);

oncrossbar = (x>maskParameters.X2 & ...
    x<maskParameters.X3 & ...
    y>maskParameters.Y4 & ...
    y<maskParameters.Y3);

% Find centerlines and crossbar ends:
crossbarCL = median(x(oncrossbar));
stemCL = median(y(onstem));
p = 2;
crossbarstart = prctile(y(oncrossbar), p);
crossbarend = prctile(y(oncrossbar), 100-p);

% Measurements in cm from "T maze coordinates for linearizer.doc":
crossbarlength = 73.5;
x0_x1 = 20.2;
x0_turn = 109;
x0_stimon = 68.8;

% Measurements in cm from "T maze markers for rat from Hu Dan.doc":
trackwidth = 7.7;
turn_x2 = 10.2;
goal_end = 3.5;

cmperpixel = crossbarlength / (crossbarend - crossbarstart);
mazespec.cmperpixel = cmperpixel;
mazespec.trackwidth = trackwidth/cmperpixel;
mazespec.X0 = crossbarCL - (trackwidth/2 + turn_x2 + x0_turn)/cmperpixel;
mazespec.X1 = mazespec.X0 + x0_x1/cmperpixel;
mazespec.X2 = crossbarCL - trackwidth/(2*cmperpixel);
mazespec.X3 = crossbarCL + trackwidth/(2*cmperpixel);
mazespec.Y0 = stemCL - trackwidth/(2*cmperpixel);
mazespec.Y1 = stemCL + trackwidth/(2*cmperpixel);
mazespec.StimOn = mazespec.X0 + x0_stimon/cmperpixel;
mazespec.RGoal = crossbarend - goal_end/cmperpixel;
mazespec.LGoal = crossbarstart + goal_end/cmperpixel;

