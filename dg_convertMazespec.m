function mazespec2 = dg_convertMazespec(mazespec, evtIDs);
% Converts a <mazespec> as returned by lfp_findMazespec to a <mazespec2> as
% needed by lfp_linearizeRodentTracker2 (i.e. containing median pixel
% coordinates), based on the measurements in "T maze coordinates for
% linearizer.doc" and the assumption that the scale factor (cmperpixel) is
% independent of position and orientation.
% The transformation formulae are thus:
%   xpix = X0 + xcoord/cmperpixel
%   ypix = Y3 - ycoord/cmperpixel
% where we must further calculate
%   Y3 = (Y0 + Y1)/2 + ((LTurnOffCoord + RTurnOffCoord)/2)/cmperpixel
% The coordinates missing from "T maze coordinates for linearizer.doc" are
% considered to be the centerline of the relevant section of T maze.  The
% fields returned are 'evtIDs' and 'medians'; no 'n' field is returned.
%NOTES
% Requires the following mazespec fields: 'cmperpixel', 'Y0', 'Y1', 'X0'.
% WARNING: the assumption re: scale factor is not very good!  In reality,
% scale factor varies by about 2/3 from the close end to the far end of the
% maze.

%$Rev: 42 $
%$Date: 2009-10-10 01:00:18 -0400 (Sat, 10 Oct 2009) $
%$Author: dgibson $

if nargin < 2
    evtIDs = { 13 6 [23 31 38 21 22 20] 9 14 15 16 7 8 17 18 };
end

% Cartesian coordinates from the notorious "T maze coordinates for
% linearizer.doc".  Event IDs in col 1, x position in col 2, y position in
% col 3.
mazecoords = [
    13  33.8    NaN
    6   51.2    NaN
    23  68.8    NaN
    31  68.8    NaN
    38  68.8    NaN
    21  68.8    NaN
    22  68.8    NaN
    20  68.8    NaN
    9   86.6    NaN
    14  109.0   NaN
    15  NaN     25.1
    16  NaN     48.7
    7   NaN     14.3
    8   NaN     59.2
    17  NaN     3.6
    18  NaN     69.9
    ];
TurnOn = 14;
RTurnOff = 15;
LTurnOff = 16;
TurnOnCoord = mazecoords(mazecoords(:,1) == TurnOn, 2);
RTurnOffCoord = mazecoords(mazecoords(:,1) == RTurnOff, 3);
LTurnOffCoord = mazecoords(mazecoords(:,1) == LTurnOff, 3);

% From the incomparable "T maze markers for rat from Hu Dan.doc" as
% further annotated Friday, January 16, 2009 7:31:25 PM by DG:
turn2crossbar = 10.2;
crossbarwidth = 7.7;

cmperpixel = mazespec.cmperpixel;
stemCLpix = (mazespec.Y0 + mazespec.Y1)/2;
crossbarCLpix = mazespec.X0 + ...
    (TurnOnCoord + turn2crossbar + crossbarwidth/2)/cmperpixel;
Y3 = stemCLpix + ...
    ((LTurnOffCoord + RTurnOffCoord)/2)/cmperpixel;
mazespec2.evtIDs = evtIDs;
for e = 1:length(evtIDs)
    row = find(mazecoords(:,1) == evtIDs{e}(1));
    if isnan(mazecoords(row, 3))
        mazespec2.medians(e,:) = ...
            [mazespec.X0 + mazecoords(row, 2)/cmperpixel, stemCLpix];
    else
        mazespec2.medians(e,:) = ...
            [crossbarCLpix, Y3 - mazecoords(row, 3)/cmperpixel];
    end
end

