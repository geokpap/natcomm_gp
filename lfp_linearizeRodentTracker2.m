function [result, units] = lfp_linearizeRodentTracker2(fnums, mazespec2, varargin)
%lfp_linearizeRodentTracker2 is intended for use with lfp_createWave.
% result = lfp_linearizeRodentTracker2(fnums, mazespec2)
%   Convert T-maze (x,y) to linear distance in cm from starting position.
%   Starting position is defined as the Out of Start photobeam.
%   Conversion is done simple-mindedly by using x position relative to
%   start until the animal reaches the "turn box", which is bounded by the
%   Turn On photobeam and the two Turn Off photobeams.  Within the turn
%   box, distance is defined as the sum of the x and y distances from the
%   point where the Turn On line was crossed, and is normalized so that the
%   total distance from that point to the point where Turn Off is crossed
%   is equal to the distance measured along the centerline of the runway.
%   Position after crossing Turn Off is just the y position relative to the
%   Turn Off line.  Conversion from pixels to cm is done in a piecewise
%   linear fashion based on the available information in <mazespec2>.  In
%   the case of positions that fall outside of the range of the photobeams
%   in <mazespec2>, the same scale factor is used as for positions between
%   the closest pair of photobeams (i.e. the first scale factor is used for
%   positions before the first photobeam).
% <fnums> should contain the x channel first and the y channel second.
%   For each trial, the data are assumed to start at 2 s before
%   lfp_NominalTrialStart and to end at lfp_NominalTrialEnd.
% <mazespec2> is a structure of the type returned by lfp_measureTMaze2.
%   Fields: 'evtIDs', 'medians' (median position over trials), 'n'. 
%OPTIONS
% 'mazecoords', mazecoords - <mazecoords> overrides the default values for
%   the physical coordinates of the photobeams.  It is a 3 column array
%   with event IDs in col 1, x position in col. 2, and y position in col.
%   3.  Note that each photobeam has a valid position in only one of cols.
%   2 or 3, and NaN in the other one.  Coordinates are "ij", not Cartesian.
% 'minN', minN - overrides the default value of <minN>, which is 3.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

TurnOn = 14;
RTurnOff = 15;
LTurnOff = 16;

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
mazecoords(:,3) = 100 - mazecoords(:,3);    % convert from Cartesian

if ~isequal(size(fnums), [1 2])
    error('lfp_velo:badfilenums', ...
        '<fnums> must have exactly 2 elements on one row' );
end
argnum = 1;
mazespec = [];
minN = 3;
while argnum <= length(varargin)
    if strcmp(class(varargin{argnum}), 'char')
        switch varargin{argnum}
            case 'mazecoords'
                argnum = argnum + 1;
                mazecoords = varargin{argnum};
            case 'minN'
                argnum = argnum + 1;
                minN = varargin{argnum};
            otherwise
                error('lfp_linearizeRodentTracker2:badoption', ...
                    ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
        end
    end
    argnum = argnum + 1;
end

% Find which photobeam events we have pixel coordinates for in mazespec2
isusable = false(size(mazecoords(:,1)));
for k = 1:length(mazespec2.evtIDs)
    isinmazespec2k = ismember(mazecoords(:,1), mazespec2.evtIDs{k});
    if any(isinmazespec2k) && mazespec2.n(k) < minN
        warning('lfp_linearizeRodentTracker2:toofew', ...
            'There are two few instances of evtID %s; ignoring', ...
            mat2str(mazespec2.evtIDs{k}) );
    else
        isusable = isusable | isinmazespec2k;
    end
end
evtIDs2use = mazecoords(isusable,1);
% Construct <conversiontable> with evt IDs and physical coords as in
% mazecoords, pixel coords in the non-NaN dimension in column 4.
conversiontable = sortrows( ...
    mazecoords(ismember(mazecoords(:,1), evtIDs2use), :), [2 3] );
is_y = isnan(conversiontable(:, 2));
for row = 1:size(conversiontable, 1)
    % Assume each event ID appears in only one element of mazespec2.evtIDs
    for evtIDidx = 1:length(mazespec2.evtIDs)
        if ismember(conversiontable(row,1), mazespec2.evtIDs{evtIDidx})
            % This MUST be true for some evtIDidx, since conversiontable
            % was constructed from evtIDs2use
            break
        end
    end
    conversiontable(row, 4) = mazespec2.medians( evtIDidx, ...
        is_y(row) + 1 );
end
% Check to be sure that col. 4 is sorted the same as cols. 2,3:
xpixels = conversiontable(~is_y, 4);
ypixels = conversiontable(is_y, 4);
if any(xpixels(1:end-1) > xpixels(2:end)) || ...
        any(ypixels(1:end-1) > ypixels(2:end))
    error('lfp_linearizeRodentTracker2:mazespec2', ...
        'Disordered conversion table: %s', mat2str(conversiontable) );
end

% Mark some important rows:
first_y = find(is_y);
first_y = first_y(1);
TurnOnRow = find(conversiontable(:,1) == TurnOn);
RTurnOffRow = find(conversiontable(:,1) == RTurnOff);
LTurnOffRow = find(conversiontable(:,1) == LTurnOff);

% Calculate turn box parameters:
turnboxYcm = ( conversiontable(RTurnOffRow, 3) - ...
    conversiontable(LTurnOffRow, 3) )/2;
turnboxXcm = 14;    % sadly, this cannot be calculated from mazecoords
turnboxtotcm = turnboxXcm + turnboxYcm;
turnctrYcoord = ( conversiontable(LTurnOffRow, 3) + ...
    conversiontable(RTurnOffRow, 3) ) / 2;

result = NaN(size(lfp_Samples{fnums(1)}));
units = 'cm';

% We assume here that conversiontable is sorted with all x-coordinate rows
% first, which happens to be the way that sortrows does it (see code above)
% in Matlab R2007b, but this is not actually specified in the Matlab help.

% Convert zone-by-zone.  The zones at the ends of the stem and arms are
% special because they extend indefinitely; the "turn box" is also special
% (see header comments). 

% left of stem of T
row = 2;
numcm = conversiontable(row, 2) - conversiontable(row-1, 2);
fnumidx = 1;
cmperpixel = numcm/(conversiontable(row, 4) - conversiontable(row-1, 4));
inzone = lfp_Samples{fnums(fnumidx)} < conversiontable(row, 4);
if ~isempty(inzone)
    samples = cmperpixel * ...
        (lfp_Samples{fnums(fnumidx)}(inzone) - conversiontable(row-1, 4));
    result(inzone) = samples + ...
        conversiontable(row-1, 2) - conversiontable(1, 2);
end

% stem of T
for row = 3:TurnOnRow
    [samples, inzone] = linearize(row, conversiontable, fnums, result);
    if ~isempty(inzone)
        result(inzone) = samples + ...
            conversiontable(row-1, 2) - conversiontable(1, 2);
    end
end

% On the L arm, where direction of travel is opposite to direction of Y
% axis, travelling a negative distance along the Y axis translates into a
% positive linearized distance from the start.  Therefore, the sign of
% <samples> must be reversed when linearizing.  Similarly, the distance
% between the reference photobeam on row number (<row>-1) and LTurnOff will
% also be negative, and must also therefore be have its sign reversed.

% above L arm
row = first_y+1;
numcm = conversiontable(row, 3) - conversiontable(row-1, 3);
fnumidx = 2;
cmperpixel = numcm/(conversiontable(row, 4) - conversiontable(row-1, 4));
inzone = lfp_Samples{fnums(fnumidx)} < conversiontable(row, 4) ...
    & isnan(result);
if ~isempty(inzone)
    samples = cmperpixel * ...
        (lfp_Samples{fnums(fnumidx)}(inzone) - conversiontable(row-1, 4));
    result(inzone) = -samples ...
        + conversiontable(TurnOnRow, 2) - conversiontable(1, 2) ...
        + turnboxtotcm ...
        - conversiontable(row-1, 3) + conversiontable(LTurnOffRow, 3);
end

% L arm
for row = (first_y+2) : LTurnOffRow
    [samples, inzone] = linearize(row, conversiontable, fnums, result);
    if ~isempty(inzone)
        result(inzone) = -samples ...
            + conversiontable(TurnOnRow, 2) - conversiontable(1, 2) ...
            + turnboxtotcm ...
            - conversiontable(row-1, 3) + conversiontable(LTurnOffRow, 3);
    end
end

% R arm
% No inversion problems here.  Linear distance is just <samples> plus the
% total distance from start to the reference photobeam on row number
% (<row>-1), which in turn is start-to-TurnOn plus <turnboxtotcm> plus
% RTurnOff-to-ref.
for row = (RTurnOffRow+1):(size(conversiontable,1) - 1)
    [samples, inzone] = linearize(row, conversiontable, fnums, result);
    if ~isempty(inzone)
        result(inzone) = samples ...
            + conversiontable(TurnOnRow, 2) - conversiontable(1, 2) ...
            + turnboxtotcm ...
            + conversiontable(row-1, 3) - conversiontable(RTurnOffRow, 3);
    end
end        

% below R arm
row = size(conversiontable,1);
numcm = conversiontable(row, 3) - conversiontable(row-1, 3);
fnumidx = 2;
cmperpixel = numcm/(conversiontable(row, 4) - conversiontable(row-1, 4));
inzone = lfp_Samples{fnums(fnumidx)} >= conversiontable(row-1, 4) ...
    & isnan(result);
if ~isempty(inzone)
    samples = cmperpixel * ...
        (lfp_Samples{fnums(fnumidx)}(inzone) - conversiontable(row-1, 4));
    result(inzone) = samples ...
        + conversiontable(TurnOnRow, 2) - conversiontable(1, 2) ...
        + turnboxtotcm ...
        + conversiontable(row-1, 3) - conversiontable(RTurnOffRow, 3);
end

% turn box must be done trial-by-trial because linear scale factor varies.
% Each "zone" to be scaled thus has both temporal and spatial criteria.
for trial = 1:length(lfp_SelectedTrials)
    % First we find the sample numbers and pixel coordinates of the entry
    % to the turn box and exit from the turn box.  The pixel coordinates
    % that we need to measure are the ones in the direction perpendicular
    % to the direction of motion, which presumably vary slowly enough so
    % that temporal slippage between the video tracker and event recorder
    % is negligible.  We assume that the pixel coordinates parallel to the
    % direction of motion are exactly as specified in mazespec2.
    firstevtidx = lfp_TrialIndex(trial,1);
    lastevtidx = lfp_TrialIndex(trial,2);
    turnevtidx = find(lfp_Events(firstevtidx:lastevtidx,2) == TurnOn) ...
        + firstevtidx - 1;
    turnidx = lfp_time2index(lfp_Events(turnevtidx,1));
    turnoffevtidx = find(ismember( ...
        lfp_Events(firstevtidx:lastevtidx,2), [RTurnOff LTurnOff] )) ...
        + firstevtidx - 1;
    turnoffidx = lfp_time2index(lfp_Events(turnoffevtidx,1));
    turnX = conversiontable(TurnOnRow,4);
    turnY = lfp_Samples{fnums(2)}(turnidx);
    turnoffX = lfp_Samples{fnums(1)}(turnoffidx);
    if lfp_Events(turnoffevtidx,2) == RTurnOff
        turnoffY = conversiontable(RTurnOffRow,4);
    else
        turnoffY = conversiontable(LTurnOffRow,4);
    end
    pixelsX = turnoffX - turnX;
    pixelsY = turnoffY - turnY; % Note that this can be negative
    totpixels = pixelsX + abs(pixelsY); % This cannot be negative
    cmperpixel = turnboxtotcm/totpixels; % This cannot be negative
    firstidx = lfp_TrialIndex(trial, 3);
    lastidx = lfp_TrialIndex(trial, 4);
    inzone = find( ...
        lfp_Samples{fnums(1)}(firstidx:lastidx) >= ...
        conversiontable(TurnOnRow, 4) ...
        & lfp_Samples{fnums(2)}(firstidx:lastidx) >= ...
        conversiontable(LTurnOffRow, 4) ...
        & lfp_Samples{fnums(2)}(firstidx:lastidx) < ...
        conversiontable(RTurnOffRow, 4) ) + firstidx - 1;
    if ~isempty(inzone)
        result(inzone) = ...
            conversiontable(TurnOnRow, 2) - conversiontable(1, 2) + ...
            (lfp_Samples{fnums(1)}(inzone) - turnX) * cmperpixel + ...
            (lfp_Samples{fnums(2)}(inzone) - turnY) * cmperpixel ...
            * sign(pixelsY);
    end
end
for trial = 1:length(lfp_SelectedTrials)
    if any(isnan(result(lfp_TrialIndex(trial,3):lfp_TrialIndex(trial,4))));
        nanidx = find(isnan( ...
            result(lfp_TrialIndex(trial,3):lfp_TrialIndex(trial,4)) )) ...
            + lfp_TrialIndex(trial,3) - 1;
        warning('lfp_linearizeRodentTracker2:NaN', ...
            'Trial %d contains unassigned values between %.6f and %.6f s', ...
            trial, lfp_index2time(nanidx(1)), lfp_index2time(nanidx(end)));
    end
end

end


function [samples, inzone] = linearize(row, conversiontable, fnums, result)
% Returns <samples> scaled to cm relative to the photobeam position on
% row number (<row>-1) of <conversiontable>.  <inzone> is true for each
% element of lfp_Samples that is in the zone from (<row>-1) to <row> of
% <conversiontable>.
global lfp_Samples
if isnan(conversiontable(row,2))
    numcm = conversiontable(row, 3) - conversiontable(row-1, 3);
    fnumidx = 2;
else
    numcm = conversiontable(row, 2) - conversiontable(row-1, 2);
    fnumidx = 1;
end
cmperpixel = numcm/(conversiontable(row, 4) - conversiontable(row-1, 4));
inzone = lfp_Samples{fnums(fnumidx)} < conversiontable(row, 4) & ...
    lfp_Samples{fnums(fnumidx)} >= conversiontable(row-1, 4) ...
    & isnan(result);
if isempty(inzone)
    samples = [];
else
    samples = cmperpixel * ...
        (lfp_Samples{fnums(fnumidx)}(inzone) - conversiontable(row-1, 4));
end
end
