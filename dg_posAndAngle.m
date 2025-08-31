function [x, y, theta] = dg_posAndAngle(targ, colors)
%INPUT
% <targ> is a targets array as returned by Nlx2MatVT_v3.
% <colors> is optional; default is to ignore target colors.
% GOOD INTENTIONS - NOT IMPLEMENTED
% If <colors> is
% given, it is a 2 x 3 bit array where the first row represents the color
% of the LED on the animal's left, and the second row is the color on the
% right.
%OUTPUTS
% All the outputs are column vectors of length size(targ,2).
% <x>, <y> are the positions of the midpoint of the first two targets as
% recorded in <targ>.  theta is the angle in radians of the two targets
% with respect to vertical (12 o'clock), positive to the left
% (counter-clockwise).  Thus theta is the angle with respect to the
% horizontal (standard engineering representation) of the rat's head
% direction when the LEDs are mounted to the left and right of the
% headstage.  Note that if <colors> is not specified, then <theta> ranges
% only -pi/2 to +pi/2, whereas if <colors> is given, <theta> ranges from
% -pi to +pi.

if nargin < 2
    usecolors = false;
else
    usecolors = true;
end

x = NaN(size(targ,2), 1);
y = NaN(size(targ,2), 1);
theta = NaN(size(targ,2), 1);

for frame = 1:size(targ, 2)
    for t = 1:2
        [r(t), g(t), b(t), inten(t), thisframe_x(t), thisframe_y(t)] = ...
            vt_bitfield_decode(targ(t, frame));
    end
    x(frame) = mean(thisframe_x);
    y(frame) = mean(thisframe_y);
    if usecolors
        error('colors not implemented');
    else
        % The lower y value is above, and so is arbitrarily assigned as the
        % left LED; if both LEDs are at the same y, the first target is
        % assigned to be left.
        [c, left] = min(thisframe_y);
        if length(left) > 1
            left = 1;
        end
    end
    right = mod(left,2) + 1;    % i.e., the other one
    theta(frame) = atan( ...
        (thisframe_x(left) - thisframe_x(right)) / ...
        (thisframe_y(left) - thisframe_y(right)) );
end


