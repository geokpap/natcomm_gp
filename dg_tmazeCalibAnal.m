function mazespec = dg_tmazeCalibAnal(calib)
%mazespec = dg_tmazeCalibAnal(calib)
% Reports analysis results on a set of T maze calibrations.
%mazespec = dg_tmazeCalibAnal(calib)
% Converts the set of T maze calibrations to a mazespec, does not display
% analysis report.

%$Rev: 42 $
%$Date: 2009-10-10 01:00:18 -0400 (Sat, 10 Oct 2009) $
%$Author: dgibson $

showstuff = true;
if nargout > 0
    showstuff = false;
end

x2_3 = 119.5;   % point 2 to point 3, cm
x2_13 = 20.2 + 48.5;   % point 2 to point 13, cm
x2_start = 20.2;    % point 2 to front of start box, cm
y6_8 = 73.8;    % point 6 to point 8, cm
y20_21 = 66;    % point 20 (LGoal) to point 21 (RGoal), cm
overhang = 10;   % number of pixels by which markers exceed track width

% scale is calculated from the long distances, and so is averaged over
% changes due to keystoning
xscale1 = x2_3 / (calib.x(3) - calib.x(2));
xscale2 = x2_3 / (calib.x(4) - calib.x(1));
xscale = x2_3 / (mean(calib.x([3 4 5 7])) - mean([calib.x([1 2])]));
pixelprecision_xs = ((calib.x(3) - calib.x(2)) + 1) ...
    / (calib.x(3) - calib.x(2)) - 1;
yscale1 = y6_8 / (mean(calib.y(7:8)) - mean(calib.y(5:6)));
yscale2 = y20_21 / (calib.y(21) - calib.y(20));
% Twice as many points go into yscale1 as yscale2, so it is weighted more
% heavily in the average:
yscale = (2*yscale1 + yscale2)/3;
pixelprecision_ys = ((calib.y(8) - calib.y(6)) + 1) ...
    / (calib.y(8) - calib.y(6)) - 1;
% yscale is calculated over only about half the distance of the xscale, so
% should be given less weight than xscale in average:
mazespec.cmperpixel = (yscale + 2*xscale)/3;
trackwidth = mean([
    calib.y(1)-calib.y(2)
    calib.y(4)-calib.y(3)
    calib.x(6)-calib.x(5)
    calib.x(8)-calib.x(7)
    ]);
x2 = mean(calib.x([3 4 5 7]));
mazespec.X0 = x2 - x2_3/mazespec.cmperpixel;
mazespec.X1 = mazespec.X0 + x2_start / mazespec.cmperpixel;
stemcenterline = mean(calib.y(1:4));
mazespec.Y0 = stemcenterline - trackwidth/2;
mazespec.Y1 = stemcenterline + trackwidth/2;

% keystoning: a simple linear approximation in the x direction for camera
% over the base or the crossbar, based on track width in y direction. Since
% there is no width measurement at the middle of the stem, we cannot
% calculate keystoning for a camera over the center of the maze.
ywidth_lowX = calib.y(1) - calib.y(2);
ywidth_highX = calib.y(4) - calib.y(3);
keystoning = ywidth_highX / ywidth_lowX - 1;
pixelprecision_k = (ywidth_lowX + 1) / ywidth_lowX - 1;

% angleY defined as the angle between the crossbar and the Y axis, angleX
% is the angle between the stem and the X axis.  Positive =
% counterclockwise.
dpr = 360 / (2*pi);
angleY1 = dpr * asin((calib.x(8) - calib.x(6))/(calib.y(8) - calib.y(6)));
angleY2 = dpr * asin((calib.x(7) - calib.x(5))/(calib.y(7) - calib.y(5)));
ypixelangle = dpr * asin(1/(calib.y(7) - calib.y(5)));
angleX1 = dpr * asin((calib.y(2) - calib.y(3))/(calib.x(3) - calib.x(2)));
angleX2 = dpr * asin((calib.y(1) - calib.y(4))/(calib.x(4) - calib.x(1)));
xpixelangle = dpr * asin(1/(calib.x(4) - calib.x(1)));
angleX = mean([angleX1 angleX2]);
angleY = mean([angleY1 angleY2]);

if showstuff
    disp(sprintf('\nScale:'));
    disp(sprintf('       1     2   avg'));
    disp(sprintf('x  %5.3f %5.3f %5.3f', xscale1, xscale2, xscale));
    disp(sprintf('y  %5.3f %5.3f %5.3f', yscale1, yscale2, yscale));
    disp(sprintf('Grand Average  %5.3f', mazespec.cmperpixel));
    disp(sprintf('Measured Scale Precision: %.1f%%', 100*max(abs( ...
        ([xscale1 xscale2 yscale1 yscale2] - mazespec.cmperpixel) ) ...
        / mazespec.cmperpixel )));
    disp(sprintf('x scale pixel precision %.1f%%', 100*pixelprecision_xs));
    disp(sprintf('y scale pixel precision %.1f%%', 100*pixelprecision_ys));
    disp(sprintf('Keystoning: %.1f%%', 100*keystoning));
    disp(sprintf('Keystoning pixel precision: %.1f%%', 100*pixelprecision_k));
    disp(sprintf('\nRotation Angle:'));
    disp(sprintf('       1     2   avg'));
    disp(sprintf('x  %5.3f %5.3f %5.3f', angleX1, angleX2, angleX));
    disp(sprintf('y  %5.3f %5.3f %5.3f', angleY1, angleY2, angleY));
    disp(sprintf('Grand Average  %5.3f', mean([angleX angleY])));
    disp(sprintf('x pixel angle precision: %5.3f', xpixelangle));
    disp(sprintf('y pixel angle precision: %5.3f', ypixelangle));
    disp('');
    
    % [9 10 13 20 21] are used only for plotting
    for k = [9 10 13]
        if isnan(calib.y(k))
            calib.y(k) = (calib.y(1) + calib.y(2)) / 2;
        end
    end
    for k = 20:21
        if isnan(calib.x(k))
            calib.x(k) = (calib.x(5) + calib.x(6)) / 2;
        end
    end
    figure;
    plot( calib.x([1 2 3 5 6 8 7 4 1]), ...
        calib.y([1 2 3 5 6 8 7 4 1]), 'Marker', '.' );
    hold on;
    plot(calib.x(1), calib.y(1), ...
        'g', 'Marker', '.');
    plot(calib.x(4), calib.y(4), ...
        'r', 'Marker', '.');
    plot(calib.x(20), calib.y(20), ...
        'g', 'Marker', '.');
    plot(calib.x(21), calib.y(21), ...
        'r', 'Marker', '.');
    plot(calib.x(13), calib.y(13), ...
        'm', 'Marker', '.');
    plot(calib.x(9), calib.y(9), ...
        'm', 'Marker', '.');
    plot(calib.x(10), calib.y(10), ...
        'm', 'Marker', '.');
    calc13 = mazespec.X0 + x2_13 / mazespec.cmperpixel;
    plot([calc13 calc13], ...
        [mazespec.Y0 - overhang mazespec.Y1 + overhang], 'k');
    calcstart = mazespec.X0 + x2_start / mazespec.cmperpixel;
    plot([calcstart calcstart], ...
        [mazespec.Y0 - overhang mazespec.Y1 + overhang], 'k');
    title(calib.calibname, 'Interpreter', 'none');
    grid on;
    axis equal;
    xlim([0 650]);
    ylim([0 450]);
    axis ij;
end
    
