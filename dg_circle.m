function dg_circle(x,y,r,s)
%dg_circle(x,y,r,s)
%   <x>, <y>: position
%   <r>: radius
%   <s>: linespec

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

angle = linspace(0, 2*pi, 101);
x_coords = x + r * cos(angle);
y_coords = y + r * sin(angle);
plot(x_coords, y_coords, s);