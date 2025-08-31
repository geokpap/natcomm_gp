function dg_rbcolormapdemo
% This is a demonstration of the red-blue colormap contained in
% 'dg_rbcolormap.mat'.  'dg_rbcolormap.mat' contains just one variable
% named 'rbcolormap'.  It is a 64-color colormap useful for creating
% pseudocolor images (e.g. using 'imagesc') that look similar both to
% trichromats and to red-green colorblind people.  It prints well in native
% RGB color mode and when transformed to CMYK color mode. The resolution of
% extreme values is ever so slightly better in RGB mode when printed.

%$Rev: 202 $
%$Date: 2014-07-24 15:15:09 -0400 (Thu, 24 Jul 2014) $
%$Author: dgibson $

load('dg_rbcolormap.mat');
hF = figure;
imagesc(peaks);
caxis([-8 8]);
colorbar;
set(hF, 'colormap', rbcolormap);
