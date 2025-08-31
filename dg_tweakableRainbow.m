function mycmap = dg_tweakableRainbow
%dg_tweakableRainbow produces a color map based on the idea of 'jet', but
%   with an effort to prevent brightness from distorting the perception of
%   the quantitative scale (e.g. when pure yellow sticks out more than red,
%   but corresponds to a lower value).  Rather than trying to provide
%   meaningful input parameters and defalut values, I have simply tweaked
%   the params a, b, c, d in situ, and you can do the same.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

blues = [ zeros(8,2) (9:16)'/16 ];
rainbow = hsv2rgb([ (2*(1 - (1:56)/56)/3)' ones(56, 2) ]);
a = 2;
normrainbow = rainbow ./ repmat((sum(rainbow.^a,2)).^(1/a), 1, 3);
mycmap = [ blues; normrainbow ];
b = .7; % yellow emphasis
c = -.5; % green emphasis
d = -.2; % blue emphasis
fudge = ones(64,1);
fudge = fudge + [ zeros(40,1); b*dpss(24,3,1) ];
fudge = fudge + [ zeros(8,1); c*dpss(56,4,1) ];
fudge = fudge + [ d*dpss(16,4,1); zeros(48,1) ];
mycmap = mycmap .* repmat(fudge, 1, 3);
mycmap = min(mycmap, ones(size(mycmap)));
