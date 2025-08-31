function dg_plothht(imf, Ts)
%dg_plothht(imf, Ts); for Empirical Mode Decomposition (Hilbert-Huang
%   Transform)
% plot one IMF together with various aspects of its Hilbert analytic
% function Ts: time, s step in seconds

%$Rev: 48 $
%$Date: 2010-03-18 13:54:03 -0400 (Thu, 18 Mar 2010) $
%$Author: dgibson $

h = hilbert(imf);
N = length(imf);
c = linspace(0,(N-1)*Ts,N);
th = unwrap(angle(h));
freq = [NaN diff(th)/Ts/(2*pi)];
mag = abs(h);
figure;
subplot(4,1,1);
plot(c, imf);
grid on
title('component');
xlabel('time, s');
subplot(4,1,2);
plot(c, th/pi);
grid on
title('phase');
xlabel('time, s');
ylabel('*pi');
subplot(4,1,3);
plot(c, freq);
grid on
title('frequency');
xlabel('time, s');
ylabel('Hz');
subplot(4,1,4);
plot(c, mag);
grid on
title('magnitude');
xlabel('time, s');
