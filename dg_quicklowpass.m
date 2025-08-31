function result = dg_quicklowpass(Fc, Fs, samples)
% Performs a lowpass filtering operation on the data in <samples> using a
% FIR approximation to the impulse response of the 8th-order zero phase
% Butterworth filter implemented in lfp_bandpass2 (4th-order filtering
% followed by 4th-order reverse-time filtering).  Normalized for unity gain
% at DC.
%INPUTS
% Fc: nominal cutoff frequency (-6 dB point)
% Fs: sample rate
% samples: waveform to filter.
%OUTPUTS
% result: array of the same size as <samples>, containing the result of
%   filtering <samples>.  
%NOTES
%   The waveform of <result> tapers toward zero at the beginning and end due
% to the fact that the computation is done using the Matlab 'conv'
% function.
%   The finite impulse response is stored in the data file
% dg_lowpassFIR.mat, and is the center portion of the impulse response
% produced by lfp_bandpass2 for a cutoff frequency of 80 Hz and a sample
% rate of 40000.6888968566 Hz.  At 100 Hz, the attenuation was 16.85 dB,
% and the rolloff was 48 dB/octave thereafter.

%$Rev: 258 $
%$Date: 2017-06-21 17:36:42 -0400 (Wed, 21 Jun 2017) $
%$Author: dgibson $

% <scalefac> is the number by which the X axis of the FIR must be
% multiplied to achieve the specified <Fc>.
scalefac = (80 / Fc) * (Fs / 40000.6888968566);
s = load('dg_lowpassFIR.mat');
numpts = round(scalefac * length(s.myFIR));
FIR = interp1( 1:length(s.myFIR), s.myFIR, ...
    linspace(1, length(s.myFIR), numpts) );
result = conv(samples(:), FIR, 'same') / sum(FIR);
