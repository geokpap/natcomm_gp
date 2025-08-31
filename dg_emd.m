function [imf, h, f] = dg_emd(x, varargin)
% [imf, h, f] = dg_emd(x) Empirical Mode Decomposition (Hilbert-Huang
%   Transform)
%INPUTS
% x: the time series to decompose
%OUTPUTS
% imf: an array containing IMFs in components x points form.  The last row
%   is whatever is left over after subtracting the previous rows, and so is
%   less likely than the others to be a true and good IMF.
% h: if 'hilbert' option is invoked, the Hilbert analytic function for each
%   component in <imf> except the last; empty otherwise. 
% f: the frequencies calculated from <h> if 'hilbert' option is invoked;
%   empty otherwise.  The frequency at a given point is calculated from the
%   phase difference between that point and its predecessor, with the first
%   point receiving the value NaN.  Trimmed endpoints are also given the
%   value NaN.  <f> is in units of radians per sample; to convert to Hz,
%   fHz = f / (Ts*2*pi), where <Ts> is the sample period.
%OPTIONS
% 'chatty' - but not so much so as 'verbose'; good for long-running
%   processes.  THe "chatty" message is also sent to the lfp_lib log file
%   if there is one, as determined by the value of lfp_LogFileName.
%   (Strictly speaking, that might make this function a member of lfp_lib
%   rather than dg_lib, but historical inertia prevails.)
% 'hilbert', trimming - the quitting criterion for sifting is that the
%   phase of the Hilbert analytic function be monotonic increasing in an
%   interval defined by <trimming>.  If <trimming> is 'extrema', then the
%   interval is the largest one that contains none of the first maximum,
%   first minimum, last maximum, and last minimum in the raw waveform <x>.
%   If trimming is a number, then this number of samples is trimmed off of
%   each end.  If it is a two-element vector, then the first element is the
%   number of samples trimmed off of the beginning and the second is the
%   number trimmed off the end.  If 'hilbert' is combined with 'minsift',
%   then sifting continues until IMFs that yield monotonic increasing phase
%   and have a constant number of extrema and zero crossings are found <n>
%   times in a row.
% 'hilbertSD', n - equivalent to "'hilbert', 'extrema'" plus the
%   additional criterion that no frequency value be more than <n> standard
%   deviations from the mean.  Can be combined with 'hilbert' to specify
%   some other value for <trimming>.
% 'maxsift', n - similar to S in Huang et al, but provides an arbitrary
%   quitting criterion in case an IMF that meets other criteria cannot be
%   found (default is 5000), i.e. when the total number of sifts done
%   reaches <n>
% 'maxp', maxp - specifies the maximum power that will be permitted
%   at any point in the mean envelope of an IMF, expressed as a fraction of
%   the average power per point in the current residual waveform from which
%   the next IMF is being sought.  Default = .01.
% 'method', method - for dg_emd, ultimately passed through to Matlab
%   function 'interp1'.  Default: 'spline'.
% 'minsift', n - same as S in Huang et al pp 2321-2322: specifies the
%   number of sifts to do after finding an IMF with a constant number of
%   extrema and zero crossings (default is 2)
% 'monotonous', n - sets quitting criterion for outer loop s.t. in the last
%   putative IMF the product of the numbers of peaks and valleys is less
%   than or equal to <n> (default 0, i.e. no more than one full cycle of
%   oscillation). 
% 'verbose' - reports progress of each iteration and summary stats for each
%   IMF

%$Rev: 79 $
%$Date: 2010-08-27 16:20:25 -0400 (Fri, 27 Aug 2010) $
%$Author: dgibson $

% References
%
% Liang H, Bressler SL, Buffalo EA, Desimone R, and Fries P, "Empirical
% mode decomposition of field potentials from macaque V4 in visual spatial
% attention", Biol. Cybern. (2005) 92: 380-392
%
% Huang NE, Wu M-LC, Long SR, Shen SSP, Qu W, Gloersen P, Fan KL, "A
% confidence limit for the empirical mode decomposition and Hilbert
% spectral analysis", Proc Roy Soc Lond A (2003) 459, 2317-2345

% Starting from code in zip file at
% http://www.mathworks.com/matlabcentral/fileexchange/19681-hilbert-huang-transform
%
% 9-Mar-2010 DG
%   Changed quitting criterion in inner loop to match Liang et al.
%   Subject to the interpretation of "has no turning points" as "is
%       monotonic", and "too small" as "less than one part in a million",
%       changed quitting criterion of outer loop to match Liang et al.
% 12-Mar-2010 DG change return value from cell vector to 2-D array in
%   components x points form.
% 16-Mar-2010 DG replaced isimf function call with inline code and changed
%   quitting criterion in inner loop to match Huang et al "second
%   criterion" pp 2321-2322.
% 17-Mar-2010 DG added options 'minsift', 'maxsift', 'hilbert'
% 18-Mar-2010 DG added option 'hilbertSD'; changed function ismonotonic to
%   function ismonotonous

global lfp_LogFileName

chattyflag = false;
hilbertflag = false;
hilbertSD = 0;
minsift = 2;
maxp = .01;
maxsift = 5000;
method = 'spline';
mono_n = 0;
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'hilbert'
            hilbertflag = true;
            argnum = argnum + 1;
            trimming = varargin{argnum};
            if ~( isequal(trimming, 'extrema') || isnumeric(trimming) )
                error('dg_emd:trimming', ...
                    '"%s" is not a legal value for <trimming>', ...
                    dg_thing2str(trimming) );
            end
        case 'chatty'
            chattyflag = true;
        case 'hilbertSD'
            argnum = argnum + 1;
            hilbertSD = varargin{argnum};
        case 'maxp'
            argnum = argnum + 1;
            maxp = varargin{argnum};
        case 'maxsift'
            argnum = argnum + 1;
            maxsift = varargin{argnum};
        case 'method'
            argnum = argnum + 1;
            method = varargin{argnum};
        case 'minsift'
            argnum = argnum + 1;
            minsift = varargin{argnum};
        case 'monotonous'
            argnum = argnum + 1;
            mono_n = varargin{argnum};
        case 'verbose'
            verboseflag = true;
        otherwise
            error('dg_emd:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if hilbertSD
    if ~hilbertflag
        hilbertflag = true;
        trimming = 'extrema';
    end
end

x = reshape(x, 1, []);
imf = [];
h = [];
f = [];
startidx = [];
rawsumabs = sum(abs(x));
s1 = getenvelope(x, method);
s2 = -getenvelope(-x, method);
s3 = (s1+s2)/2;
while ~ismonotonous(x, mono_n) && sum(abs(x)) > 1e-6 * rawsumabs
   x1 = x;
   nsifts = 0;  % number of sifts w/ same numbers of extrema & zero xings
   totalsifts = 0;
   u1prev = NaN; % previous number of 0 xings; NaN => previous iter not IMF
   u2prev = NaN; % previous number of extrema; NaN => previous iter not IMF
   maxpower = maxp * sum(x1.^2)/length(x1);
   while nsifts < minsift && totalsifts < maxsift
       % do the sift
       x1 = x1-s3;
       totalsifts = totalsifts + 1;
       % test the result of the sift
       N  = length(x1);
       u1 = sum(x1(1:N-1).*x1(2:N) < 0);
       peaks = dg_findpks(x1);
       troughs = dg_findpks(-x1);
       if hilbertflag && isempty(startidx)
           % We are finding the first component; find the interval to test
           % for the 'hilbert', 'extrema' option
           if isequal(trimming, 'extrema')
               if isempty(peaks) || isempty(troughs) ...
                       || peaks(1) == 1 || troughs(1) == 1 ...
                       || peaks(end) == length(x) || troughs(end) == length(x)
                   error('Burma!');
               end
               startidx = max(peaks(1), troughs(1)) + 1;
               endidx = min(peaks(end), troughs(end)) - 1;
           elseif isequal(size(trimming), [1 1])
               startidx = trimming + 1;
               endidx = length(x) - trimming - 1;
           else
               startidx = trimming(1);
               endidx = trimming(2);
           end
       end
       s1 = getenvelope(x1, method);
       s2 = -getenvelope(-x1, method);
       s3 = (s1+s2)/2;
       u2 = length(peaks)+length(troughs);
       % mean envelope must not contain any individual points that exceed
       % 1/1000 of the average power in x1, and numbers of zero xings and
       % extrema may not differ by more than one:
       p3 = s3.^2;
       isIMF = all(p3 <= maxpower) && abs(u1-u2) <= 1;
       if isIMF
           % x1 is an IMF; see if this iteration counts towards minsift
           gooditer = false;
           if nsifts == 0 || (u1 == u1prev && u2 == u2prev)
               if hilbertflag
                   hfunc = hilbert(x1);
                   th = unwrap(angle(hfunc));
                   freq = [NaN diff(th)];
                   if hilbertSD
                       fmean = mean(freq(startidx:endidx));
                       fsd = std(freq(startidx:endidx));
                       gooditer = all(freq(startidx:endidx) > 0) ...
                           && all( abs(freq(startidx:endidx) - fmean) ...
                           < hilbertSD * fsd );
                   else
                       gooditer = all(freq(startidx:endidx) > 0);
                   end
               else
                   gooditer = true;
               end
           end
           u1prev = u1;
           u2prev = u2;
       else
           % x1 is not an IMF.
           gooditer = false;
           u1prev = NaN;
           u2prev = NaN;
       end
       if gooditer
           % this iteration counts towards minsift
           nsifts = nsifts + 1;
       else
           % In case it's possible for additional sifts to turn an IMF into
           % a non-IMF, we reset the successive sift counter:
           nsifts = 0;
       end
       if verboseflag
           fprintf( ...
               '%4d %2d meanp3=%.3g maxp3=%.3g gooditer=%d isIMF=%d\n', ...
               totalsifts, nsifts, mean(p3), max(p3), gooditer, isIMF);
           if ~gooditer
               % gooditer is true if:
               % ~(nsifts == 0 || (u1 == u1prev && u2 == u2prev)) 
               % i.e., nsifts ~= 0 && (u1 ~= u1prev || u2 ~= u2prev)
               fprintf('nsifts=%d; u1 %d %d; u2 %d %d\n', ...
                   nsifts, u1, u1prev, u2, u2prev);
           end
           if ~isIMF
               fprintf( ...
                   'all(p3 <= maxpower): %d\nabs(u1-u2) <= 1: %d\n', ...
                   all(p3 <= maxpower), abs(u1-u2) <= 1);
           end
       end
   end
   if nsifts < minsift
       warning('dg_emd:maxsift', ...
           'Terminated search for IMF %d due to maxsift = %d',  ...
           size(imf,1) + 1, maxsift);
   end
   
   if verboseflag || chattyflag
       msg = sprintf('Found IMF %d after %d sifts,  0xings:%d, extrema:%d', ...
           size(imf,1) + 1, totalsifts, u1, u2);
       disp(msg);
       if ~isempty(lfp_LogFileName)
           lfp_log(sprintf('PID %d %s', dg_pid, msg));
       end
   end
   
   imf(end+1,:) = x1;
   if hilbertflag
       h(end+1,:) = hfunc;
       f(end+1,:) = [NaN(1, startidx-1) freq(startidx:endidx) ...
           NaN(1, length(x) - endidx)];
   end
   x = x-x1;
end
imf(end+1,:) = x;

% FUNCTIONS

function u = ismonotonous(x, n)
% More flexible than the old "ismonotonic" (which was misnamed anyway,
% since it returned true for waves having a single extremum). Returns true
% if the product of the numbers of peaks and valleys is less than or equal
% to <n>. Degenerates into "ismonotonic" when n=0.
u1 = length(dg_findpks(x))*length(dg_findpks(-x));
u = u1 <= n;


function s = getenvelope(x, method)

N = length(x);
p = dg_findpks(x);
s = interp1([0 p N+1],[0 x(p) 0], 1:N, method);

