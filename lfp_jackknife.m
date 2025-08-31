function [avg lowCL upCL] = lfp_jackknife(funch, p, trials, varargin)
% Computes result of evaluating <funch> using jackknife method.
% INPUTS
%   funch - function handle to a function that returns a scalar, vector, or
%       2-D array, and that accepts <trials> as its first argument.
%   p - statistical confidence level to use for confidence limits
%   trials - as usual
%   varargin - put any additional arguments that are needed by <funch> here
% OUTPUTS
%   avg - an array of the same shape as the result returned by <funch>.
%   