function num = dg_tabread2(filename)
% Extracts numbers from ASCII spreadsheets.
%INPUTS
% filename: absolute or relative path to file to read.
%OUTPUTS
% num: Behaves similarly to num = xlsread(filename), except that <filename>
%   must be the pathname of a text file:  dg_tabread ignores rows or
%   columns of nothing but text; however, if a cell not in a leading row or
%   column is empty or contains text, dg_tabread puts a NaN in its place in
%   the return array, num.
%NOTES
% The presence of even a single numerical field in a row or column causes
% that raw or column to be included in <num>, i.e. NOT to be treated as a
% header.  Runs 5 or 6 times faster than 'dg_tabread'.

%$Rev: 278 $
%$Date: 2021-10-04 18:41:51 -0400 (Mon, 04 Oct 2021) $
%$Author: dgibson $

T = readtable(filename);
A = table2array(T);
num = cellfun(@str2double, A);
cols2del = all(isnan(num), 1);
rows2del = all(isnan(num), 2);
num(rows2del, :) = [];
num(:, cols2del) = [];


