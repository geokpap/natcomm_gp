function test_dg_canonicalSeries

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

inout = {[1 4 5 100] '[ 1 4:5 100 ]'
    [1 4 5 100 110 111 112] '[ 1 4:5 100 110:112 ]'
    [1 4 5 100 110 111] '[ 1 4:5 100 110:111 ]'
    [4 5 100 110 111] '[ 4:5 100 110:111 ]'
    [1 2] '[ 1:2 ]'
    3 '3'
    [1 2 3 4 5 6 7 8 9] '[ 1:9 ]'
    [1 3 5 7 9] '[ 1 3 5 7 9 ]'
};

for r = 1: size(inout,1)
    if ~isequal(dg_canonicalSeries(inout{r,1}), inout{r,2})
        disp(['Failed: dg_canonicalSeries(' mat2str(inout{r,1}) ')']);
        return
    end
end
disp('Test completed successfully');
