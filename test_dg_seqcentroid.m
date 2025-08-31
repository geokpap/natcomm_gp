function test_dg_seqCentroid

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

testseqs = {[1 1 1 1]
[1 2 1 3]
[1 1 1 3]
[1 2 2 1]
[1 2 3 3]
[2 2 3 3]
[4 4 4 4]
[4 5 6 7]
[5 5 6 7]
[5 4 4 4]
[4 5 4 5]
[8 9]
[9 10]
[8 8 9]
[8 9 8 9]
[1 2 6 7]
[1 2 8 9]
[4 5 8 9]};
match = zeros(11);
match(find(eye(11)~=1)) = 1;
indel(1:11) = 1;

% col 1=seqs, 2=minidx, 3=mindist
inout = {
    testseqs(1:6) 2 8
    testseqs(7:11) 1 10
    testseqs(12:15) 1 5
    testseqs([1:11 16]) 2 30
    testseqs([7:15 18]) 10 24
    testseqs([1:6 12:15 17]) 11 25
    testseqs(1:17) 17 47
    testseqs([1:15 18]) 2 48
    testseqs 17 49
};

for r = 1: size(inout,1)
    [a b] = dg_seqCentroid(inout{r,1}, match, indel);
    if ~isequal(a, inout{r,2}) ...
            ||  ~isequal(b, inout{r,3})
        disp(['Failed r ' num2str(r)]);
        return
    end
end
disp('Test completed successfully');
