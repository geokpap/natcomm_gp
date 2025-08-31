function test_dg_aggregate

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

%%%%%%%%%%%%%%%%%%%%%%%%%%%% test cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1:
nums{1}= [
    1 2 3 4
    2 2 3 4 ];
cells{1} = {
    [] 2 [] 4
    1  2 [] [] };
cols{1} = 4;
numsout{1}= [
    1 2 3 4
    2 2 3 4 ];
cellsout{1} = {
    [] 2 [] 4
    1  2 [] [] };
%2:
nums{2}= [
    1 2 3 4
    2 2 3 4 ];
cells{2} = {
    [] 2 [] 4
    1  2 [] [] };
cols{2} = 1;
numsout{2}= [
    NaN 2 3 4 ];
cellsout{2} = {
    1  [2;2] [] 4 };
%3:
nums{3}= [
    1 2 3 4
    1 2 4 4
    2 2 3 4 ];
cells{3} = {
    [] 2 [] 4
    7  8  9 []
    1  2 [] [] };
cols{3} = 1;
numsout{3}= [
    1   2 4 4
    NaN 2 3 4 ];
cellsout{3} = {
    7  8  9 []
    1  [2;2] [] 4 };
%4:
r=4;
nums{r}= [
    1 2 3 4
    1 2 4 4
    2 2 3 4 ];
cells{r} = {
    [] 2 [] 4
    7  8  9 []
    1  2 [] [] };
cols{r} = [1 4];
numsout{r}= [
    1   2 4 4
    NaN 2 3 NaN ];
cellsout{r} = {
    7  8  9 []
    1  [2;2] [] 4 };
%5:
r=5;
nums{r}= [
    1 2 3 4
    1 2 4 4
    2 2 3 4 ];
cells{r} = {
    [] 2 [] 4
    7  8  9 []
    1  2 [] [] };
cols{r} = [4];
numsout{r}= [
    1 2 3 4
    1 2 4 4
    2 2 3 4 ];
cellsout{r} = {
    [] 2 [] 4
    7  8  9 []
    1  2 [] [] };
%6:
r=6;
nums{r}= [
    1 2 3 4
    1 2 4 5
    1 2 4 6
    2 2 3 7 ];
cells{r} = {
    []                  {'ack'} {[]}    {'phffft'}
    {'foo';'bar';'baz'} []      {'boom'}     []
    {'one';'two'}       []      []           []
    {'last'}            {'end'} []       {'g';'mung'}
    };
cols{r} = [1 4];
numsout{r}= [
    NaN 2 3 NaN
    NaN 2 4 NaN
    ];
cellsout{r} = {
    {'last'}            {'ack';'end'} {[]} {'phffft';'g';'mung'}
    {'foo';'bar';'baz';'one';'two'} []         {'boom'}    []
};
%7:
r=7;
nums{r}= [
    1 2 3 4
    1 2 4 5
    1 2 4 6
    2 2 3 7 ];
cells{r} = {
    []                  {'ack'} {[]}    {'phffft'}
    {'foo';'bar';'baz'} []      {'boom'}     []
    {'one';'two'}       []      []           []
    {'last'}            {'end'} []       {'g';'mung'}
    };
cols{r} = [1:4];
numsout{r}= [
    NaN NaN NaN NaN
    ];
cellsout{r} = {
    {'foo';'bar';'baz';'one';'two';'last'} {'ack';'end'} {[];'boom'} {'phffft';'g';'mung'}
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%% end test cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%


for r = 1: length(nums)
    [a b] = dg_aggregate(nums{r}, cells{r}, cols{r});
    if ~isequalwithequalnans(a, numsout{r}) ...
            ||  ~isequalwithequalnans(b, cellsout{r})
        disp([sprintf('Failed test %d: ', r) dg_thing2str(nums{r}) ', ' ...
            dg_thing2str(cells{r}) ', ' dg_thing2str(cols{r})]);
        return
    end
end
disp(sprintf('All %d tests completed successfully',r));
