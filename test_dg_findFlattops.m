function test_dg_findFlattops

%$Rev: 127 $
%$Date: 2011-08-12 16:54:36 -0400 (Fri, 12 Aug 2011) $
%$Author: dgibson $

r = 1;
samples{r} = [
    1 1 1 2 3 4 3 2 3 4 5 5 4 3 4 3 2 2
    ];
vargs{r} = {};
vmaxidxout{r} = [
     6
    11
    15
    ];
startendsout{r} = [
    11 12
    ];

r = 2;
samples{r} = [
    1 1 1 2 3 4 3 2 3 4 5 5 4 3 4 3 2 2 3 3 3 0 5 9
    ];
vargs{r} = {};
vmaxidxout{r} = [
     6
    11
    15
    19
    ];
startendsout{r} = [
    11 12
    17 18
    19 21
    ];

r = 3;
samples{r} = [
    3 4 3 2 3 4 5 5 4 3 4 3 2 2 3
    ];
vargs{r} = {};
vmaxidxout{r} = [
     2
     7
    11
    ];
startendsout{r} = [
     7 8
    13 14
    ];

r = 4;
samples{r} = [
    4 3 2 2 3 4 5 5 4 3 4 3 2 2 3 2
    ];
vargs{r} = {};
vmaxidxout{r} = [
     7
    11
    15
    ];
startendsout{r} = [
     3 4
     7 8
    13 14
    ];

r = 5;
samples{r} = -[
    4 3 2 2 3 4 5 5 4 3 4 3 2 2 3 2
    ];
vargs{r} = {};
vmaxidxout{r} = [
     3
    10
    13
    ];
startendsout{r} = [
     3 4
     7 8
    13 14
    ];

r = 6;
samples{r} = [
    4 3 2 1 0 0 0 5 4 3 4 3 2 2 1 0
    ];
vargs{r} = {};
vmaxidxout{r} = [
     8
    11
    ];
startendsout{r} = [
     5 7
    13 14
    ];

r = 7;
samples{r} = -[
    4 3 2 1 0 0 0 5 4 3 4 3 2 2 1 0
    ];
vargs{r} = {};
vmaxidxout{r} = [
     5
    10
    ];
startendsout{r} = [
     5 7
     13 14
    ];

r = 8;
samples{r} = [
    4 3 4 3 2 2 1 0 4 3 2 1 0 0 0 5
    ];
vargs{r} = {};
vmaxidxout{r} = [
     3
     9
    ];
startendsout{r} = [
     5 6
    13 15
    ];

r = 9;
samples{r} = -[
    4 3 4 3 2 2 1 0 4 3 2 1 0 0 0 5
    ];
vargs{r} = {};
vmaxidxout{r} = [
     2
     8
    13
    ];
startendsout{r} = [
     5 6
    13 15
    ];

r = 10;
samples{r} = [
    0 0 0 0 0 0 0 0 0 0 2 4 3 4 3 2.5 2 1 0 4 3 2 1 0 0.1 0.2 5
    ];
vargs{r} = {};
vmaxidxout{r} = [
    12
    14
    20
    ];
startendsout{r} = [
    ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%% end test cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1: length(samples)
    [vmaxidx, runidx, endrunidx] = dg_findFlattops(samples{r});
    if ~isequalwithequalnans(vmaxidx, vmaxidxout{r}) || ...
            ~isempty(startendsout{r}) && ...
            ~isequalwithequalnans([runidx, endrunidx], startendsout{r}) ...
            || isempty(startendsout{r}) && ( ...
            ~isempty(runidx) || ~isempty(endrunidx) )
        fprintf('Failed test %d: [vmaxidx, runidx, endrunidx] = dg_findFlattops(', r);
        fprintf('%s', dg_thing2str(samples{r}));
        if ~isempty(vargs{r})
            fprintf(', %s', dg_thing2str(vargs{r}));
        end
        fprintf(')\n');
        return
    end
end
fprintf('All %d tests completed successfully\n',r);
