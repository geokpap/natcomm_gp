function test_dg_findScale

%$Rev: 105 $
%$Date: 2011-04-29 16:25:17 -0400 (Fri, 29 Apr 2011) $
%$Author: dgibson $

m = [1 2 5 1];
r = 0;
for e = [-5 0 5]
    for k = 1:3
        for d = [0 -.00001]
            r = r + 1;
            nums{r} = (m(k)+d)*10^e; %#ok<*AGROW>
            mantissa{r} = m(k);
            exp{r} = e;
            opts{r} = {};
        end
        d = .00001;
        r = r + 1;
        opts{r} = {};
        nums{r} = (m(k)+d)*10^e;
        mantissa{r} = m(k+1);
        if k==3
            exp{r} = e+1;
        else
            exp{r} = e;
        end
    end
end

r=28;
nums{r} = [rand(1,3) 1];
mantissa{r} = 1;
exp{r} = 0;
opts{r} = {};

r=29;
nums{r} = 2 * [rand(1,3) 1];
mantissa{r} = 2;
exp{r} = 0;
opts{r} = {};

r=30;
nums{r} = 5 * [rand(1,3) 1];
mantissa{r} = 5;
exp{r} = 0;
opts{r} = {};

r=31;
nums{r} = [1.6 9.7 9.6 4.9 8.0];
opts{r} = {};
mantissa{r} = 1;
exp{r} = 1;

r=32;
nums{r} = [1.6 9.7 9.6 4.9 8.0];
opts{r} = {'levels', 2};
mantissa{r} = [1 5];
exp{r} = [1 0];

r=33;
nums{r} = [1.6 9.7 9.6 4.9 8.0];
opts{r} = {'levels', 3};
mantissa{r} = [1 5 2];
exp{r} = [1 0 0];

r=34;
nums{r} = [1.6 9.7 9.6 4.9 8.0];
opts{r} = {'levels', 4};
mantissa{r} = [1 5 2 1];
exp{r} = [1 0 0 0];

r=35;
nums{r} = [1.6 9.7 9.6 4.9 8.0];
opts{r} = {'levels', 5};
mantissa{r} = [1 5 2 1];
exp{r} = [1 0 0 0];

r=36;
nums{r} = [1.6 9.7 9.6 4.9 8.0];
opts{r} = {'levels', 0};
mantissa{r} = [1 5 2 1];
exp{r} = [1 0 0 0];

r=37;
nums{r} = [1.6 9.7 9.6 4.9 8.0 0.9];
opts{r} = {'levels', 0};
mantissa{r} = [1 5 2 1 5];
exp{r} = [1 0 0 0 -1];

r=38;
nums{r} = [1.6 9.7 9.6 4.9 8.0 0.25];
opts{r} = {'levels', 0};
mantissa{r} = [1 5 2 1 5 2];
exp{r} = [1 0 0 0 -1 -1];

r=39;
nums{r} = [1.6];
opts{r} = {'levels', 0};
mantissa{r} = [2];
exp{r} = [0];



%%%%%%%%%%%%%%%%%%%%%%%%%%%% end test cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%


for r = 1: length(nums)
    [a b] = dg_findScale(nums{r}, opts{r}{:});
    if ~isequalwithequalnans(a, mantissa{r}') ...
            ||  ~isequalwithequalnans(b, exp{r}')
        disp([sprintf('Failed test %d: ', r) dg_thing2str(nums{r})]);
        return
    end
end
fprintf('All %d tests completed successfully\n',r);