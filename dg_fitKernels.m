function [kernels, Sigma, E, CovB, logL] = dg_fitKernels(datary, ...
evtable, evtwindows, coeffs, maxiter)
%   Step 1 of fitting the "Lak...Harris_Carandini_2018_Biorxiv" model to
% the data produced by makeLakHarris.  The model is defined (in
% pseudo-code) as
%   R(t,j) = sum(A(e,j)*K(:,t-TS(e,j))) + C
% where <t> is the time point, <j> is the trial number, <e> is the event
% number, A is the kernel amplitude factor, and K is the array of
% event-aligned kernels.  The sum is taken over all events <e> for a given
% <j> and <t>. C = 0.  
%   On step 1 of each iteration, we take the kernel amplitude factors
% A(e,j) to be fixed a priori values and fit the "response" values in
% <datary> by fitting the values in the event-aligned kernels.  Thus the
% "predictors" on step 1 are the amplitude factors A(e,j).
%INPUTS
% datary: a cell column vector each element of which contains a column
%   vector time series of some measured value (e.g. spike density)
%   representing one trial.
% evtable: a 2D array with the same number of rows as <datary>. Each
%   column contains the sample index into <datary{rownum}> for one
%   event.  Address as evtable(trialnum, evtnum).
% evtwindows: a 2 row array of integers with the same number of columns as
%   <evtable>.  Each column specifies the relative starting index (row 1)
%   and ending index (row 2) of the corresponding event's time window for
%   computing the event-locked response kernel.  The relative index of the
%   alignment event for the window is defined to be 0, and samples before
%   the alignment event have negative relative indices.  Address as
%   evtwindows(start|end, evtnum).
% coeffs: a 2D array with the same number of rows as <datary> and the same
%   number of columns as <evtable>.  <coeffs(j,e)> contains the values of
%   A(e,j) referenced by the model; those values are taken as givens on
%   step 1 of the analysis.
% maxiter: value of the 'maxiter' parameter to 'mvregress'.
%OUTPUTS
% kernels: a cell row vector with the same number of columns as <evtable>.
%   Each element contains a column vector that specifies the time series
%   for the corresponding event-related kernel.  The "time 0" sample (i.e.
%   the alignment sample) for event kernel <e> is
%   kernels{e}(-evtwindows(1,e) + 1).
% C: scalar, the constant term in the model.

%$Rev:  $
%$Date:  $
%$Author: dgibson $

numevents = size(evtable, 2);
numtrials = length(datary);
kernels = cell(0, numevents);
evtmin = min(evtable);  % row vector
evtmax = max(evtable);  % row vector
% evtwinstartmin: scalar sample index into datary{rownum} at which the
% modeled data <Y> begins (see below).
evtwinstartmin = min(evtmin + evtwindows(1,:));
% evtwinendmax: scalar sample index into datary{rownum} at which the
% modeled data <Y> ends (see below).
evtwinendmax = max(evtmax + evtwindows(2,:)); % scalar

% Set up data ('response' or 'Y') values, Y(trialnum, sampidx). (See "Set
% Up Multivariate Regression Problems" in the Matlab help file for further
% explanation of 'X', 'Y', 'd', 'K', and 'n', where 'n' corresponds to our
% <numtrials>.)
%   numpts: number of points modeled; points not modeled in any given trial
% are assigned the value 0 in both the data set <Y> and the model output.
numpts = evtwinendmax - evtwinstartmin + 1; 
Y = zeros(numtrials, numpts);
for trialnum = 1:numtrials
    startidx = min(evtable(trialnum,:) + evtwindows(1,:));
    if any(startidx < 1)
        error('analLakHarris1:startidx', ...
            'Trial %d cannot accommodate start of event window %s.\n', ...
            trialnum, dg_thing2str( ...
            find(evtable(trialnum,:) + evtwindows(1,:) < 1) ));
    end
    endidx = max(evtable(trialnum,:) + evtwindows(2,:));
    if any(endidx > length(datary{trialnum}))
        error('analLakHarris1:endidx', ...
            'Trial %d cannot accommodate end of event window %s.\n', ...
            trialnum, dg_thing2str( ...
            find( evtable(trialnum,:) + evtwindows(2,:) ...
            > length(datary{trialnum}) )));
    end
    Y(trialnum, 1 : endidx - startidx + 1) = ...
        datary{trialnum}(startidx:endidx);
end

% Set up 'predictor' or 'X' design matrix, X{trialnum}(sampidx, coeffidx).
% The number of dimensions 'd' of each observation (trial) is <numpts>.
% The number of linear regression coefficients 'K' is the sum of the
% numbers of points in all the event kernels.  Coefficients
% contain the kernel values concatenated in the same order in which the
% event indices are listed across the columns of <evtable>.
%   kernelidces: the start (first row) and end (second row) indices into
% the vector of linear regression coefficients <beta> that correspond to
% the start and end sample of each event kernel.  Address as
% kernelidces(start|end, evtnum).
kernelidces = zeros(2, numevents);
kernelidces(:, 1) = 1 + [0; evtwindows(2,1) - evtwindows(1,1)];
for evtnum = 2:size(evtwindows,2)
    kernelidces(:, evtnum) = kernelidces(2, evtnum-1) + 1 ...
        + [0; evtwindows(2,evtnum) - evtwindows(1,evtnum)]; 
end
% Because the alignment of each event kernel is generally different from
% trial to trial, we must use the cell vector method for specifying the
% design matrices <X>.  Each value X{trialnum}(sampidx, coeffidx) connects
% a regression coefficient to a <sampidx> in Y(trialnum, sampidx) via the a
% priori amplitude factor for one of the kernels.
X = cell(size(datary));
for trialnum = 1:numtrials
    % The total number of coefficients in the model is kernelidces(2, end).
    X{trialnum} = zeros(numpts, kernelidces(2, end));
    % On trial <trialnum>, we require that <sampidx> = <evtable(trialnum,
    % evtnum)> in Y(trialnum, sampidx) correspond to relative index 0 of
    % the event-locked response kernel for event <evtnum>.  This implies
    % that regression coefficient <kernelidces(1, evtnum)> corresponds to
    % <sampidx> = <evtable(trialnum, evtnum) + evtwindows(1, evtnum)>, and
    % regression coefficient <kernelidces(2, evtnum)> corresponds to
    % <sampidx> = <evtable(trialnum, evtnum) + evtwindows(2, evtnum)>. We
    % therefore need to assign the "predictor" value
    % <coeffs(trialnum,evtnum)> to the main diagonal of the square subarray
    % of X{trialnum} running from 
    %   sampidx = evtable(trialnum, evtnum) + evtwindows(1, evtnum)
    % through
    %   sampidx = evtable(trialnum, evtnum) + evtwindows(2, evtnum)
    % and from
    %   coeffidx = kernelidces(1, evtnum)
    % through
    %   coeffidx = kernelidces(2, evtnum).
    for evtnum = 1:numevents
        coeffidx = kernelidces(1, evtnum);
        for sampidx = evtable(trialnum, evtnum) + evtwindows(1, evtnum) ...
                : evtable(trialnum, evtnum) + evtwindows(2, evtnum)
            X{trialnum}(sampidx, coeffidx) = coeffs(trialnum,evtnum);
            coeffidx = coeffidx + 1;
        end
        if coeffidx ~= kernelidces(2, evtnum) + 1
            error('Oops!');
        end
    end
end

% Run 'mvregress':
[beta, Sigma, E, CovB, logL] = mvregress( X, Y, ...
    'algorithm', 'mvn', ...
    'tolbeta', 1e-3, 'tolobj', 1e-3, ...
    'maxiter', maxiter, ...
    'outputfcn', @analLakHarris1_outfun );
for evtnum = 1:numevents
    C = beta(1);
    kernels{1, evtnum} = ...
        beta(kernelidces(1, evtnum):kernelidces(2, evtnum));
end

end


function rc = analLakHarris1_outfun(beta, thing, status)
% Seems to me that this method of stopping iterations means that you never
% get to see the results of the final iteration; but, OK.
persistent hF hA prevbeta prevfval tolbeta tolobj hA2
rc = false;
if isempty(hF)
    hF = figure;
    hA = subplot(2,1,1);
    hA2 = subplot(2,1,2);
end
if isvalid(hF)
    if thing.iteration == 0
        % First iteration is useless because beta = 0.
        prevbeta = beta;
        prevfval = NaN;
        tolbeta = NaN(0,1);
        tolobj = NaN(0,1);
        set(hA2, 'NextPlot', 'replace');
        plot(hA2, beta, 'x');
    else
        set(hA2, 'NextPlot', 'add');
        if isequal(status, 'done')
            title(hA, 'Done');
        else
            tolbeta(end+1,1) = norm(beta - prevbeta) / ...
                (sqrt(length(beta)) * (1 + norm(beta)));
            tolobj(end+1,1) = abs(thing.fval - prevfval) / ...
                (1 + abs(thing.fval));
            plot(hA, 1:thing.iteration, log10([tolbeta tolobj]), '-+');
            title(hA, 'Iterating');
            legend(hA, {'tolbeta' 'tolobj'});
            prevbeta = beta;
            prevfval = thing.fval;
            plot(hA2, beta, 'x');
        end
        drawnow limitrate
    end
else
    rc = true;
    hF = [];
end
end


