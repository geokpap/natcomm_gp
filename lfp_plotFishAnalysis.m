function lfp_plotFishAnalysis(EPdata, extrema, filenames, opts, blockdirs)
%INPUTS
% EPdata: a cell array of 3D arrays where each plane, given by
%   EPdata{refidx}(:,:,chanidx,colnum), contains the data returned by the
%   call to lfp_disp for the channel number <chanidx> of channel group
%   <refidx> to be plotted in column number <colnum>.
% extrema: a cell array containing one cell for each channel group (as
%   defined by the value of <tackle.localavgrefs> submitted to
%   lfp_fishAnalysis), where each cell contains an array in channels X
%   blocks format of a structure containing fields 'sigminTS' and
%   'sigmaxTS' as calculated by lfp_fishAnalysis>analyzeOneChannel.
% filenames: a cell array of cell arrays of strings, where each top-level
%   cell contains the names of all the non-empty files in one channel
%   group.

%$Rev: 266 $
%$Date: 2012-03-07 17:26:42 -0500 (Wed, 07 Mar 2012) $
%$Author: dgibson $

% persistent hF probably no longer needed

% Plot one channel (i.e. one row of plotz)
copydata = [];

for refidx = 1:length(EPdata)
    for chanidx = 1:size(EPdata{refidx}, 3)
        ncols = size(EPdata{refidx}, 4);
        if mod(chanidx-1, 4) == 0
            hF = figure;
        end
        rownum = mod(chanidx-1, 4) + 1;
        for colnum = 1:ncols
            blocknum =  str2double(regexprep(blockdirs(colnum).name, ...
                '^block(\d+)$', '$1'));
            filename = filenames{refidx}{chanidx};
            [p, namestem] = fileparts(filename); %#ok<ASGLU>
            hA = subplot(4, ncols, ...
                (rownum-1)*ncols + colnum, 'Parent', hF);
            dg_plotShadeCL(hA, ...
                [EPdata{refidx}(:,1,chanidx,colnum) ...
                EPdata{refidx}(:,2,chanidx,colnum) - ...
                opts.numSEMs*EPdata{refidx}(:,3,chanidx,colnum)/2 ...
                EPdata{refidx}(:,2,chanidx,colnum) + ...
                opts.numSEMs*EPdata{refidx}(:,3,chanidx,colnum)/2 ...
                EPdata{refidx}(:,2,chanidx,colnum) ], 'Color', [0 0 0]);
            if opts.sumflag
                if colnum == 1
                    sumdata = EPdata{refidx}(:,1:2,chanidx,colnum);
                elseif colnum < ncols
                    sumdata(:,2) = sumdata(:,2) + ...
                        EPdata{refidx}(:,2,chanidx,colnum);
                else
                    % last column
                    set(hA, 'NextPlot', 'add');
                    plot(sumdata(:,1), sumdata(:,2));
                end
            end
            if opts.sumcol
                cols2sum = setdiff(1:ncols, opts.sumcol);
                if colnum == cols2sum(1)
                    sumdata = EPdata{refidx}(:,1:2,chanidx,colnum);
                elseif colnum ~= opts.sumcol
                    sumdata(:,2) = sumdata(:,2) + ...
                        EPdata{refidx}(:,2,chanidx,colnum);
                end
            end
            if opts.copycol(1) && colnum == opts.copycol(1)
                copydata = EPdata{refidx}(:,1:2,chanidx,colnum);
            end
            ylabel(hA, namestem, 'Interpreter', 'none');
            title(hA, sprintf('block %d', blocknum), 'Interpreter', 'none');
            plot(get(hA,'XLim'), [0 0], ...
                'Color', [1 1 1]*0.5);
            
            % plot extrema (if there are any)
            if ~isempty(extrema{refidx}(chanidx,colnum).sigmaxTS) || ...
                    ~isempty(extrema{refidx}(chanidx,colnum).sigminTS)
                current_ylim = get(hA, 'YLim');
                for k=1:length(extrema{refidx}(chanidx,colnum).sigmaxTS);
                    plot([1 1]*extrema{refidx}(chanidx,colnum).sigmaxTS(k), ...
                        current_ylim, 'Color', [.6 0 .6]);
                end
                for k=1:length(extrema{refidx}(chanidx,colnum).sigminTS);
                    plot([1 1]*extrema{refidx}(chanidx,colnum).sigminTS(k), ...
                        current_ylim, 'Color', [0 .7 0]);
                end
                set(hA, 'YLim', current_ylim);  % shouldn't be needed, BUT...
            end
        end
        % Add comparison trace(s)
        if opts.sumcol
            hA = subplot(4, ncols, ...
                (rownum-1)*ncols + opts.sumcol, 'Parent', hF);
            set(hA, 'NextPlot', 'add');
            plot(sumdata(:,1), sumdata(:,2));
        end
        if opts.copycol(1)
            if isempty(copydata)
                warning('lfp_fishAnalysis:copydata', ...
                    'There were no data in column %d, can''t copy', ...
                    opts.copycol(1));
            else
                hA = subplot(4, ncols, ...
                    (rownum-1)*ncols + opts.copycol(2), 'Parent', hF);
                set(hA, 'NextPlot', 'add');
                plot(copydata(:,1), copydata(:,2));
            end
        end
        % Set common vertical scale for row
        ylims = NaN(2, ncols);
        for k = 1:ncols
            hA = subplot(4, ncols, (rownum-1)*ncols + k);
            ylims(:,k) = get(hA, 'YLim');
        end
        commonylim = [min(ylims(1,:)) max(ylims(2,:))];
        for k = 1:ncols
            hA = subplot(4, ncols, (rownum-1)*ncols + k);
            set(hA, 'YLim', commonylim);
        end
    end
end

