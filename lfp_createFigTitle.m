function result = lfp_createFigTitle(hA, figtype, trials, win, ...
    extraline, clickstr, varargin)
%result = lfp_createFigTitle(hA, figtype, trials, win, extraline, clickstr)
%
%INPUTS
% figtype: a string that identifies what type of analysis was done.
% trials: a list of trialnums, or a string describing the trial set.  If it
%   is a string, then current values of lfp_lib globals are not consulted
%   in generating the fig title.
% win: the time window relative to lfp_AlignmentRef that was analyzed, or
%   [] if there was no defined window.
% extraline: a string appended as an extra line to the title, or '' if no
%   extra line is needed
% clickstr: a string containing text to display in the command window when
%   the title is clicked, or '' if no such text is needed
% hA: the handle to the 'axes' which is to be titled, or [] to return a
%   figtitle structure for later use with 'figtitle' option
%  
%OUTPUT
% result: handle to the text object which is the 'Title' property of <hA>.
%   The first line contains the identification of the session, trials, and
%   <figtype>.
%   The second line contains lfp_AlignmentRef, and <win>.  If
%   <hA> is [], the <result> is a 'figtitle' structure.
%
%OPTIONS
% 'alignment', alignment - overrides the current value of lfp_AlignmentRef
%   in the title text.  <alignment> is formatted as text by dg_thing2str.
% 'append', str - appends the string <str> to the second line, after <win>.
% 'figtitle', figtitle - <figtitle> must be a value returned by a previous
%   call to lfp_createFigTitle with hA=[].  The arguments figtype, trials,
%   win, extraline, clickstr are all ignored in this case, and no lfp_lib
%   globals are referenced.
% 'filenames', filenames - any number of file names go right after the
%   sessionID(s). Can be a string or cell string array.
% 'oneline' - displays only <figtype> plus the message 'click for details',
%   and prepends the entire text of the full title to <clickstr>.  Only has
%   an effect when hA is not empty.
% 'sessionstr', sessionstr - string that overrides the current value of
%   lfp_SessionNames as the source of <sessionstr>, which is the first item
%   on the second line of the title.  By default, all the session names in
%   lfp_SessionNames are strung together with delimiting spaces in between.

%$Rev: 394 $
%$Date: 2019-02-14 18:06:15 -0500 (Thu, 14 Feb 2019) $
%$Author: dgibson $

global lfp_TrialStyle lfp_SessionNames lfp_AlignmentRef lfp_EventNames ...
    lfp_AlignStyle

alignment = [];
appendstr = '';
figtitle = [];
filenames = [];
onelineflag = false;
sessionstr = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'alignment'
            argnum = argnum + 1;
            alignment = varargin{argnum};
        case 'append'
            argnum = argnum + 1;
            appendstr = varargin{argnum};
        case 'figtitle'
            argnum = argnum + 1;
            figtitle = varargin{argnum};
        case 'filenames'
            argnum = argnum + 1;
            filenames = varargin{argnum};
        case 'oneline'
            onelineflag = true;
        case 'sessionstr'
            argnum = argnum + 1;
            sessionstr = varargin{argnum};
        otherwise
            error('lfp_createFigTitle:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

if isempty(alignment) && isnumeric(trials)
    if isequal(lfp_AlignStyle, 'ID')
        alignment = lfp_AlignmentRef;
    else
        if length(lfp_AlignmentRef) == 1
            alignment = lfp_EventNames{lfp_AlignmentRef};
        else
            alignment = '[';
            for k = 1:length(lfp_AlignmentRef)
                alignment = sprintf( '%s%s', alignment, ...
                    lfp_EventNames{lfp_AlignmentRef(k)} );
                if k == length(lfp_AlignmentRef)
                    alignment(end+1) = ']'; %#ok<AGROW>
                else
                    alignment(end+1) = ' '; %#ok<AGROW>
                end
            end
        end
    end
end

if isempty(figtitle)
    if isnumeric(trials) && length(trials) == 1
        firstline = sprintf('trial %s (%d)', ...
            lfp_getTrialID(trials), trials);
    else
        if isnumeric(trials)
            firstline = sprintf('%s', ...
                lfp_getTrialsLabel(trials, lfp_TrialStyle));
        else
            firstline = trials;
        end
    end
    if isempty(win)
        winstr = '[]';
    elseif isnumeric(win) && length(win) == 2
        if isequal(fix(win),win)
            winstr = mat2str(win);
        else
            winstr = sprintf('[%.3f %.3f]', win(1), win(2));
        end
    else
        winstr = '?';
    end
    if ~ischar(sessionstr)
        sessionstr = lfp_SessionNames{1};
        for k = 2:length(lfp_SessionNames)
            sessionstr = [ sessionstr ' ' lfp_SessionNames{k} ]; %#ok<AGROW>
        end
    end
    if ~isempty(filenames)
        if ischar(filenames)
            sessionstr = [ sessionstr ': ' filenames ];
        else
            if ~iscell(filenames)
                error('lfp_createFigTitle:filenames', ...
                    'Option value <filenames> must be either a string or a cell array.');
            end
            sessionstr = [ sessionstr ':' ];
            for k = 1:length(filenames)
                sessionstr = [ sessionstr ' ' filenames{k} ]; %#ok<AGROW>
            end
        end
    end
    titlestr = sprintf('%s %s\n%s align=%s win=%s', ...
        firstline, figtype, sessionstr, dg_thing2str(alignment), winstr);
    if ~isempty(appendstr)
        titlestr = sprintf('%s %s', titlestr, appendstr);
    end
    if ~isempty(extraline)
        titlestr = sprintf('%s\n%s', titlestr, extraline);
    end
else
    if isempty(hA)
        error('lfp_createFigTitle:hA', ...
            '<hA> cannot be empty when using ''figtitle''');
    end
    titlestr = figtitle.titlestr;
    figtype = figtitle.figtype;
    clickstr = figtitle.clickstr;
end    

if isempty(hA)
    result.titlestr = titlestr;
    result.figtype = figtype;
    result.clickstr = clickstr;
else
    result = get(hA, 'Title');
    set(result, 'Interpreter', 'none');
    if ~isempty(clickstr)
        set(result, 'HandleVisibility', 'on');
    end
    if onelineflag
        set(result, 'String', sprintf('%s  click for details', figtype));
        set(result, 'ButtonDownFcn', @dg_showUserData, ...
                'UserData', titlestr);
    else
        set(result, 'String', titlestr);
        if ~isempty(clickstr)
            set( result, 'ButtonDownFcn', @dg_showUserData, ...
                'UserData', clickstr );
        end
    end
end

