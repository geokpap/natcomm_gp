function dg_showUserData(hObject, eventdata) %#ok<INUSD>
% This function is for use with graphics objects (or similar) that have
% already-formatted textual information stored in the 'UserData' property.
% It does nothing more than display it on the command window.

%$Rev: 220 $
%$Date: 2015-05-28 12:27:28 -0400 (Thu, 28 May 2015) $
%$Author: dgibson $

disp(get(hObject, 'UserData'));
