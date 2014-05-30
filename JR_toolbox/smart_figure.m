function [h] = smart_figure(name,show_now)
%   SMART_FIGURE
%
%       [H] = SMART_FIGURE(NAME,SHOW_NOW)
% 
%       Returns the handler for the figure with a given name or creates
%       it if it does not exist
%
%   Created by Alexandre Gramfort on 2008-06-25.
%   Copyright (c) 2007-2009 Alexandre Gramfort. All rights reserved.

% % $Id: smart_figure.m 171 2009-10-22 13:23:06Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-10-22 15:23:06 +0200 (Thu, 22 Oct 2009) $
% $Revision: 171 $

me = 'SMART_FIGURE';

if nargin == 1
    show_now = true;
end

VISIBLE = true;
% VISIBLE = false;

hw = get(0,'children');
hw = sort(hw);

for i = 1:length(hw),
  s = get(hw(i),'Name');
  if(strcmp(deblank(s),deblank(name))),
    h = hw(i);
    if show_now
        % figure(h);
        set(0,'CurrentFigure',h);
    end
    return;
  end
end

% we exited out without a match, make the window
if VISIBLE
    h = figure;
else
    h = figure('Visible','off');
end
set(h,'Name',name);

end %  function

