function [] = progressbar(n,N,w)
%   PROGRESSBAR   Display an ascii progress bar
%       [] = PROGRESSBAR(N,N,W)
%
% displays the progress of n out of N.
% n should start at 1.
% w is the width of the bar (default w=20).
%
%   Created by Alexandre Gramfort on 2008-06-26.
%   Copyright (c) 2007-2009 Alexandre Gramfort. All rights reserved.
%
% $Id: progressbar.m 171 2009-10-22 13:23:06Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-10-22 15:23:06 +0200 (Thu, 22 Oct 2009) $
% $Revision: 171 $

me = 'PROGRESSBAR';

if nargin == 0
    eval(['help ',lower(me)])
    return
elseif nargin == 1
    N = 1;
end

if nargin<3
    w = 20;
end

% progress char
cprog = '.';
cprog1 = '*';
% begining char
cbeg = '[';
% ending char
cend = ']';

p = min( floor(n/N*(w+1)), w);

global pprev;
if isempty(pprev)
    pprev = -1;
end

if not(p==pprev)
    ps = repmat(cprog, [1 w]);
    ps(1:p) = cprog1;
    ps = [cbeg ps cend];
    if n>1
        % clear previous string
        fprintf( repmat('\b', [1 length(ps)]) );
    end
    fprintf(ps);
end
pprev = p;
if n==N
    fprintf('\n');
end