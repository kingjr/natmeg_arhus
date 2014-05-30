function a = getoption(options,key,default_value)
%   GETOPTION   get option parameter
%       [A] = GETOPTION(OPTIONS,KEY,DEFAULT_VALUE)
%
%   a = getoption(options, 'key', val);
%
%   is equivalent to the code:
%
%   if isfield(options, 'key')
%       a = options.key;
%   else
%       a = val;
%   end
%
%   this function is part of the EMBAL toolbox, see COPYING for license
%
%   Copyright (c) 2009 Alexandre Gramfort. All rights reserved.

% $Id: getoption.m 196 2009-11-03 14:59:49Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-11-03 15:59:49 +0100 (Tue, 03 Nov 2009) $
% $Revision: 196 $

if isfield(options, key)
    a = eval(['options.' key ';']);
else
    if nargin < 3
        error(['parameter ',key,' is mandatory'])
    end
    a = default_value;
end

end %  function