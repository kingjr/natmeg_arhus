function xticknb(handle,n,keep_tick)
if nargin < 3, keep_tick = false; end
if nargin < 2, n = handle; handle = gca; end
xtick = get(handle, 'xtick');
if ~keep_tick, xtick = xlim; end
set(handle,'xtick', linspace(min(xtick), max(xtick),n));
