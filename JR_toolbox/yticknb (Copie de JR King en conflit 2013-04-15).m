function yticknb(handle,n,keep_tick)
if nargin < 3, keep_tick = false; end
if nargin < 2, n = handle; handle = gca; end
ytick = get(handle, 'ytick');
if ~keep_tick, ytick = ylim; end
set(handle,'ytick', linspace(min(ytick), max(ytick),n));
