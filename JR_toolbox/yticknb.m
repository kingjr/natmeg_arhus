function yticknb(varargin)
switch length(varargin)
    case 3,
        handle = varargin{1};
        n = varargin{2};
        keep = varargin{3};
    case 2
        handle = varargin{1};
        n = varargin{2};
        keep = false;
    case 1
        handle = gca;
        n = varargin{1};
        keep = false;
end
if keep
    ytick = get(handle, 'ytick');
else
    ytick = ylim;
end
if n == 0
    set(handle,'ytick', []);
else
    set(handle,'ytick', linspace(min(ytick), max(ytick),n));
end
