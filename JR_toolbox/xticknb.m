function xticknb(varargin)
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
    xtick = get(handle, 'xtick');
else
    xtick = xlim;
end
if n == 0
    set(handle,'xtick', []);
else
    set(handle,'xtick', linspace(min(xtick), max(xtick),n));
end
