function  [handle_line handle_fill] = plot_eb(X,Y,varargin)
% [handle_line handle_fill] = plot_eb(X,Y,[options])
%--------------------------------------------------------------------------
% plot error bar , deals with with NaN and discontinuous segments 
%--------------------------------------------------------------------------
% Inputs: 
%   X: vector x of size m
%   Y: n x m matrix where n is the different instances
%   ['type']: 'sem', 'std', 'prcentile'. Default: 'sem'
%   ['color']: default: [.2 .3 .9]
%   [plotting argument ('FaceColor', 'EdgeAlpha', 'linestyle' etc) can be
%   passed through.
% Outputs:
%   handles
%--------------------------------------------------------------------------
% (c) JeanRémi King, jeanremi.king+matlab [àt] gmail [d0t] com
%--------------------------------------------------------------------------

%% options
if nargin < 3
    varargin = {};
else
    if length(varargin)==1
        varargin = {'color', varargin{:}};  % compatibility wrapper
    end
end
for ii = 1:2:length(varargin)
    eval([varargin{ii} '=varargin{ii+1};']);
end
if ~exist('color', 'var'), color   = [.2 .3 .9]; end
if ~exist('type', 'var'), type = 'sem'; end

%% check errors
if size(X,2)==1, X = X';end
if size(X,2)~=size(Y,2), error('X and Y must have identical dimension 2'); end
if length(size(Y))>2, error('Y must be an n*m matrix'); end

%% compute nanmean
m=nanmean(Y);

s = NaN;
switch type
    case 'sem'
        for x = 1:length(m)
            s(x)=nanstd(Y(:,x)) ./ sqrt(size(Y,1) - sum(isnan(Y(:,x))));
        end
        area = [m-s;m+s];
    case 'std'
        for x = 1:length(m)
            s(x)=nanstd(Y(:,x));
        end
        area = [m-s;m+s];
    case 'prctile'
        for x = length(m):-1:1
            s1(x)=prctile(Y(:,x),10);
            s2(x)=prctile(Y(:,x),90);
        end
        area = [s1;s2];
end

%% find contiguous segments
stops = [0 find(mod(isnan(nanmean(Y(:,2:end)))+~isnan(nanmean(Y(:,1:(end-1)))),2)==0) size(Y,2)];

%% passing variable in plotting function function
cleanvar = @(rm) reshape(cat(1,...
    varargin(...
    setdiff(1:2:length(varargin),...
    (find(ismember(lower(varargin(1:2:end)),lower(rm)))-1)*2+1)),...
    varargin(...
    setdiff(2:2:length(varargin),...
    find(ismember(lower(varargin(1:2:end)),lower(rm)))*2))),1,[]);

%% plot line
vars = cat(2,{'linewidth',2,'Color',color},cleanvar({ 'type', 'rmnan', 'color', 'FaceColor', 'FaceAlpha', 'EdgeColor', 'EdgeAlpha'}));
for ii = (length(stops)-1):-1:1
    toi = (stops(ii)+1):stops(ii+1);
    handle_line{ii}=plot(X(toi),m(toi),vars{:});
    if ~exist('h','var'), 
        h = ishold;
        hold on;
    end
end

%% plot error bar
for ii = (length(stops)-1):-1:1
    toi = (stops(ii)+1):stops(ii+1);
    vars = cat(2,{'edgecolor',color, 'FaceAlpha', .3, 'EdgeAlpha', .5},cleanvar({ 'type', 'rmnan', 'color'}));
    handle_fill{ii}=fill([X(toi) X(toi(end:-1:1))],[area(1,toi) area(2,toi(end:-1:1))],color);
    set(handle_fill{ii},vars{:});
end

%% put hold back
if ~h, hold off; end