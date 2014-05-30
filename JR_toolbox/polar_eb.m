function  [handle_line handle_fill] = polar_eb(TH,R,varargin)
% [handle_line handle_fill] = polar_eb(TH,R,[options])
%--------------------------------------------------------------------------
% plot error bar , deals with with NaN and discontinuous segments 
%--------------------------------------------------------------------------
% Inputs: 
%   TH: angle: vector x of size m
%   R: radius: n x m matrix where n is the different instances
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
if size(TH,2)==1, TH = TH';end
if size(TH,2)~=size(R,2), error('X and Y must have identical dimension 2'); end
if length(size(R))>2, error('Y must be an n*m matrix'); end



%% compute nanmean
m=nanmean(R);

s = NaN;
switch type
    case 'sem'
        for x = 1:length(m)
            s(x)=nanstd(R(:,x)) ./ sqrt(size(R,1) - sum(isnan(R(:,x))));
        end
        area = [m-s;m+s];
    case 'std'
        for x = 1:length(m)
            s(x)=nanstd(R(:,x));
        end
        area = [m-s;m+s];
    case 'prctile'
        for x = length(m):-1:1
            s1(x)=prctile(R(:,x),10);
            s2(x)=prctile(R(:,x),90);
        end
        area = [s1;s2];
end

%% find contiguous segments
stops = [0 find(mod(isnan(nanmean(R(:,2:end)))+~isnan(nanmean(R(:,1:(end-1)))),2)==0) size(R,2)];

%% passing variable in plotting function function
cleanvar = @(rm) reshape(cat(1,...
    varargin(...
    setdiff(1:2:length(varargin),...
    (find(ismember(lower(varargin(1:2:end)),lower(rm)))-1)*2+1)),...
    varargin(...
    setdiff(2:2:length(varargin),...
    find(ismember(lower(varargin(1:2:end)),lower(rm)))*2))),1,[]);

%% transform into polar
[X Y] = pol2cart(repmat(TH,size(R,1),1),R);
[X m] = pol2cart(TH,m);
[areaX areaY] = pol2cart(repmat(TH,2,1),area);

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
    vars = cat(2,{'edgecolor','none', 'FaceAlpha', .3, 'EdgeAlpha', .5},cleanvar({ 'type', 'rmnan', 'color'}));
    handle_fill{ii}=fill([areaX(1,toi) areaX(2,toi(end:-1:1))],[areaY(1,toi) areaY(2,toi(end:-1:1))],color);
    set(handle_fill{ii},vars{:});
end

%% put hold back
if ~h, hold off; end