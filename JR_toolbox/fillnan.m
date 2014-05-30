function handle = fillnan(X,Y,C,varargin)
% handle = fillnan(X,Y,C,varargin)
% same as fill but deals with NaN as multiple objects
%--------------------------------------------------------------------------
% (c) JeanRémi King, jeanremi.king+matlab [àt] gmail [d0t] com
%--------------------------------------------------------------------------

% options
if nargin < 3
    varargin = {};
end

% find contiguous segments
stops = [1 find(mod(isnan(mean(Y(:,2:end)))+~isnan(mean(Y(:,1:(end-1)))),2)==0) length(Y)];
% plot continuous segments
for ii = (length(stops)-1):-1:1
    toi = (stops(ii)+1):stops(ii+1);
    handle{ii} =fill([X(toi) X(toi(end:-1:1))],[Y(1,toi) Y(2,toi(end:-1:1))],C,varargin{:});
    if ~exist('h','var'), 
        h = ishold;
        hold on;
    end
end
% put hold back
if ~h, hold off; end