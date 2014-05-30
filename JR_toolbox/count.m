function out = count(m, varargin)
% [element_number] = count(matrix, varargin)
% options:
%       - e: vector of values to be counted. By default, takes all
%           possible values.
%       - out: 'en', 'e', 'n'
if nargin == 1, varargin = {}; end
for ii = 1:2:length(varargin);
    eval([varargin{ii} '=varargin{ii+1};']);
end
%-- user didn't specifed variable to search: use all unique values
if ~exist('e', 'var')
    e = unique(m(:));
    if ~isempty(find(isnan(e(:)), 1))
        e = cat(1,e(~isnan(e)), NaN);
    end
    e = reshape(e,[],1);
end
%-- what to output
if ~exist('out', 'var')
    out = 'en';
end
%-- count
n = NaN(length(e),1);
for ii = 1:length(e)
    if isnan(e(ii))
        n(ii) = sum(isnan(m(:)));
    else
        n(ii) = sum(m(:)==e(ii));
    end
end
%-- output
switch out
    case 'en', out = cat(2,e,n);
    case 'e',  out = e;
    case 'n',  out = n;
end