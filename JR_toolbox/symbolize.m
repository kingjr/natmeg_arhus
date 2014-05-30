function [symbols_dec dict] = symb(signal,varargin)
% [symbols_dec dictionary] = symb(signal)
% ---------------------------------------------------------------------
%  Transforms a continuous signal into a sequence of ranking symbols.
%
% input:
%   signal:     continuous signal
%   [kernel=3]: size of the symbole
%   [tau=1]:    number of samples separating two consecutive time points
%   [delta=1]:  number of samples separating two consecutive symbols
%   [dim=1]:    dimension in which the symbol should be transform
%   [symbol_value=1]:    dimension in which the symbol should be transform
%   [keep_order=false]:  for dictionary order
% output:
%   symbols_dec: sequence of symbols 
%   dictionary: all possible symbols
%
% ---------------------------------------------------------------------
% Jean-Rémi King, jeanremi.king [at] gmail.com
% ---------------------------------------------------------------------

%% default parameters
if nargin == 1, varargin = {}; end
for ii = 1:2:length(varargin)
    eval([varargin{ii} '=varargin{ii+1};']);
end

if ~exist('kernel', 'var'),         kernel = 3; end
if ~exist('tau', 'var'),            tau = 1; end
if ~exist('delta', 'var'),          delta = 1; end
if ~exist('dim', 'var'),            dim = 1; end
if ~exist('keep_order', 'var'),     keep_order = false; end
if ~exist('symbol_value', 'var'),   symbol_value = 'continuous'; end

%% arrange dimensions
signal = permute(signal,[dim setdiff(1:length(size(signal)),dim)]);

%% make sure first dimension is used
if length(size(signal)) == 1 && size(signal,1) == 1, signal = signal'; end
%% concatenate across last+1 dimension
sz = size(signal);
if sz(end) == 1 && length(sz)>1, sz(end) = [];end
add_dim = length(sz)+1;
add_dims = repmat(',:',1,add_dim-2);

%% transform into symbol space
m = max([tau delta]);
builder = 'srank=cat(add_dim' ;
for k = 1:kernel
    builder= cat(2,builder,[', signal(' num2str((k-1)*m+1) ':delta:(end-' num2str((kernel-k)*m) ')' add_dims ')']);
end
builder = [builder ');'];
eval(builder);
clear signal;
[~, order] = sort(srank,add_dim);
if keep_order
    [~, order] = sort(order,add_dim);
end
clear srank

%% build dictionary
dict = sortrows(perms(1:kernel),kernel:-1:1);
dict_base = sbase(permute(dict,[1 add_dim 2:(add_dim-1)]),kernel);

%% change base
symbols = sbase(order,kernel);
clear order

%% transform into 1:kernel values
symbols_dec = symbols;
if strcmp(symbol_value, 'continuous')
    for symb = 1:length(dict)
        symbols_dec(symbols_dec==dict_base(symb)) = symb;
    end
end
clear symbols;

%% arrange dimensions back
if length(sz)>1
    symbols_dec = permute(symbols_dec,[2:dim 1 (dim+1):length(sz)]);
end
return

function symbols=sbase(index,kernel)
%% computing symbols: fast computation
sz = size(index);
index = reshape(index,[],kernel);
symbols = zeros(length(index),1);
for ii = 1:kernel
    symbols = symbols+(index(:,ii)-1) * (kernel^(ii-1));
end
if length(sz)>2
    symbols = reshape(symbols,sz(1:(end-1)));
end