function folds = stratified_folding(labels,n_fold,varargin)
% folds = stratified_folding(labels,n_fold,[randomize_fold])
% ------------------------------------------------------------------------- 
% inputs:   labels: vector of n x 1 integers
%           n_fold: number of fold
%           [randomize_fold] : default = true; avoids that same folds gets
%           least number of trials of each condition in case of imbalanced
%           datasets.
% output:   folds:  cell with index of each trial in each fold
% 
% example:  labels = randi(4,[500,1]);
%           folds = stratified_folding(labels, 10);
%           for f = 1:length(folds), 
%               disp(count(labels(folds{f})));
%           end
% -------------------------------------------------------------------------
% Jean-Rémi King, jeanremi.king+matlab@gmail.com, 22/12/2013
% -------------------------------------------------------------------------

%% options
if nargin >2
    for ii = 1:length(varargin)
        eval([varargin{2+ii} '= varargin{3+ii};']);
    end
end
if ~exist('randomize_fold', 'var'), randomize_fold = true; end

%% checkup   
n_label = count(labels); % count number of instances of each number found

if min(n_label(:,2))<n_fold,
    error('too many folds for the least numbered category');
end

%% distribute across folds
used = false(length(labels),1);
for fold = n_fold:-1:1
    folds{fold} = [];
    for condition = 1:size(n_label,1)
        % find trial of right condition that are not already in a fold
        index = find(labels==n_label(condition,1) & ~used);
        % Number of trial of each condition in each fold
        index = index(1:floor(n_label(condition,2)/n_fold));
        % save in fold
        folds{fold} = cat(1,folds{fold}, index);
        % removed from remaining set
        used(index) = true;
    end
end

%% redistributre remaining trials in case indivisible numbers per fold
for condition = 1:size(n_label,1)
    if randomize_fold
        folds_order = randperm(n_fold);
    else
        folds_order = 1:n_fold;
    end
    for fold = folds_order
        % find trial of right condition that are not already in a fold
        index = find(labels==n_label(condition,1) & ~used);
        % Number of trial of each condition in each fold
        if length(index)>=1
            index = index(1);
            % save in fold
            folds{fold} = cat(1,folds{fold}, index);
            % removed from remaining set
            used(index) = true;
        end
    end
end

return

function out = count(m, varargin)
% [element_number] = count(matrix, varargin)
% options:
%       - e: vector of values to be counted. By default, takes all
%           possible values.
%       - out: 'en', 'e', 'n'

%% options
if nargin == 1, varargin = {}; end
for ii = 1:2:length(varargin);
    eval([varargin{ii} '=varargin{ii+1};']);
end
%% user didn't specifed variable to search: use all unique values
if ~exist('e', 'var')
    e = unique(m(:));
    if ~isempty(find(isnan(e(:)), 1))
        e = cat(1,e(~isnan(e)), NaN);
    end
    e = reshape(e,[],1);
end
%% what to output
if ~exist('out', 'var')
    out = 'en';
end
%% count
n = NaN(length(e),1);
for ii = 1:length(e)
    if isnan(e(ii))
        n(ii) = sum(isnan(m(:)));
    else
        n(ii) = sum(m(:)==e(ii));
    end
end
%% output
switch out
    case 'en', out = cat(2,e,n);
    case 'e',  out = e;
    case 'n',  out = n;
end