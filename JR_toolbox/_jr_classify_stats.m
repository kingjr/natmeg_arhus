function results = jr_classify_stats(file,varargin)
%  jr_classify_stats(file/results [,varargin])

if nargin == 1
    varargin = {};
end

for ii = 1:2:length(varargin), eval([varargin{ii} '= varargin{ii+1}']); end




%% Load data
if ischar(file), results = load(file);
else results = file; end

%% Parameters
%-- stats
if ~exist('method','var'), method = 'ttest2'; end
%-- data type
if ~exist('datatype','var'), 
    if isfield(results, 'probas'), datatype = 'probas'; 
    elseif isfield(results, 'distance'), datatype = 'distance'; 
    end
end
eval(['data = results.' datatype ';']);
eval(['datag= results.' datatype 'g;']);

g= ~isempty(results.yg);

n_splits    = size(data,1);
n_folds     = size(results.all_folds,2);
classes     = unique(results.y);
n_classes   = length(classes);
gclasses     = unique(results.yg);
n_gclasses   = length(gclasses);
%-- find all 2x2 comparisons
cs          = nchoosek(1:n_classes,2); 

%% each split
[mpg mpgk spg spgk] = deal([]);
for s = size(results.all_folds,1):-1:1
    %% mean and std probabilities
    for c1 = 1:n_classes
        %-- for training categories
        for c2 = 1:n_classes
            mp(s,c1,c2,:,:) = nanmean(data(s,results.y==classes(c2),:,:,c1),2);
            sp(s,c1,c2,:,:) = nanstd(data(s,results.y==classes(c2),:,:,c1),[],2);
        end
        %-- for generalization categories
        for c2 = 1:n_gclasses
            mpg(s,c1,c2,:,:) = nanmean(nanmean(datag(s,results.yg==gclasses(c2),:,:,c1,:),6),2);
            spg(s,c1,c2,:,:) = nanstd(nanmean(datag(s,results.yg==gclasses(c2),:,:,c1,:),6),[],2);
        end
    end
    %% mean probabilities for each fold
    for k = size(results.all_folds,2):-1:1
        test = find(results.all_folds(s,k,:)==0);
        for c1 = 1:n_classes
            %-- for training categories
            for c2 = 1:n_classes
                mpk(s,k,c1,c2,:,:) = nanmean(data(s,intersect(find(results.y==classes(c2)), test),:,:,c1),2);
                spk(s,k,c1,c2,:,:) = nanstd(data(s,intersect(find(results.y==classes(c2)), test),:,:,c1),[],2);
            end
            %-- for generalization categories
            for c2 = 1:n_gclasses
                mpgk(s,k,c1,c2,:,:) = nanmean(datag(s,results.yg==gclasses(c2),:,:,c1,k),2);
                spgk(s,k,c1,c2,:,:) = nanstd(datag(s,results.yg==gclasses(c2),:,:,c1,k),[],2);
            end
        end    
    end
    %% each 2x2 comparison
    for c = size(cs,1):-1:1
        %-- overall
        sel = [find(results.y==classes(cs(c,1)));find(results.y==classes(cs(c,2)))]; % fix 2012 11 29: more than 2 categories
        auc(s,c,:,:,:) = reshape(colAUC(sq(data(s,sel,:,:,cs(c,1))),results.y(sel)),size(data,3),size(data,4));
        %-- each fold
        for k = size(results.all_folds,2):-1:1
            fprintf('*');
            %-- find trials not used in training sets
            test = find(results.all_folds(s,k,:)==0);
            %-- subselect categories for 2x2 comparison
            y1 = intersect(test,find(results.y==classes(cs(c,1))));
            y2 = intersect(test,find(results.y==classes(cs(c,2))));
            %-- compute auc
            try
            auck(s,c,k,:,:) = colAUC(...
                sq(cat(2,...
                data(s,y1,:,:,cs(c,1)),...
                data(s,y2,:,:,cs(c,1)))),...
                cat(1,results.y(y1),results.y(y2)));
            catch
                auck(s,c,k,:) = NaN;
            end
            if g
%                 aucgk(s,c,k,:,:) = colAUC(...
%                 sq(cat(2,...
%                 datag(s,:,:,:,cs(c,1),k),...
%                 datag(s,:,:,:,cs(c,1),k))),...
%                 cat(1,results.yg,results.yg));
            end
            %-- compute prediction score
            %predictk(s,c,k,:) = nanmean(sq(results.predict(s,test,:,:))==repmat(results.y(test),[1 size(results.predict,3)]),1);
            %-- avoid 1D bug in matlab
            if length(y1) == 1, y1 = [y1 y1]; end
            if length(y2) == 1, y2 = [y2 y2]; end
%             predictbk(s,c,k,:) = (...
%                 nanmean(sq(results.predict(s,y1,:,:))==repmat(results.y(y1),[1 size(results.predict,3)]),1)+...
%                 nanmean(sq(results.predict(s,y2,:,:))==repmat(results.y(y2),[1 size(results.predict,3)]),1))/2;
        end
%         %-- compute stats across folds
        switch method
            case 'ttest',  
                p = [];
%                 [h pk(s,c,:,:)]= ttest(sq(auck(s,c,:,:,:)),.5);
            case 'ttest2',  
                [h p(s,c,:,:)]  = ttest2(...
                    sq(data(s,results.y==classes(cs(c,1)),:,:,cs(c,1))),...
                    sq(data(s,results.y==classes(cs(c,2)),:,:,cs(c,1))));
%                 [h pk(s,c,:,:)]= ttest2(sq(mp1(s,c,:,1,:,:)),sq(mp2(s,c,:,2,:,:)));
            case 'ranksum', 
                [p(s,c,:,:)]  = ranksum_matrix(...
                    sq(data(s,results.y==classes(cs(c,1)),:,:,cs(c,1))),...
                    sq(data(s,results.y==classes(cs(c,2)),:,:,cs(c,1))));
%                 [pk(s,c,:,:)]  = ranksum_matrix(sq(mp(s,c,:,1,:,:)),sq(mp(s,c,:,2,:,:)));
        end
    end
end
%-- mean proba
results.mp = mp;
results.mpg = mpg;
results.mpk = mpk;
results.mpgk = mpgk;
%-- std proba
results.sp = sp;
results.spg = spg;
results.spk = spk;
results.spgk = spgk;
%-- area under the curve
results.auc = auc;
results.auck = auck;
results.p = p; % stats across trials
% results.pk = pk; % stats across folds
results.cs = cs;
% results.aucg = aucg;
% results.aucgk = aucgk;
% results.predictbk = predictbk;
% results.predictk = predictk;

fprintf('\n');