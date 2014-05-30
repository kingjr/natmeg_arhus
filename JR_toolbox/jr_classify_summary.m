function mvpa=jr_classify_summary(data)
%% load data
if ischar(data)
    load(data, 'probas', 'predict', 'y', 'all_folds', 'probasg', 'predictg', 'yg', 'yall');
else
    fields = fieldnames(data);
    for f = 1:length(fields)
        eval([fields{f} '=data.' fields{f} ';']);
    end
end
try mvpa.probas = probas; end
try mvpa.predict = predict; end
try mvpa.all_folds = all_folds; end
try mvpa.y = y; end
try mvpa.probasg = probasg;end
try mvpa.predictg = predictg; end
try mvpa.yg = yg; end
try mvpa.yall = yall; end
%% classification parameters
classes = unique(y);
nc = length(classes);
nchans = size(probas,3);

if isempty(predict), predict = NaN*probas; end
%% condition probability
x = NaN(length(y),nchans);
for c = 1:nc
    % for each class
    mvpa.probas_class(c,:)  = nanmean(nanmean(probas(:,y==classes(c),:,c),1),2);
    mvpa.predict_class(c,:) = nanmean(squeeze(round(nanmean(predict(:,y==classes(c),:),1))==c));
    % across class
    mvpa.probas_summary(c,:) = nansum(nanmean(probas(:,y==classes(c),:,c),1),2);
    mvpa.predict_summary(c,:) = nansum(squeeze(round(nanmean(predict(:,y==classes(c),:),1))==c));
    % sorted
    mvpa.probas_sorted(:,y==classes(c),:,:) =  probas(:,y==classes(c),:,[c:nc 1:(c-1)]);
end
mvpa.probas_summary = nansum(mvpa.probas_summary)/length(y);
mvpa.predict_summary = nansum(mvpa.predict_summary)/length(y);

%% statistical comparison across trials
cs = nchoosek(1:nc,2); % comparisons
mvpa.comparison = cs;
for c = 1:size(cs,1)
    %-- normalize proba
    x = nanmean(probas,1); % mean across splits
    x = x(:,:,:,:,[cs(c,1) cs(c,2)]);
    x = squeeze(x ./ repmat(sum(x,5),[1 1 1 1 2]));
    %-- stat
    switch 'ttest2'
        case 'ttest2'
            [unused p unused stats] = ttest2(...
                x(y==classes(cs(c,1)),:,1),...
                x(y==classes(cs(c,2)),:,1));
            mvpa.p(c,:) = p;
            mvpa.tstat(c,:) = stats.tstat;
        case 'ranksum'
            [p h stats] = ranksum_matrix(...
                x(y==classes(cs(c,1)),:,1),...
                x(y==classes(cs(c,2)),:,1));
            mvpa.p(c,:) = p;
            mvpa.tstat(c,:) = [stats.zval];
    end
    %-- mean proba diff
    mvpa.proba_diff(c,:) = nanmean(x(y==classes(cs(c,1)),:,1),1) - nanmean(x(y==classes(cs(c,2)),:,1),1);
end


%% statistical comparison across folds
cs = nchoosek(1:nc,2); % comparisons
mvpa.comparison = cs;
for c = 1:size(cs,1)
    for s = size(probas,1):-1:1
        for f = size(all_folds,2):-1:1
            %-- normalize proba
            x = probas(s,:,:,:,[cs(c,1) cs(c,2)]);
            x = x ./ repmat(sum(x,5),[1 1 1 1 2]);
            %-- avt regroup within folds
            x1(f,:) = squeeze(nanmean(x(:,intersect(find(y==classes(cs(c,1))),find(all_folds(s,f,:)==0)),:,:,1),2));
            x2(f,:) = squeeze(nanmean(x(:,intersect(find(y==classes(cs(c,2))),find(all_folds(s,f,:)==0)),:,:,1),2));
            %-- mean proba diff
            mvpa.proba_diff_fold(c,s,f,:) = x1(f,:)-x2(f,:);
        end
        %-- stat
        [unused p unused stats] = ttest2(x1,x2);
        mvpa.p_fold(c,s,:) = p;
        mvpa.tstat_fold(c,s,:) = stats.tstat;
    end
end