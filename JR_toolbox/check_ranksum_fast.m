x1 = repmat(1:8,4,1)';
x2 = 100+repmat(1:9,4,1)';
x1(1:5,4) = NaN;
x2(1:5,3) = NaN;
x2(1:2,2) = NaN;

x3=x1;
x1=x2;
x2=x3;

[p h stats] = ranksum_fast(x1,x2);

for ii = 1:4
    p2(ii) = ranksum(x1(~isnan(x1(:,ii)),ii),x2(~isnan(x2(:,ii)),ii));
    auc2(ii) = colAUC(cat(1,...
        x1(~isnan(x1(:,ii)),ii),...
        x2(~isnan(x2(:,ii)),ii)),...
        [ones(sum(~isnan(x1(:,ii))),1);...
        2*ones(sum(~isnan(x2(:,ii))),1)]);
end

for ii = 1:100
    clear x1 x2 p p2 auc auc2
    fprintf('*')
    %% build data
    x1 = randn(100,10);
    x2 = randn*randn(100,10);
    for ii =1:10
        x1(round(rand(100,1))==1,ii) = NaN;
        x2(round(rand(100,1))==1,ii) = NaN;
    end
    %% fast
    [p h stats] = ranksum_fast(x1,x2);
    auc = stats.AUC;
    %% slow
    for ii = 1:10
        p2(ii) = ranksum(x1(~isnan(x1(:,ii)),ii),x2(~isnan(x2(:,ii)),ii));
        auc2(ii) = colAUC(cat(1,...
            x1(~isnan(x1(:,ii)),ii),...
            x2(~isnan(x2(:,ii)),ii)),...
            [ones(sum(~isnan(x1(:,ii))),1);...
            2*ones(sum(~isnan(x2(:,ii))),1)]);
    end
    
    %% compare
    if sum(abs(p-p2)>10^-14) || sum(abs(auc-auc2)>10^-14)
        break
    end 
end