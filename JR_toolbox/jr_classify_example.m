%% decoding example using the scikit learn pipeline for Stan
% 1. Install scikit learn toolbox
% 1.1 add Neuro-debian repesitory
% go to http://neuro.debian.net/
% > How to use this repository => select a release
% > inter terminal:
%       paste provided command (wget ...)
%       sudo apt-get update
%       sudo apt-get install python-sklearn
%% 1.2 fast example
clear;
X               = randn(50,20,100);           % matrix of 200 trials by 300 channels by 100 time points)
y               = [ones(25,1);2*ones(25,1)];  % vector indicating trial class (i.e. labels)
X(y==1,:,50:100)= X(y==1,:,50:100) + .3;         % add information at time t=[5 10]
results         = jr_classify(X,y);         % decode across sensors at each time point
%results         = jr_classify_stats(results);   % summarize results within AUC
output          = squeeze(results.probas(1,:,:,1,1));
results.auc     = colAUC(output, results.y);
plot(squeeze(results.auc));                     % plot area under the curve

cfg = [];
cfg.wsize       = 10;
results         = jr_classify(X,y,cfg);         % decode across sensors at each time point
results         = jr_classify_stats(results);   % summarize results within AUC
plot(squeeze(results.auc));                     % plot area under the curve


%% 1.3 more detailed example
nsamples    = 150;
nchans      = 256;
ntrials     = 50;
% matrix of data
X           = zeros(ntrials,nchans,nsamples);
% /!\ matlab automatically squeezes matrices when the last dimension has
% only one column. At the moment, it is therefore necessary to have at
% least 2 time point / 2 datavector to be classified

% vector of class: (i.e. labels)
%    y==0: trials that will not be taken into account. 
%    y>0: trials that will be used in a traditional decoding method
%    y<0: trials that will be used in a generalization procedure (decode
%    pattern with y>0 and see whether this pattern generalizes to another
%    conditions).
y = [ones(ntrials/2,1); 2*ones(ntrials/2,1)];

%-- add differential information between t = [30 50] (i.e. chain of
%processing) on a fourth of the channels
for t = 30:50
    for chan = 1:4:nchans
        X(y==1,chan,t) = randn;
        X(y==2,chan,t) = randn;
    end
end
%-- add common information between t = [100:130] on a fourth of the channels
for chan = 1:4:nchans
    X(y==1,chan,100:130) = randn;
    X(y==2,chan,100:130) = randn;
end

%-- plot
figure(1);clf;set(gcf,'color','w');
subplot(3,4,1);
imagesc(squeeze(mean(X(y==1,:,:)))); title('class A');ylabel('chans');xlabel('time');
subplot(3,4,2);
imagesc(squeeze(mean(X(y==2,:,:)))); title('class B');ylabel('chans');xlabel('time');
subplot(3,4,3);
imagesc(squeeze(nanmean(X(y==1,:,:)))-squeeze(nanmean(X(y==2,:,:)))); title('A-B');ylabel('chans');xlabel('time');

%-- add gaussian noise
X = X + 4*randn(size(X));

%-- plot
subplot(3,4,5);
imagesc(squeeze(mean(X(y==1,:,:)))); title('class A');ylabel('chans');xlabel('time');
subplot(3,4,6);
imagesc(squeeze(mean(X(y==2,:,:)))); title('class B');ylabel('chans');xlabel('time');
subplot(3,4,7);
imagesc(squeeze(nanmean(X(y==1,:,:)))-squeeze(nanmean(X(y==2,:,:)))); title('A-B');ylabel('chans');xlabel('time');
subplot(3,4,8);
% parametric stat s
[h p] = ttest2(X(y==1,:,:),X(y==2,:,:));
% nonparametric stat s
p = ranksum(X(y==1,:,:),X(y==2,:,:));
imagesc(-log10(squeeze(p)));title('-log10(p)');ylabel('chans');xlabel('time');colorbar


%-- define optional parameters. See jr_classify to see all possible
%parameters
cfg         = [];
cfg.n_folds = 10; % number of folds
cfg.fs      = .50; % percentage of feature selection
cfg.C       = 1;  % SVM criterion 
% save data once
cfg.nameX           = 'mydata'; % name of data file
cfg.run_svm         = false;
cfg.load_results    = false;
jr_classify(X,y,cfg);
% run decoding 
cfg.run_svm         = true;
cfg.load_results    = true;
cfg.namey           = 'myclf'; % name of classification parameters file
jr_classify('mydata',y,cfg);
results = load('mydata_myclf_results.mat');

%-- compute area under the curve (i.e. the effect size)
auc = [];
for split = 1:size(results.probas,1) % by default split = 1
    for decoding_time = 1:size(results.probas,3) 
        for generalization_time = 1:size(results.probas,4) % by default there is no generalization time
            auc(split,decoding_time,generalization_time) = colAUC(...
                squeeze(results.probas(split,:,decoding_time,generalization_time,1))',...
                results.y);
        end
    end
end
% note that colAUC works on matrices, so you can remove these loops...

% run summarizing statistics function on decoding results
results=jr_classify_stats(results);
auc = squeeze(results.auc);
% deltap = squeeze(nanmean(results.mpk,3));
% subplot(3,4,9);plot(deltap);title('P(A|A) - P(A|B)');
subplot(3,4,10);plot(auc);title('AUC');
subplot(3,4,11);plot(-log10(squeeze(results.p)));title('-log10(p)');

%% generalization across conditions
%% generalization across time