function results = jr_classify2(Xpath,Xdims,y,cfg)
% results = jr_classify2(Xm,y,cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # Implementation of a multivariate pattern analysis based on  the
% scikitlearn toolbox (http://scikit-learn.org/stable/). It reads a matlab
% file containing :
% #     Xm:      a matrix of trials x chans x timepoint.
% #     y:       a vector indicating the class of each trial
% # The classification algorithm is based on a support vector machine.
% # (c) Jean-Remi King 2012, jeanremi.king [at] gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 3,                     cfg = [];                   end

% get name;
if ~isfield(cfg,'nameX'),           cfg.nameX = Xpath;          end %
if ~isfield(cfg,'namey'),           cfg.namey = 'default';      end % classification names for outputs
if ~isfield(cfg,'load_results'),    cfg.load_results= true;     end %
if ~isfield(cfg,'run_svm'),         cfg.run_svm= true;          end % save a new data set

%% save classification parameters
% get dimensionality

cfg.time = 1:size(Xdims,3);

disp('Prepare classification parameters...');
if ~isfield(cfg,'folding'),         cfg.folding= 'stratified';  end % folding type
if strcmp(cfg.folding,'leaveoneout')
    cfg.n_folds = sum(y>0);
else
    if ~isfield(cfg,'n_folds'),     cfg.n_folds = 5;        end % number of K stratified folding
end
if ~isfield(cfg,'compute_distance'),cfg.compute_distance = false;end % get fast continuous output
if ~isfield(cfg,'compute_probas'),  cfg.compute_probas = true;  end % get slow continuous output with probability estimation
if ~isfield(cfg,'compute_predict'), cfg.compute_predict= true;  end % get discrete output
if ~isfield(cfg,'C'),               cfg.C = 1;                  end % svm criterion
if ~isfield(cfg,'fs'),              cfg.fs = 1;                 end % percentile of feature selection
if ~isfield(cfg,'n_splits'),        cfg.n_splits= 1;            end % number of n shuffle splits (only for k folding, k < n_trials)
if ~isfield(cfg,'y2'),              cfg.y2 = y;                 end % sample weighting based on something different to y proportions
if ~isfield(cfg,'n_jobs'),          cfg.n_jobs = -1;            end % number of cores recruited, by default leave one available
if ~isfield(cfg,'dims'),
    cfg.dims = [1:length(cfg.time)];
end 
if ~isfield(cfg,'generalize_time'), cfg.generalize_time= 'none'; end % classify on all time points?

switch cfg.generalize_time
    case 'none'
        cfg.dims_tg = cfg.dims';
    case 'all'
        cfg.dims_tg = repmat(cfg.dims,length(cfg.dims),1);
    otherwise
        if ~isfield(cfg,'dims_tg'), error('need specific time generalization matrix'); end
end

%-- diplay parameters
cfg.Xdims = Xdims;
fields = fieldnames(cfg);
for f = 1:length(fields)
    eval([fields{f} '=squeeze(cfg.' fields{f} ');']);
end
disp('Save classification parameters...');
save([cfg.nameX '_' cfg.namey '_y.mat'], 'y', fields{:}, '-v7');


cfg.pathy = [cfg.nameX '_' cfg.namey '_y.mat'];

%% run classification
if length(cfg.C) == 1
    script_path = which('skl_king_parallel2.py');
else
    script_path = which('skl_king_parallel_gs2.py');
end
% script_path = which('skl_king_parallel_ica.py');
command = sprintf(['nice python "' script_path '" "' cfg.nameX '" "' cfg.pathy '" ' num2str(cfg.n_jobs)]);
disp(command);

if cfg.run_svm
    system(command);
end

%% load results
if cfg.load_results
    path_results = [cfg.nameX '_' cfg.namey '_y_results.mat'];
    if exist(path_results,'file')
        disp(['load ' path_results '...']);
        results = load(path_results);
    else
        disp([path_results 'does not exist.']);
    end
else
    results = cfg;
end
return