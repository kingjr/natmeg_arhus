function results = jr_classify(Xm,y,cfg)
% results = jr_classify(Xm,y,cfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # Implementation of a multivariate pattern analysis based on  the
% scikitlearn toolbox (http://scikit-learn.org/stable/). It reads a matlab
% file containing :
% #     Xm:      a matrix of trials x chans x timepoint.
% #     y:       a vector indicating the class of each trial
% # The classification algorithm is based on a support vector machine.
% # (c) Jean-Remi King 2012, jeanremi.king [at] gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2,                     cfg = [];                   end
if ischar(Xm),                      cfg.nameX = Xm;             end
if ~isfield(cfg,'nameX'),           cfg.nameX = 'default';      end % classification names for outputs
if ~isfield(cfg,'namey'),           cfg.namey = 'default';      end % classification names for outputs
if ~isfield(cfg,'path'),            cfg.path = [pwd '/'];       end % classification path for output
if strcmp(cfg.nameX(1), '/')
    pathX = [cfg.nameX];
else
    pathX = [cfg.path cfg.nameX];
end
if strcmp(cfg.nameX((end-4):end), '.mat');
    pathX = [pathX '.mat'];
end
if ~isfield(cfg,'saveX'),           cfg.saveX= false;           end % save a new data set
if ~isfield(cfg,'savey'),           cfg.savey= true;            end % save a classification parameters
if ~isfield(cfg,'load_results'),    cfg.load_results= true;     end %
if ~isfield(cfg,'run_svm'),         cfg.run_svm= true;          end % save a new data set
 


%% save data
if isempty(dir(pathX)) || (cfg.saveX && ~ischar(Xm))
    if ~isfield(cfg,'all_y'),       cfg.all_y= ones(size(Xm,1),1);end 
    if ~isfield(cfg,'time'),        cfg.time= 1:size(Xm,3);     end 
    
    %-- replace NaNs
    index = isnan(squeeze(nanmean(Xm,1)));
    if ~isempty(index), warning('Changing NaNs to 0!'); end
    for t = 1:size(Xm,1);
        Xm(t,index) = 0;
    end
    disp(['save ' pathX '...']);
    
    save(pathX, 'Xm', 'time', 'all_y');
end 

%% save classification parameters
if cfg.savey
    % get dimensionality
    if ~isfield(cfg,'time'),
        cfg.time=load(pathX,'time');cfg.time = cfg.time.time;
    end

    disp('Prepare classification parameters...');
    if ~isfield(cfg,'folding'),         cfg.folding= 'stratified';  end % folding type
    if ~isfield(cfg,'compute_probas'),  cfg.compute_probas = true;  end % get continuous output
    if ~isfield(cfg,'compute_predict'), cfg.compute_predict= true;  end % get discrete output
    if ~isfield(cfg,'C'),               cfg.C = 1;                  end % svm criterion
    if ~isfield(cfg,'fs'),              cfg.fs = 99;                end % percentile of feature selection
    if ~isfield(cfg,'n_folds'),         cfg.n_folds = 5;            end % number of K stratified folding
    if ~isfield(cfg,'n_splits'),        cfg.n_splits= 1;            end % number of n shuffle splits (only for k folding, k < n_trials)
    if ~isfield(cfg,'y2'),              cfg.y2 = y;                 end % sample weighting based on something different to y proportions
    if ~isfield(cfg,'dims'),            
        cfg.dims = [1:length(cfg.time)];
    end % classify on all time points
    if ~isfield(cfg,'generalize_time'), cfg.generalize_time= 'none'; end % folding type
    
    switch cfg.generalize_time
        case 'none'
            cfg.dims_tg = cfg.dims';
        case 'all'
            cfg.dims_tg = repmat(cfg.dims,length(cfg.dims),1);
        otherwise 
            if ~isfield(cfg,'dims_tg'), error('need specific time generalization matrix'); end
    end
    
    %-- diplay parameters
    disp(cfg);
    
    fields = fieldnames(cfg);
    for f = 1:length(fields)
        eval([fields{f} '=squeeze(cfg.' fields{f} ');']);
    end
    disp('Save classification parameters...');
    save([cfg.path cfg.nameX '_' cfg.namey '_y.mat'], 'y', fields{:}, '-v7');
end

%% run classification
if cfg.run_svm
    script_path = which('skl_king.py');
    command = ['python ' script_path ' ' cfg.path cfg.nameX '.mat ' cfg.path cfg.nameX '_' cfg.namey '_y.mat'];
    disp(command);
    system(command);
end

%% load results
if cfg.load_results
    path_results = [cfg.path cfg.nameX '_' cfg.namey '_results.mat'];
    disp(['load ' path_results '...']);
    results = load(path_results);
else
    results = cfg;
end