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

% get name;
if ~isfield(cfg,'namey'),           cfg.namey = 'default';      end % classification names for outputs
if ~isfield(cfg,'path'),            cfg.path = [pwd '/'];       end % classification path for output
if ~isfield(cfg,'saveX'),           cfg.saveX= true;            end % save a new data set
if ~isfield(cfg,'savey'),           cfg.savey= true;            end % save a classification parameters
if ~isfield(cfg,'load_results'),    cfg.load_results= true;     end %
if ~isfield(cfg,'run_svm'),         cfg.run_svm= true;          end % save a new data set
if ~isfield(cfg,'Xformat'),         cfg.Xformat= '.dat';        end % 

if ischar(Xm),
    % default name identical to initial file
    if ~isfield(cfg,'nameX'),
        n = Xm;
        p = strfind(n, '/');
        if ~isempty(p), n = n((p(end)+1):end); end
        p = strfind(n, cfg.Xformat);
        if ~isempty(p), n = n(1:(p(end)-1)); end
        cfg.nameX = n;
    end
    pathX = Xm;
else
    if ~isfield(cfg,'nameX'),        cfg.nameX = 'default';    end
    pathX = [cfg.path cfg.nameX cfg.Xformat];
end

%% save data
if cfg.saveX && ~ischar(Xm) %|| ~exist(pathX, 'file')
    if ~isfield(cfg,'all_y'),
        if ~isfield(cfg,'y')
            if ~exist('y', 'var'), cfg.all_y= ones(size(Xm,1),1);
            else cfg.all_y = y;
            end
        else y = cfg.y;
        end
    end
    if ~isfield(cfg,'time'),        cfg.time= 1:size(Xm,3);     end 
    %-- info
    time = cfg.time;
    all_y = cfg.all_y;
    %-- replace NaNs
    index = sum(isnan(squeeze(nanmean(Xm,1))));
    if index~=0, warning('Changing NaNs to 0!');
        for t = 1:size(Xm,1);
            Xm(t,index) = 0;
        end
    end
    disp(['save ' cfg.Xformat ': ' pathX '...']);
    switch cfg.Xformat
        case '.dat',
            binsave(pathX, Xm);
            %-- save dimensions
            Xdim = size(Xm);
            if length(Xdim) == 2, Xdim(3) = 1; end
            save([pathX(1:(end-4)) '_dims.mat'], 'Xdim');
        case '.mat'
            save(pathX, 'Xm');
    end
end 

%% save classification parameters
if cfg.savey
    % get dimensionality
    if ~isfield(cfg,'time'),
        if ~ischar(Xm)
            cfg.time = 1:size(Xm,3);
        else
            try
                cfg.time=load(pathX,'time');cfg.time = cfg.time.time;
            catch
                load(pathX,'Xm')
                cfg.time = 1:size(Xm,3);
            end
        end
    end

    disp('Prepare classification parameters...');
    if ~isfield(cfg,'folding'),         cfg.folding= 'stratified';  end % folding type
    if strcmp(cfg.folding,'leaveoneout')
        cfg.n_folds = sum(y>0);
    else
        if ~isfield(cfg,'n_folds'),         cfg.n_folds = 5;        end % number of K stratified folding
    end
    if ~isfield(cfg,'compute_distance'),cfg.compute_distance = false;end % get fast continuous output
    if ~isfield(cfg,'compute_probas'),  cfg.compute_probas = true;  end % get slow continuous output with probability estimation
    if ~isfield(cfg,'compute_predict'), cfg.compute_predict= true;  end % get discrete output
    if ~isfield(cfg,'C'),               cfg.C = 1;                  end % svm criterion
    if ~isfield(cfg,'fs'),              cfg.fs = 1;                 end % percentile of feature selection
    if ~isfield(cfg,'n_splits'),        cfg.n_splits= 1;            end % number of n shuffle splits (only for k folding, k < n_trials)
    if ~isfield(cfg,'y2'),              cfg.y2 = y;                 end % sample weighting based on something different to y proportions
    if ~isfield(cfg,'wsize'),           cfg.wsize = 1;              end % window slide
    if ~isfield(cfg,'dims'),            
        cfg.dims = [1:(length(cfg.time)-cfg.wsize+1)];
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

pathy = [cfg.path cfg.nameX '_' cfg.namey '_y.mat'];


%% run classification
if ~isfield(cfg,'n_jobs'), n_jobs = -1;else n_jobs = cfg.n_jobs;end
if ~isfield(cfg, 'C') || length(cfg.C) == 1
    script_path = which('skl_king_parallel.py');
else
    script_path = which('skl_king_parallel_gs.py');
end
% script_path = which('skl_king_parallel_ica.py');
command = sprintf(['nice python "' script_path '" "' pathX '" "' pathy '" ' num2str(n_jobs)]);
disp(command);

if cfg.run_svm
    system(command);
end

%% load results
if cfg.load_results
    path_results = [cfg.path cfg.nameX '_' cfg.namey '_results.mat'];
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


function binsave(filename, matrix)
% binsave(filename, matrix)
    fid = fopen(filename,'w');
    fwrite(fid, matrix, 'float64');
    fclose(fid);
return