function results = jr_classify(Xm,y,cfg)
% results = jr_classify(Xm,y,cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # Implementation of a multivariate pattern analysis based on  the
% scikitlearn toolbox (http://scikit-learn.org/stable/). It reads a matlab
% file containing :
% #     Xm:      a matrix of trials x chans x timepoint.
% #     y:       a vector indicating the class of each trial
% # The classification algorithm is based on a support vector machine.
% Optional inputs
% cfg.dims = [1:nt]';
% cfg.dims_tg = repmat(cfg.dims',length(cfg.dims),1)
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
if ~isfield(cfg,'void'),            cfg.void = true;            end % provide continuous feedback

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
    cfg.Xdim = size(Xm);
end

%% save data
if cfg.saveX && ~ischar(Xm) %|| ~exist(pathX, 'file')
    if cfg.void, disp(['save ' cfg.Xformat ': ' pathX '...']);end
    switch cfg.Xformat
        case '.dat',
            cfg.Xdim = size(Xm);
            binsave(pathX, Xm);
        case '.mat'
            save(pathX, 'Xm');
    end
end 

%% save classification parameters
if cfg.savey
    if length(cfg.Xdim)==2
        cfg.Xdim(3) = 1;
    end
%         cfg.dims = [1:(length(cfg.time)-cfg.wsize+1)];
%         cfg.dims_tg = cfg.dims';
%         cfg.dims_tg = repmat(cfg.dims,length(cfg.dims),1);

    %-- diplay parameters
    if cfg.void
        disp(cfg);
    end
    
    fields = fieldnames(cfg);
    for f = 1:length(fields)
        eval([fields{f} '=squeeze(cfg.' fields{f} ');']);
    end
    if cfg.void, disp('Save classification parameters...'); end
    save([cfg.path cfg.nameX '_' cfg.namey '_y.mat'], 'y', fields{:}, '-v7');
end

pathy = [cfg.path cfg.nameX '_' cfg.namey '_y.mat'];


%% run classification
if ~isfield(cfg,'n_jobs'), n_jobs = -1;else n_jobs = cfg.n_jobs;end
script_path = which('skl_king_parallel.py');

if isunix
    command = sprintf([...
        'nice python "' script_path...
        '" "' pathX ...
        '" "' pathy ...
        '" ' num2str(n_jobs)]);
else
    command = [...
        'python "' script_path...
        '" "' pathX ...
        '" "' pathy ...
        '" ' num2str(1)];
    warning('no parallel computing is supported with windows');
end
if cfg.void 
    disp(command);
end

if cfg.run_svm
    if cfg.void, system(command);
    else evalc('system(command)'); end
end

%% load results
if cfg.load_results
    path_results = [cfg.path cfg.nameX '_' cfg.namey '_results.mat'];
    if exist(path_results,'file')
        if cfg.void 
            disp(['load ' path_results '...']);
        end
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
    fwrite(fid, matrix, 'float32');
    fclose(fid);
return