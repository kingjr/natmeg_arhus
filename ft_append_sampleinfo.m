function data_concatenated = ft_append_sampleinfo(cfg,data)
if ~isfield(cfg,'add'), cfg.add = [0 0]; end % number of sample to add between run to the sampleinfo field
% take cfg and data{n} containing n ft structure to be appended
% adapt final sample info by adding the time spent on the previous run(s)
% (c) JR King 2011

%-- concatenate data
cfg.eval_runs        = '';
for run = 1:length(data)
    cfg.eval_runs    = [cfg.eval_runs ', data{' num2str(run) '}'];
end
data_concatenated           = eval(['ft_appenddata(cfg' cfg.eval_runs ');']);

if isfield(data{1},'grad'), data_concatenated.grad = data{1}.grad; end 

%-- concatenate sample info
try
data_concatenated.sampleinfo = data{1}.sampleinfo;
for run = 2:length(data)
    data_concatenated.sampleinfo = cat(1,...
        data_concatenated.sampleinfo,...
        data{run}.sampleinfo + ...
        repmat(max(data_concatenated.sampleinfo) + cfg.add,size(data{run}.sampleinfo,1),1));
end
catch
warning('could not get sampleinfo');
end
return