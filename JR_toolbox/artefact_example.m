% get correct paths
cd '/home/umr975/Desktop/Nas_Backups/JR/consciousness_indexes/';
tlb_path    = '~/Desktop/Nas_Backups/toolbox_tmp/';
data_path   = '~/Desktop/Nas_data/LOCAL_GLOBAL/';
addpath([tlb_path 'JR_toolbox/']);
addpath([tlb_path 'fieldtrip/fieldtrip-20110822/']);fieldtripdefs;
% get raw files
files       = dir([data_path 'controls_RAW/baptiste/b*.raw']);
% define vertical eye movements
cfg.artchan = [46 241 10 238]; 
% make files paths
cfg.dataset = cellfun(@(x) cat(2,[data_path 'controls_RAW/baptiste/'], x), {files.name}, 'UniformOutput', false);
% find artefact via PCA
comp        = artefact_find(cfg);
% make fake data
data        = rand(257,200,1000);
% correct data
data_art    = artefact_correct(comp,data);