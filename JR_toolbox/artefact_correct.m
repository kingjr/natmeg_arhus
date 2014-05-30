function data = artefact_correct(cfg,data)
% data = correct_artefact(cfg,data)
% -------------------------------------------------------------------------
% apply artefact correction from PCA (see find_artefact function)
% -------------------------------------------------------------------------
% input:
%       - data:             nchan x nsample x ntrial
%
%       - cfg.all_comp      nchan x nchan matrix of all components
%       - cfg.clear_comp    nchan x nchan matrix of non artefacted components
% output:
%       - data
% -------------------------------------------------------------------------
% (c) Jean-Rémi King 2011, all rights reserved
% jeanremi.king+matlab@gmail.com
% -------------------------------------------------------------------------

if ~nargin == 2,                                            error('missing fields');        end
if ~isfield(cfg,{'chantypes', 'all_comps', 'clear_comps'}), error('missing cfg fields');    end

[nchan nsample ntrial] = size(data);
% reshape to fasten operation
data = reshape(data,nchan, []);
% apply correction per chan type
for chantype = 1:length(cfg.chantypes)
    data(cfg.chantypes{chantype},:) = ...
        cfg.all_comp{chantype}*...
        cfg.clear_comp{chantype}*...
        data(cfg.chantypes{chantype},:);
end
% give data back in its initial format
data = reshape(data,nchan, nsample, ntrial);
return