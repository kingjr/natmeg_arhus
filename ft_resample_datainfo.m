function post_sampledata = ft_resample_datainfo(pre_sampledata,post_sampledata)
% post_sampledata = ft_resample_datainfo(pre_sampledata,post_sampledata)
% adjust sampleinfo after ft_resampledata function
% (c) Imen & JR
post_sampledata.sampleinfo = round(pre_sampledata.sampleinfo .* (post_sampledata.fsample / pre_sampledata.fsample));