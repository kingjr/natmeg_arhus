close;
cfg             = [];
cfg.elecfile    = 'EGI_net_256.mat'; %b for rotation
layout          = ft_prepare_layout(cfg);
layout.outline  = {}; % remove outline
perimeter       = [241 242 243 247 251 256 91 102 111 120 133 145 165 174 187 199 208 216 229 233 237 240 239 238 234 10 18 25 31 32 37 46 244 241];
ear             = {[67 68 69 74 83 93 92 82 73 67], [210 219 218 217 209 201 191 192 202 210]};
layout.mask{1}  = [layout.pos(perimeter,1), layout.pos(perimeter,2)];
layout.mask{2}  = [layout.pos(ear{1},1), layout.pos(ear{1},2)];
layout.mask{3}  = [layout.pos(ear{2},1), layout.pos(ear{2},2)];
layout.outline  = layout.mask;
ft_plot_lay(layout)
save('/media/DATA/Pro/Toolbox/JR_toolbox/my_EGI_net.mat', 'layout');