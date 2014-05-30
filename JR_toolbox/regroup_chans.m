function [newdata rois cluster] = regroup_chans(data,options)
% [newdata rois] = regroup_chans(data,options)
% regroups channels together
% cfg.nb
% cfg.badchans
% cfg.reorder
% cfg.method
% cfg.rois
% cfg.see_cluster
% cfg.void
% cfg.index
if nargin  == 1, options.tmp  = []; end
if ~isfield(options,'nb'),              options.nb               = 16;       end
if ~isfield(options,'badchans'),        options.badchans         = [];       end
if ~isfield(options,'reorder'),         options.reorder          = false;    end % do not mean
if ~isfield(options,'method'),          options.method           = 'nanmean';   end
if ~isfield(options,'rois'),            options.rois             = 'auto';    end
if ~isfield(options,'see_cluster'),     options.see_cluster      = false;    end
if ~isfield(options,'void'),            options.void             = 0;       end
if ~isfield(options,'pnt'),           load('my_EGI_net.mat', 'layout'),options.pnt= layout.pnt; end

options.badchans        = [options.badchans 257];
%-- manual regrouping of 16 channels
rois(1).name            = 'face_left';
rois(1).index             = [241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256];
rois(2).name            = 'face_right';
rois(2).index             = [225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240];
rois(3).name            = 'front';
rois(3).index             = [38 34 28 37 33 27 32 26 31 20 25 13 19 18 12 11];
rois(4).name            = 'ftemporal_left';
rois(4).index             = [46 39 35 36 40 47 54 48 49 55 56 57 61 62 63 64 67];
rois(5).name            = 'ftemporal_right';
rois(5).index             = [10 1 2 3 4 219 220 221 222 223 211 212 213 203 204 194 224];
rois(6).name            = 'frontal';
rois(6).index             = [24 30 16 23 29 15 22 21 14 7 6 5 207 215];
rois(7).name            = 'vertex';
rois(7).index             = [9 186 132 45 81 8 198 185 144 131 90 80 53 44 17];
rois(8).name            = 'vertex_left';
rois(8).index             = [41 42 43 50 51 52 58 59 60 65 66 72 77 78 79 87];
rois(9).name            = 'vertex_right';
rois(9).index             = [197 206 214 184 196 205 183 195 182 155 164 173 143 154 163 153];
rois(10).name           = 'ptemporal_left';
rois(10).index            = [68 69 70 71 74 75 76 83 84 85 86  94 95 96 97 98];
rois(11).name           = 'ptemporal_right';
rois(11).index            = [152 161 170 178 190 162 171 179 191 172 180 192 181 193 202 210];
rois(12).name           = 'parietal';
rois(12).index            = [88 89 99 100 101 109 110 118 119 127 128 129 130 142 141 140];
rois(13).name           = 'occtemporal_left';
rois(13).index            = [92 93 103 104 105 106 107 108 112 113 114 115 116 121 122 123 134];
rois(14).name           = 'occtemporal_right';
rois(14).index            = [201 209 200 189 177 169  160 151 150 159 168 176 188 175 167 158 166];
rois(15).name           = 'occipital';
rois(15).index            = [117 126 139 125 138 124 137 149 136 148 135 147 157 146 156];
rois(16).name           ='around';
rois(16).index            = [73 82 91 102 111 120 133 145 165 174 187 199 208 216 217 218];


%-- remove bad channels
for bad_chan  = options.badchans
    for roi  = 1:length(rois)
        rois(roi).index(ismember(rois(roi).index,options.badchans))=[];
    end
end
if strcmp(options.rois, 'auto')
%-- select k clusters
switch options.nb
    case 1,
        rois(17).index    = [rois(1).index rois(2).index rois(3).index rois(4).index rois(5).index rois(6).index rois(7).index rois(8).index rois(9).index rois(10).index rois(11).index rois(12).index rois(13).index rois(14).index rois(15).index rois(16).index ];
        rois(17).index    = [rois(6).index rois(7).index rois(8).index rois(9).index rois(12).index];
        rois            = rois(17);
    case 10
        rois(17).index    = [rois(1).index rois(2).index];    % face
        rois(18).index    = [rois(4).index rois(5).index];    % ftemporal
        rois(19).index    = [rois(8).index rois(9).index];    % htemporal
        rois(20).index    = [rois(10).index rois(11).index];  % ptemporal
        rois(21).index    = [rois(13).index rois(14).index];  %.occtemporal
        rois            = rois([3 6 7 12 15 17:21]);
    case 16
        % keep manual
    case 256
        for roi        = 1:256
            rois(roi).index= roi;
        end
        rois            = rois(setdiff(1:256,options.badchans));
    otherwise
        rois            = [];
        if options.void,warning('automatic clustering of channels');end;
        idx             = kmeans(options.index,options.nb); % automatic clustering
        ids             = unique(idx);
        for roi = 1:length(ids)
            rois(roi).index = setdiff(find(idx==roi),options.badchans)';
        end
end
else
    rois = options.rois;
end

%-- duplicate one rois group to correct for stupid matlab mean function
for roi  = 1:length(rois)
    if length(rois(roi).index) == 1;
        rois(roi).index = [rois(roi).index rois(roi).index];
    end
end

%-- clean empty rois
roi = 1;
while roi <= length(rois)
    if isempty(rois(roi).index)
        rois = rois(setdiff(1:length(rois),roi));
    else
        roi = roi+1;
    end
end

if options.reorder %-- reorder
    index            = [rois.index];
    switch length(size(data))
        case 1, newdata  = data(index);
        case 2, newdata  = data(index,:);
        otherwise
            error(['data contains too many dimensions'])
    end
else                % -- regroup
    for roi  = 1:length(rois)
        switch length(size(data))
            case 1, eval(['newdata(roi)  = ' options.method '(data(rois(roi).index));']);
            case 2, eval(['newdata(roi,:)  = ' options.method '(data(rois(roi).index,:));']);
        end
    end
end

cluster = NaN(size(options.pnt,1),1);
for r = 1:length(rois)
    for c = 1:length(rois(r).index)
        cluster(rois(r).index(c)) = r;
    end
end