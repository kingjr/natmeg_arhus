function out = WSMI(X,dict,varargin)
% out = wsmi(X,dict,[weight_type],[[across])
% input: 
%       X: data matrix: channels x symbols x trials
%       dict: symbol dictionary. 
%       [across]: vector of pairs of signals to compute
%       [weight_type]: weighting matrix: 'identical', 'opposite'
%       [order]: mask order
%       [compute]: 'slow'. 'fast' requires very high RAM.
%       [keep_PXY]: false
% output:
%       out.across: pairs of signals computed
%       out.PXY: conjunction matrix (pairs x trials x symbol x symbol)
%       out.PX: probability of each symbol (signal x trials x symbols)
%       out.wsmi: weighted symbolic mutual information (pairs x trials)

%% default parameters
if nargin == 1, varargin = {}; end
for ii = 1:2:length(varargin)
    eval([varargin{ii} '=varargin{ii+1};']);
end
if ~exist('across','var'),  across = nchoosek(1:size(X,1),2); end

%% build weighting mask
[ns kernel] = size(dict);
if ~exist('mask', 'var')
    if ~exist('weight_type','var'),weight_type = 'opposite'; end
    switch weight_type
        case {'opposite' , 'identical'}
            mask = NaN(ns);
            for s1 = 1:ns
                for s2 = s1:ns
                    [mask(s1,s2) mask(s2,s1)] = deal(sum(abs(dict(s1,:)-dict(s2,:)).^2));
                end
            end
            if strcmp(mask,'opposite')
                % compute max permutation
                m = sqrt(sum(abs((1:kernel)-(kernel:-1:1)).^2));
                % apply inverse
                mask = squeeze(min(cat(3,mask,m-mask'),[],3));
            end
        otherwise
            mask = inf(ns);
    end
    if ~exist('order','var'),   order = 1; end
    if ~isempty(order)
        mask = mask>order;
    end
end

%% collapse repetition dimension
sz = size(X);
X = reshape(X,size(X,1),size(X,2),[]);
nr = size(X,3); % number of repetition
nsig = size(X,1); % number of signals
nt = size(X,2); % number of time samples
na = size(across,1); % number of pair of signal
nsym = size(dict,1); % number of distinct symbols

%% Compute probabilities
%wsmi = uint(zeros(na,nr));
wsmi = uint8(zeros(na,nr));

%% expand matrix for faster search of conjunctions
xy = false(nsig,nr,nt,ns);
PX = zeros(nsig,nr,ns);
for r = nr:-1:1
    for sig = nsig:-1:1
        for sym = nsym:-1:1
            sel = X(sig,:,r)==sym;
            % expand
            xy(sig,r,sel,sym) = true;
            % prior probability
            PX(sig,r,sym) = sum(sel);
        end
    end
end
PX = double(PX./repmat(sum(PX,3),[1 1 nsym]));

%% Build conjunction matrix
PXY = zeros(na,nr,ns,ns);
if ~exist('compute', 'var'), compute = 'slow'; end
switch compute
    case 'slow'
        for s1 = ns:-1:1
            for s2 = ns:-1:1
                PXY(:,:,s1,s2) = sum(xy(across(:,1),:,:,s1) & xy(across(:,2),:,:,s2),3);
            end
        end
    case 'medium'
        for s2 = ns:-1:1
            PXY(:,:,:,s2) = sum(...
                repmat(xy(across(:,1),:,:,s1),[1 1 1 ns]) & ...
                xy(across(:,2),:,:,:),3);
        end
    case 'fast'
        PXY = uint(sum(...
                repmat(permute(xy(across(:,1),:,:,:),[1 2 3 5 4]), [1 1 1 ns 1]).*...
                repmat(xy(across(:,2),:,:,:), [1 1 1 1 ns]),3));
end
% normalize
PXY = double(PXY/nt);

%% Compute WMSI
wsmi = zeros(na,nr);
smi = zeros(na,nr);
for s1 = 1:ns
    for s2 = 1:ns
        smi_ = PXY(:,:,s1,s2) .*...
            log(...
            PXY(:,:,s1,s2) ./ ...
            squeeze(...
            PX(across(:,1),:,s1).*PX(across(:,2),:,s2)));
        smi_(isnan(smi_(:))) = 0;
        smi = smi + smi_;
        wsmi = wsmi + mask(s1,s2)*smi;
    end
end
wsmi = wsmi/log(ns);
smi = smi/log(ns);


%% out
out.across = across;
out.PXY = reshape(PXY,[na,sz(3:end) ns,ns]);
out.PX = reshape(PX,[nsig,sz(3:end),ns]);
out.smi = reshape(smi,[na,sz(3:end),1]);
out.wsmi = reshape(wsmi,[na,sz(3:end),1]);