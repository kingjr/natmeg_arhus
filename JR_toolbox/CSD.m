% CSD - Current Source Density (CSD) transformation based on spherical spline
%       surface Laplacian as suggested by Perrin et al. (1989, 1990)
%
% Usage: [X, Y] = CSD(Data, G, H, cfg);
%
% input:    Data = nchan x nsample
%           G & H: nchan x nchan matrix
%           cfg.lamda:      default: 1.0e-5
%           cfg.head:       default: 1
%           cfg.memory:     default: 'medium'
% output:
%           X: nchan x nsample CSD data
%
% Input parameters:
%   Data = surface potential electrodes-by-samples matrix
%      G = g-function electrodes-by-electrodes matrix
%      H = h-function electrodes-by-electrodes matrix
% lambda = smoothing constant lambda (default = 1.0e-5)
%   head = head radius (default = no value for unit sphere [�V/m�])
%          specify a value [cm] to rescale CSD data to smaller units [�V/cm�]
%          (e.g., use 10.0 to scale to more realistic head size)
%
% Output parameter:
%      X = current source density (CSD) transform electrodes-by-samples matrix
%      Y = spherical spline surface potential (SP) interpolation electrodes-
%          by-samples matrix (only if requested)
%
%
%
% Updates:
%   - JR KING 2011 10 07: fix bug in Y calculation, change scaling method
%   - JR KING 2011 10 06: fix cfg.lambda bug + memory parameter
%   - JR KING 2011 10 05: fasten procedure
%--------------------------------------------------------------------------
% adapted by Jean-Remi King from:
% (published in appendix of Kayser J, Tenke CE, Clin Neurophysiol 2006;117(2):348-368)
% Implementation of algorithms described by Perrin, Pernier, Bertrand, and
% Echallier in Electroenceph Clin Neurophysiol 1989;72(2):184-187, and
% Corrigenda EEG 02274 in Electroenceph Clin Neurophysiol 1990;76:565.
%
% Copyright (C) 2003 by J�rgen Kayser (Email: kayserj@pi.cpmc.columbia.edu)
% GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
% Updated: $Date: 2005/02/11 14:00:00 $ $Author: jk $
%        - code compression and comments
% Updated: $Date: 2007/02/07 11:30:00 $ $Author: jk $
%        - recommented rescaling (unit sphere [�V/m�] to realistic head size [�V/cm�])
%   Fixed: $Date: 2009/05/16 11:55:00 $ $Author: jk $
%        - memory claim for output matrices used inappropriate G and H dimensions
%   Added: $Date: 2009/05/21 10:52:00 $ $Author: jk $
%        - error checking of input matrix dimensions
%--------------------------------------------------------------------------

function [X, Y] = CSD(Data, G, H, cfg)
if nargin < 3, error('wrong number of argument');   end
if nargin == 3,             cfg         = [];       end
if ~isfield(cfg,'lambda'),  cfg.lambda  = 1.0e-5;   end
if ~isfield(cfg,'head'),    cfg.head    = 1.0;      end
if ~isfield(cfg,'memory'),  cfg.memory  = 'medium'; end
[nElec,nPnts] = size(Data);                                                 % get data matrix dimensions
head = cfg.head ^2;                                                         % or rescale data to head sphere [�V/cm�]
mu          = mean(Data);                                                   % get grand mean
Z           = (Data - repmat(mu,nElec,1));                                  % compute average reference
if nargin < 5; head = 1.0; end;                                             % initialize scaling variable [�V/m�]
head        = head^2;                                                       % or rescale data to head sphere [�V/cm�]
if nargin < 4; cfg.lambda = 1.0e-5; end;                                        % initialize smoothing constant
for e = 1:size(G,1);                                                        % add smoothing constant to diagonale
    G(e,e)  = G(e,e) + cfg.lambda;
end;
Gi          = inv(G);                                                       % compute G inverse
TC          = sum(Gi);
sgi         = sum(TC);                                                      % compute sum total
[X Y]       = deal(NaN(size(Data)));                                        % initialize
scale       = median(abs(Z(:)))+eps;                                        % rescale to avoid very small numbers
Z           = Z/scale;
Cp          = (Gi * Z)';                                                    % compute preliminary C vector
c0          = sum(Cp,2) / sgi;                                              % common constant across electrodes
C           = Cp - (repmat(c0,1,nElec) .* repmat(TC,nPnts,1));              % compute final C vector
switch cfg.memory
    case 'low'                                                              %----------old algorithm: slow but works on all platforms
        for p = 1:nPnts
            Cp = Gi * Z(:,p);                                               % compute preliminary C vector
            c0 = sum(Cp) / sgi;                                             % common constant across electrodes
            C = Cp - (c0 * TC');                                            % compute final C vector
            parfor e = 1:nElec;                                             % compute all CSDs ...
                X(e,p) = sum(C .* H(e,:)') / head;                          % ... and scale to head size
            end
           if nargout > 1; 
                parfor e = 1:nElec;                                         % if requested ...
                    Y(e,p) = c0 + sum(C .* G(e,:)');                        % ... compute all SPs
                end
            end
        end
    case 'medium'                                                           %---------- new algo: faster but more memory
        parfor e = 1:nElec;                                                 
            X(e,:) = sum(C .* repmat(H(e,:), [nPnts,1]),2);
        end
        X       = X'/head;                                                  % ... and scale to head size
        if nargout > 1
            parfor e = 1:nElec;                                             
                Y(e,:) = c0 + sum(C .* repmat(G(e,:),[nPnts 1]),2);         
            end
        end
    case 'high'                                                             %---------- super fast but rapidly requires very large RAM
        X       = squeeze(sum(repmat(C,[1 1 nElec]).*permute(repmat(H,[1 1 nPnts]),[3 1 2]),2))/head;
        if nargout > 1;
            Y = repmat(c0',[nElec 1]) + squeeze(sum(repmat(C,[1 1 nElec]).*permute(repmat(G,[1 1 nPnts]),[3 1 2]),2))'; % to be check!!
        end
    otherwise, error('unknown method');
end
X               = X' * scale;                                               % back to normal scale
Y               = Y * scale;
end

