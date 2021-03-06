function [p, h, stats] = ranksum_fast(x,y,varargin)
%RANKSUM Wilcoxon rank sum test for equal medians.
%   P = RANKSUM(X,Y) performs a two-sided rank sum test of the hypothesis
%   that two independent samples, in the vectors X and Y, come from
%   distributions with equal medians, and returns the p-value from the
%   test.  P is the probability of observing the given result, or one more
%   extreme, by chance if the null hypothesis ("medians are equal") is
%   true.  Small values of P cast doubt on the validity of the null
%   hypothesis.  The two sets of data are assumed to come from continuous
%   distributions that are identical except possibly for a location shift,
%   but are otherwise arbitrary.  X and Y can be different lengths.
%   The two-sided p-value is computed by doubling the most significant
%   one-sided value.
%
%   The Wilcoxon rank sum test is equivalent to the Mann-Whitney U test.
%
%   [P,H] = RANKSUM(...) returns the result of the hypothesis test,
%   performed at the 0.05 significance level, in H.  H=0 indicates that
%   the null hypothesis ("medians are equal") cannot be rejected at the 5%
%   level. H=1 indicates that the null hypothesis can be rejected at the
%   5% level.
%
%   [P,H] = RANKSUM(...,'alpha',ALPHA) returns the result of the hypothesis
%   test performed at the significance level ALPHA.
%
%   [P,H] = RANKSUM(...,'method',M) computes the p-value exactly if M is
%   'exact', or uses a normal approximation if M is 'approximate'.  If you
%   omit this argument, RANKSUM uses the exact method for small samples and
%   the approximate method for larger samples.
%
%   [P,H,STATS] = RANKSUM(...) returns STATS, a structure with one or two
%   fields.  The field 'ranksum' contains the value of the rank sum
%   statistic.  For the 'approximate' method, the field 'zval' contains the
%   value of the normal (Z) statistic.
%
%   See also SIGNTEST, SIGNRANK, KRUSKALWALLIS, TTEST2.

%   References:
%      [1] Hollander, M. and D. A. Wolfe.  Nonparametric Statistical
%          Methods. Wiley, 1973.
%      [2] Gibbons, J.D.  Nonparametric Statistical Inference,
%          2nd ed.  M. Dekker, 1985.

%   Copyright 1993-2009 The MathWorks, Inc. 
%   $Revision: 1.14.4.6 $
%--------------------------------------------------------------------------
% JR King modification for fast matrix computation + addup of AUC, R and U
% values, deal with NaN

% Check most of the inputs now
alpha = 0.05;
if nargin>2 && isnumeric(varargin{1})
   % Grandfathered syntax:  ranksum(x,y,alpha)
   alpha = varargin{1};
   varargin(1) = [];
end
oknames = {'alpha' 'method'};
dflts   = {alpha   ''};
[eid,emsg,alpha,method] = internal.stats.getargs(oknames,dflts,varargin{:});
if ~isempty(eid)
   error(sprintf('stats:ranksum:%s',eid),emsg);
end

if ~isscalar(alpha)
   error('stats:ranksum:BadAlpha','RANKSUM requires a scalar ALPHA value.');
end
if ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
   error('stats:ranksum:BadAlpha','RANKSUM requires 0 < ALPHA < 1.');
end



% Get the samples and their sizes, find the larger sample
sx = size(x);
sy = size(y);
x = reshape(x,size(x,1),[]);
y = reshape(y,size(y,1),[]);
nd = size(x,2);

% remove common NaN
x = x(mean(isnan(x(:,:)),2)<1,:);
y = y(mean(isnan(y(:,:)),2)<1,:);


if sum(sx~=sy)>1
   error('stats:ranksum:InvalidData',...
         'RANKSUM_FAST requires identical matrices');
end
nx = numel(x(:,1));
nnx = nx-sum(isnan(x),1); % count specific NaN
ny = numel(y(:,1));
nny = ny-sum(isnan(y),1); % count specific NaN

ns = min(nx,ny);
nS = max(nx,ny);

nns = min(nnx,nny);

%% always put smallest sample size first
[smsample lgsample] = deal(NaN(nS,nd)); 
smsample(1:nx,nnx<=nny)= x(:,nnx<=nny);
smsample(1:ny,nnx>nny) = y(:,nnx>nny);
lgsample(1:ny,nnx<=nny)= y(:,nnx<=nny);
lgsample(1:nx,nnx>nny) = x(:,nnx>nny);

% Now deal with the method argument
if isempty(method)
    if length(nnx)~=1 && (length(unique(nnx))>1 || length(unique(nny))>1)
        if max(nns)<10 && max((nnx+nny))<20         
            method = 'exact';
        elseif min(nns)<10 && min((nnx+nny))<20
            method = 'approximate';
        else
            warning('inappropriate method in some dimensions')
            method = 'approximate';
        end
    else
        if ns<10 && (nx+ny)<20
            method = 'exact';
        else
            method = 'approximate';
        end
    end
end

% Compute the rank sum statistic based on the smaller sample
[ranks, tieadj] = tiedrank([smsample; lgsample]);
w = nansum(ranks(1:nS,:),1);
% w = NaN(1,size(ranks,2));
% w(nnx<=nny) = nansum(ranks(1:ns,nnx<=nny),1);
% w(nnx>nny) = nansum(ranks(ns+1:end,nnx>nny),1);

% JR addup: see There is an alternative formulation of this test that yields a statistic
%commonly denoted by U. U is related to T by the formula U=T-k*(k+1)/2,
%where k is the size of the smaller sample (or either sample if both contain
%the same number of individuals). For a presentation of the U statistic, see:
% S. Siegel and N. J. Castellan, Jr., Nonparametric Statistics for the
% Behavioral Sciences, 2d ed. McGraw-Hill, New York, 1988, Section 6.4, “The
% Wilcoxon-Mann-Whitney U Test.”
%For a detailed derivation and discussion of the Mann-Whitney test as developed
%here, as well as its relationship to U, see:
% F. Mosteller and R. Rourke, Sturdy Statistics: Nonparametrics and Order
% Statistics, Addison-Wesley, Reading, MA, 1973, Chapter 3, “Ranking Methods for
% Two Independent Samples.”

%stats.U = w-ns*(ns+1)/2;
stats.U = w-nns.*(nns+1)/2;

% JR addup: see Cardillo G. (2009). MWWTEST: Mann-Whitney-Wilcoxon non parametric test for two unpaired samples.
% http://www.mathworks.com/matlabcentral/fileexchange/25830

% stats.AUC = stats.U ./ (nx.*ny);
stats.AUC = stats.U ./ (nnx.*nny);
stats.AUC(nnx>nny) = 1-stats.AUC(nnx>nny);


wmean = nns.*(nnx + nny + 1)/2;
if isequal(method,'exact')    % use the sampling distribution of W
      % For small samples, enumerate all possibilities, find smaller tail prob
      [un ulast uid] = unique(nns);
      for u = 1:length(un)
          % find corresponding dim
          dims = find(nns==un(u));
          if un(u) >0
              % non nan ranks
              nnranks = ranks(~isnan(ranks(:,dims(1))),dims(1));
              % compute distribution once per unique degree of freedom (determined
              % by the number of NaNs)
              allpos = nchoosek(nnranks,nns(dims(1)));
              sumranks = nansum(allpos,2);
              np = length(sumranks);
              for ii = length(dims):-1:1
                  plo = sum(sumranks<=w(dims(ii)))/np;
                  phi = sum(sumranks>=w(dims(ii)))/np;
                  p(dims(ii)) = min(plo,phi);
              end
          else
              p(dims) = NaN; % if no remaining sample
          end
      end
      p = min(2*p, 1);           % 2-sided, p>1 means the middle is double-counted
      stats.zval = []; %% TO BE COMPLETED!
      stats.R= [];%% TO BE COMPLETED!
else                          % use the normal approximation
   tiescor = 2 * tieadj / ((nnx+nny) .* (nnx+nny-1));
   wvar  = nnx.*nny.*((nnx + nny + 1) - tiescor)/12;
   wc = w - wmean;
   z = (wc - 0.5 .* sign(wc))./sqrt(wvar);
   p = 2*normcdf(-abs(z));
   stats.zval(nnx<=nny) = z(nnx<=nny);
   stats.zval(nnx>nny) = -z(nnx>nny);
   
   %JR addup: see http://yatani.jp/HCIstats/MannWhitney
   stats.R = stats.zval./sqrt(nnx+nny);
end

%% reshape to matrix
if length(sx)>2
    p = reshape(p,sx(2:end));
    stats.AUC = reshape(stats.AUC,sx(2:end));
    stats.R = reshape(stats.R,sx(2:end));
    stats.U = reshape(stats.U,sx(2:end));
    stats.zval = reshape(stats.zval,sx(2:end));
end

%% compute h
if nargout > 1,
   h = (p<=alpha);
   if length(sx)<=2
       stats.ranksum = w;
   else
       stats.ranksum = reshape(w,sx(2:end));
   end
end
