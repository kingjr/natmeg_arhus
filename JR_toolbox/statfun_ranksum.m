function [s, cfg] = statfun_ranksum(cfg, dat, design)
% [s, cfg] = statfun_ranksum(cfg, dat, design)
% adaptation of STATFUN_indepsamplesT with a non parametric test.
% JR King

% STATFUN_indepsamplesT calculates the independent samples T-statistic 
% on the biological data in dat (the dependent variable), using the information on 
% the independent variable (iv) in design.
%
% Use this function by calling one of the high-level statistics functions as:
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'indepsamplesT'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = statfun_indepsamplesT(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variable (iv),  Nfac x Nreplications
%
% Configuration options:
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
%
% The following options are relevant if cfg.computecritval='yes' and/or
% cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail = -1, 0, or 1, left, two-sided, or right (default=1)
%              cfg.tail in combination with cfg.computecritval='yes'
%              determines whether the critical value is computed at
%              quantile cfg.alpha (with cfg.tail=-1), at quantiles
%              cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%              quantile (1-cfg.alpha) (with cfg.tail=1).
%
% Design specification:
%   cfg.ivar        = row number of the design that contains the labels of the conditions that must be 
%                        compared (default=1). The labels are the numbers 1 and 2.
%

% Copyright (C) 2006, Eric Maris
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the defaults
if ~isfield(cfg, 'computestat'),    cfg.computestat    = 'yes'; end
if ~isfield(cfg, 'computecritval'), cfg.computecritval = 'no';  end
if ~isfield(cfg, 'computeprob'),    cfg.computeprob    = 'no';  end
if ~isfield(cfg, 'alpha'),          cfg.alpha          = 0.05;  end
if ~isfield(cfg, 'tail'),           cfg.tail           = 1;     end

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
  % probabilities can only be calculated if the test statistics are calculated
  cfg.computestat = 'yes';
end;
if isfield(cfg,'uvar') && ~isempty(cfg.uvar)
    error('cfg.uvar should not exist for an independent samples statistic');
end

% perform some checks on the design
sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);
nreplc1 = sum(~isnan(dat(:,sel1)), 2);
nreplc2 = sum(~isnan(dat(:,sel2)), 2);
nrepl   = nreplc1 + nreplc2;
if any(nrepl<size(design,2)),
  warning_once('Not all replications are used for the computation of the statistic.');
end;

df = nrepl - 2;

s.df    = df;
[p h stats] = ranksum_matrix(dat(:,sel1)',dat(:,sel2)','alpha',cfg.alpha);
s.prob = p;
s.stat = [stats.ranksum]';
if strcmp(cfg.computecritval,'yes')
    for d = length(p):-1:1
        s.critval(d) = MW2cv(p(1),[length(sel1) length(sel2)]);
    end
end

function [p, h, stats] = ranksum_matrix(x,y,varargin)
% [p, h, stats] = ranksum_matrix(x,y,varargin)
% same as ranksum but allows mutldimensional arrays as input
% removes NaN values in each dimension
% enable matlabpool for faster computations
% (c) JR KING, jeanremi.king+matlab@gmail.com
if size(x,2) ~= size(y,2), error('x and y should be of the same dimensions!');end

original_size = size(x);
x = reshape(x,size(x,1),[]);
y = reshape(y,size(y,1),[]);

[p h] = deal(NaN(size(x,2),1));
dim = 1;
[p(dim), h(dim) stats(dim)] = inloop(x(:,dim),y(:,dim),varargin);
parfor dim = 2:size(x,2)
    if mod(dim,size(x,2)/20)==1, fprintf('*');end
    [p(dim), h(dim) stats(dim)] = inloop(x(:,dim),y(:,dim),varargin);
end
%-- reshape data
if length(original_size) > 2
    p = reshape(p,original_size(2:end));
    h = reshape(h,original_size(2:end));
end

function [p h stats] = inloop(xt,yt,varargin)
p = 1;h = 0; stats = [];
stats.zval = [];
stats.ranksum= [];

%-- deal with nan
xt(isnan(xt)) = [];
yt(isnan(yt)) = [];
if ~isempty(xt) && ~isempty(yt)
    if ~isempty(varargin)
        [p h stats] = main(xt,yt,varargin);
    else
        [p h stats] = main(xt,yt);
    end
end
return


function [p, h, stats] = main(x,y,varargin)
add = [];
try
    if ~isempty(varargin)
        for arg = 1:length(varargin)
            add = [add ', ' varargin{arg}];
        end
    end
    eval(['[p h stats] = ranksum(x,y' add ');']);
catch
    eval(['[p h stats] = ranksum(x,y);']);
end
return

function x = MW2cv(P,N)
%MW1CV   critical Mann-Whitney's U associated to a p-value. 
%It obtain a Mann-Whitney's U of two random variables with continuous cumulative
%distribution associated to a p-value. This procedure is highly recommended for sample sizes
%7< nx & ny <=40. For nx & ny <=7 it is recommended to use the MW1cv function;
%otherwise, U-value may be a poor approximation.
%It works with a procedure to get the nearest cumulative distribution relative value to P.
%[Based on the Fortran77 algorithm AS 62 Appl. Statist. (1973)]
%
%   Syntax: function x = MW2cv(P,N) 
%      
%     Inputs:
%          P - cumulative probability value of interest.
%          N - 2-element vector of sample sizes for the two samples []. 
% The input quantities should be scalars.
%     Outputs:
%          x - Mann-Whitney's U statistic.
%
%    Example: For two independent samples we are interested to get the
%             Mann-Whitney's statistic U with an associated cumulative
%             probability P = 0.95. Sample sizes are n1 = 36 and n2 = 14.
%
%                              P = 0.95; N = [36,14];
%
%     Calling on Matlab the function: 
%             x = MW2cv(P,N)
%
%       Answer is:
%
%                 328    
%

%  Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%
%  May 23, 2003.
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A. and R. Hernandez-Walls. (2003). MW2cv: Critical Mann-Whitney's U 
%    associated to a p-value: nx or ny >7. A MATLAB file. [WWW document]. URL http://
%    www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=3555&objectType=FILE
%
%  References:
% 
%  Mann, H. B. and Whitney, D. R. (1947), On a test of whether one of two   
%           random variables is stochastically larger than the other. Annals
%           of Mathematical Statistics, 18: 50-60.
%  Algorithm AS 62 (1973). Journal of Applied Statistics, 22(2):1-3.
%

if nargin <  2,
   error('Requires two input arguments.');
end

nmin = min(N);  %largest sample size.
nmax = max(N);  %smallest sample size.

if (nmin <= 7) & (nmax <= 7);
   fprintf('Warning: For nx and ny <= 7, the p-value may be a poor approximation.\n'); 
   fprintf('It is recommended to use the MW1cv function you can find on the\n');
   fprintf('Matlab>File Exchange Antonio Trujillo-Ortiz'' Author Page.\n');
   disp(' ');
   cont=input('Do you want to continue anyway (y/n):','s');
   
   if (cont=='y');
      disp('Here it goes.');
      
      mn1 = prod(N)+1;
      n1 = nmax+1;
      freq = [ones(n1,1); zeros(mn1-n1,1)];
      
      lwrk = floor((mn1+1)/2 + nmin);
      work = zeros(lwrk,1);
      
% Generate successively higher-order distributions
      in = nmax;
      for i = 2:nmin
         in = in+nmax;
         n1 = in+2;
         l = 1 + in/2;
         k = i;
         
% Generate complete distribution from outside inwards
         for j = 1:l
            k = k+1;
            n1 = n1-1;
            summ = freq(j) + work(j);
            freq(j) = summ;
            work(k) = summ - freq(n1);
            freq(n1) = summ;
         end;
      end;
      
      freq = freq/sum(freq);  % Make distribution relative
      
% Cumulative frequency distribution
      cumfreq = cumsum(freq);
      
%Location of the interested Mann-Whitney's U on all the possible U's for this 2-sample sizes.
%Here we are using a procedure to get the nearest fc value to P.
      cumfreq=cumfreq-P;
      u = find(abs(cumfreq)==min(abs(cumfreq(:))));
      
      UU = [0:length(freq)-1];  %vector of all the possible Mann-Whitney's U values.
      
%Association of the interested Mann-Whitney's U with its cumulative distribution.
      x = UU(u);
   else
   end
else
   
   mn1 = prod(N)+1;
   n1 = nmax+1;
   freq = [ones(n1,1); zeros(mn1-n1,1)];
   
   lwrk = floor((mn1+1)/2 + nmin);
   work = zeros(lwrk,1);
   
%Generate successively higher-order distributions
   in = nmax;
   for i = 2:nmin
      in = in+nmax;
      n1 = in+2;
      l = 1 + in/2;
      k = i;
      
%Generate complete distribution from outside inwards
      for j = 1:l
         k = k+1;
         n1 = n1-1;
         summ = freq(j) + work(j);
         freq(j) = summ;
         work(k) = summ - freq(n1);
         freq(n1) = summ;
      end;
   end;
   
   freq = freq/sum(freq);  % Make distribution relative
   
% Cumulative frequency distribution
   cumfreq = cumsum(freq);
   
%Location of the interested Mann-Whitney's U on all the possible U's for this 2-sample sizes.
%Here we are using a procedure to get the nearest fc value to P.
   cumfreq=cumfreq-P;
   u = find(abs(cumfreq)==min(abs(cumfreq(:))));
   
   UU = [0:length(freq)-1];  %vector of all the possible Mann-Whitney's U values.
   
%Association of the interested Mann-Whitney's U with its cumulative distribution.
   x = UU(u);      
end

