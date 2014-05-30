function myeb(X,Y,mycolor,varargin)
% myeb(X,Y,mycolor,varargin)

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colors to use:
% blue [.2 .3 .9]
% red [.9 .2 .2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% myeb(Y,varargin);
%
% This function makes nice coloured, shaded error bars. Exactly what
% it does depends on Y, and on whether you give it one or two inputs. 
%
% If you only pass it Y, and no other arguments, it assuemd you're
% giving it raw data. 
%
%		myeb(Raw_Data)
%
% 	.) if Y is 2D array, it will then plot mean(Y) with errorbars given
% 	by std(Y). In this case there is only one mean vector with its
% 	errorbars. 
% 
%	.) if Y is 3D array, it will plot size(Y,3) lines with the
%	associated errorbars. Line k will be mean(Y(:,:,k)) with errorbars
%	given by std(Y(:,:,k))
%
% If you pass it 2 arguments, each has to be at most 2D. 
%
%		myeb(mu,std)
%
% 	.) if mu and std are 1D, it just plots one line given by mu with a
% 	shaded region given by std. 
%
%	.) if mu and std are 2D, it will plot size(Y,2) lines in the
%	standard sequence of colours; each line mu(:,k) will have a shaded
%	region in the same colour, but less saturated given by std(:,k)
%
%
% Quentin Huys, 2007
% Center for Theoretical Neuroscience, Columbia University
% Email: qhuys [at] n e u r o theory [dot] columbia.edu
% (just get rid of the spaces, replace [at] with @ and [dot] with .)
if nargin < 3  
    mycolor = [.2 .3 .9];
end
if nargin <= 2
      Y = X;
end

%-- alpha
alpha = 0.05 ./ size(X,2); %bonferroni correction
try
load('normEstimators.mat'); % quantiles
quantile = find(normEstimators(:,1) > (1 - alpha),1);
if isempty(quantile) 
    quantile = normEstimators(end,1);
    display('Too small alpha value!')
end
catch
    quantile = 1;
    normEstimators = [1.96 2.576];
end

mycolor = mycolor ./ 4;

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;


if isempty(varargin);
	if length(size(Y))==2 
		m=nanmean(Y);
		s=nanstd(Y) / sqrt(size(Y,1)); 
		ind1=1:length(m);
		ind2=ind1(end:-1:1);
%         a = [1.96 2.576];%quantile for alpha = 0.05, and alpha = 0.01
		hold on; g=fill([X flipdim(X,2)],[m-(normEstimators(quantile,2) * s) m(ind2)+(normEstimators(quantile,2) * s(ind2))],mycolor .* 4);
%         hold on; h=fill([X flipdim(X,2)],[m-(normEstimators(quantile,2) * s) m(ind2)+(normEstimators(quantile,2) * s(ind2))],mycolor .* 3);        
		plot(X,m,'linewidth',2,'Color',mycolor .* 2);
%         set(h,'edgecolor',mycolor .* 3, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);  
        set(g,'edgecolor',mycolor .* 4, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5); 
        
		hold off
	elseif length(size(Y))>2 
		cla; hold on; 
		ind1=1:size(Y,2);
		ind2=ind1(end:-1:1);
		if size(Y,3)>8; col=jet(size(Y,3));ccol=col+.8; ccol(ccol>1)=1;end
		for k=1:size(Y,3)
			m=nanmean(Y(:,:,k));
			s=nanstd(Y(:,:,k));
			h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],ccol(k,:));
			set(h,'edgecolor',ccol(k,:))
		end
		for k=1:size(Y,3)
			m=nanmean(Y(:,:,k));
			s=nanstd(Y(:,:,k));
			plot(ind1,m,'linewidth',2,'color',col(k,:))
		end
		hold off 
	end

elseif length(varargin)==1;

	m=Y;
	s=varargin{1};
	if length(size(Y))>2; error;
	elseif min(size(Y))==1;
		if size(m,1)>1; m=m';s=s';end
		ind1=1:length(m);
		ind2=ind1(end:-1:1);
		hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],.6*ones(1,3));
		set(h,'edgecolor',.6*ones(1,3))
		plot(ind1,m,'linewidth',2)
		hold off
	else 
		ind1=(1:size(Y,1));
		ind2=ind1(end:-1:1);
		cla; hold on; 
		if size(Y,2)>8; col=jet(size(Y,2));ccol=col+.8; ccol(ccol>1)=1;end
		for k=1:size(Y,2)
			mm=m(:,k)';
			ss=s(:,k)';
			h=fill([ind1 ind2],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
			set(h,'edgecolor',ccol(k,:))
		end
		for k=1:size(Y,2);
			mm=m(:,k)';
			ss=s(:,k)';
			plot(ind1,mm,'linewidth',2,'color',col(k,:))
		end
		hold off 
	end
end


